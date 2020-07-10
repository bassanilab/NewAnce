/*
Copyright (C) SIB - Swiss Institute of Bioinformatics, Lausanne, Switzerland
Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
*/

package newance.psmcombiner;

import newance.psmconverter.PeptideSpectrumMatch;
import newance.util.NewAnceParams;

import java.io.File;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * @author Markus Müller
 */
public class HistogramTree {

    private SmoothedScoreHistogram scoreHistogram;
    private final List<HistogramTree> children;
    private final HistogramTree parent;
    private final String id;
    private final int level;
    private final String group;

    public HistogramTree(SmoothedScoreHistogram scoreHistogram, String id, String group) {
        this.children = new ArrayList<>();
        this.scoreHistogram = scoreHistogram;
        this.parent = null;
        this.id = id;
        this.level = 0;
        this.group = group;
    }

    public HistogramTree(SmoothedScoreHistogram scoreHistogram, HistogramTree parent, String id, String group, int level) {
        this.children = new ArrayList<>();
        this.scoreHistogram = scoreHistogram;
        this.parent = parent;
        this.id = id;
        this.level = level;
        this.group = group;
    }

    public HistogramTree addHistogram(SmoothedScoreHistogram scoreHistogram, String id, String group) {

        HistogramTree histogramTree = new HistogramTree(scoreHistogram, this, id, group,level+1);
        children.add(histogramTree);

        return histogramTree;
    }

    public void add(PeptideSpectrumMatch peptideSpectrumMatch) {

        scoreHistogram.add(peptideSpectrumMatch);
    }

    public ScoreHistogram getScoreHistogram() {
        return scoreHistogram;
    }

    public void setScoreHistogram(SmoothedScoreHistogram scoreHistogram) {
        this.scoreHistogram = scoreHistogram;
    }

    public List<HistogramTree> getChildren() {
        return children;
    }

    public HistogramTree getParent() {
        return parent;
    }

    public String getId() {
        return id;
    }

    public boolean isLeaf() {
        return children.isEmpty();
    }

    public int getLevel() {
        return level;
    }

    public String getGroup() {
        return group;
    }

    public String print(float lFDRThreshold) {

        String treeString = "";
        if (level==0) treeString += String.format("lFDRThreshold = %.4f",lFDRThreshold)+"\n\n";
        String tab = "";
        for (int i=0;i<level;i++) tab += "\t";

        int tCnt = (int)scoreHistogram.getTotTargetCnt();
        int dCnt = (int)scoreHistogram.getTotDecoyCnt();
        float[] cnts = scoreHistogram.getTargetDecoyCounts(lFDRThreshold);
        String txt = String.format("totT=%d, totD=%d, selT=%d, selD=%d, grpFDR=%.4f, V=%d Pi_1=%.4f, Pi_0=%.4f",tCnt,dCnt,(int)cnts[1],(int)cnts[0],2*cnts[0]/(cnts[0]+cnts[1]),
                scoreHistogram.canCalculateFDR()?1:0,scoreHistogram.getPi_1(),scoreHistogram.getPi_0());

        treeString += tab+"{";
        if (!children.isEmpty()) {
            treeString += "\n"+tab+"\t"+id+": "+txt+"\n";
        } else {
            treeString += id+": "+txt;
        }

        for (HistogramTree node : children) {
            treeString += node.print(lFDRThreshold);
        }

        if (!children.isEmpty()) {
            treeString += tab+"}\n";
        } else {
            treeString += "}\n";
        }

        return treeString;
    }

    public String print(Map<String, Float> grpThresholdMap) {

        String treeString = "";
        if (level==0) for (String group : grpThresholdMap.keySet())  treeString += String.format("%s lFDRThreshold = %.4f",group, grpThresholdMap.get(group))+"\n";
        String tab = "";
        for (int i=0;i<level;i++) tab += "\t";

        int tCnt = (int)scoreHistogram.getTotTargetCnt();
        int dCnt = (int)scoreHistogram.getTotDecoyCnt();
        String txt = "";
        if (grpThresholdMap.containsKey(group)) {
            float[] cnts = scoreHistogram.getTargetDecoyCounts(grpThresholdMap.get(group));
            txt += String.format("totT=%d, totD=%d, selT=%d, selD=%d, grpFDR=%.4f, V=%d Pi_1=%.4f, Pi_0=%.4f", tCnt, dCnt, (int) cnts[1], (int) cnts[0], 2 * cnts[0] / (cnts[0] + cnts[1]),
                    scoreHistogram.canCalculateFDR() ? 1 : 0, scoreHistogram.getPi_1(), scoreHistogram.getPi_0());
        } else {
            txt += String.format("totT=%d, totD=%d, V=%d Pi_1=%.4f, Pi_0=%.4f", tCnt, dCnt,
                    scoreHistogram.canCalculateFDR() ? 1 : 0, scoreHistogram.getPi_1(), scoreHistogram.getPi_0());
        }

        treeString += tab+"{";
        if (!children.isEmpty()) {
            treeString += "\n"+tab+"\t"+id+": "+txt+"\n";
        } else {
            treeString += id+": "+txt;
        }

        for (HistogramTree node : children) {
            treeString += node.print(grpThresholdMap);
        }

        if (!children.isEmpty()) {
            treeString += tab+"}\n";
        } else {
            treeString += "}\n";
        }

        return treeString;
    }

    public void writeHistogram(String outputDir, String fileTag) {

        scoreHistogram.write(outputDir,fileTag,id);

        for (HistogramTree node : children) {
            node.writeHistogram(outputDir, fileTag);
        }
    }

    public void writeHistogram(String outputDir, String fileprefix, int level) {

        if (this.level<=level) {
            scoreHistogram.write(outputDir,fileprefix,id);
        }

        if (this.level==level) return;

        for (HistogramTree node : children) {
            node.writeHistogram(outputDir, fileprefix);
        }
    }

    public void setCanCalculateFDR(int minNrPsms) {

        scoreHistogram.setCanCalculateFDR(minNrPsms);

        for (HistogramTree node : children) {
            node.setCanCalculateFDR(minNrPsms);
        }
    }

    // import non-leaf histos
    public void importPriorCoreHistos() {

        if (isLeaf()) return;

        NewAnceParams params = NewAnceParams.getInstance();

        // import histos from files if provided
        if (!NewAnceParams.getInstance().getReadHistos().isEmpty() && (!scoreHistogram.canCalculateFDR() || params.isForceHistos())) {

            File histoFile = new File(params.getReadHistos() + File.separatorChar + "prior_histo_" + id + ".txt");

            if (!histoFile.exists()) {

                if (!id.equals("root") && !id.equals("Z1") && !id.equals("Z2") && !id.equals("Z3")) {
                    histoFile = new File(NewAnceParams.getInstance().getReadHistos() + File.separatorChar + "prior_histo_Z3.txt");
                    if (histoFile.exists()) {
                        System.out.println("WARNING: histogram file " + NewAnceParams.getInstance().getReadHistos() + File.separatorChar + "prior_histo_" + id + ".txt does not exist. Taking prior_histo_Z3.txt instead.");
                    } else {
                        System.out.println("ERROR: histogram file " + histoFile.getAbsolutePath() + " does not exist.  Abort.");
                        System.exit(1);
                    }
                } else {
                    System.out.println("ERROR: histogram file " + histoFile.getAbsolutePath() + " does not exist. Abort.");
                    System.exit(1);
                }
            }

            ScoreHistogram3D scoreHistogram3D = ScoreHistogram3D.read(histoFile);
            setScoreHistogram(scoreHistogram3D);
            scoreHistogram3D.setCanCalculateFDR(0);
        }

        for (HistogramTree node : children) {
            node.importPriorCoreHistos();
        }
    }


    public void calcClassProbs() {

        scoreHistogram.calcClassProb();

        for (HistogramTree node : children) {
            node.calcClassProbs();
        }
    }

    public void smoothHistogram(int degree) {

        if (degree<=0) return;

        for (int i=0;i<degree-1;i++) scoreHistogram.smoothHistogram(false);
        scoreHistogram.smoothHistogram(true);

        for (HistogramTree node : children) {
            node.smoothHistogram(degree);
        }
    }


    public void calcLocalFDR() {

        if (scoreHistogram.canCalculateFDR())
            scoreHistogram.calcLocalFDR();
        else {
            float pi1 = (float)scoreHistogram.getPi_1();
            float pi0 = (float)scoreHistogram.getPi_0();

            HistogramTree validParent = getNextValidParent();
            float ratio = (pi0==0)?10000:pi1/pi0;
            scoreHistogram.calcLocalFDR(ratio, validParent.getScoreHistogram());

            if (validParent.getScoreHistogram().canCalculateFDR()) {
                System.out.println("Not enough PSMs in group "+id+" to calculate local FDR. Parent histogram "+parent.getId()+" used instead.");
            } else {
                System.out.println("Not enough PSMs to calculate local FDR in parent histogram: "+validParent.getId()+". Results may not be reliable. Try importing prior histos.");
            }
        }

        for (HistogramTree node : children) {
            node.calcLocalFDR();
        }
    }

    public float[] getTargetDecoyCounts(float maxlFDR) {

        return scoreHistogram.getTargetDecoyCounts(maxlFDR);
    }

    private HistogramTree getNextValidParent() {

        HistogramTree p = parent;
        HistogramTree validP = (parent==null)?this:parent;

        while (p!=null && !p.scoreHistogram.canCalculateFDR()) {
            validP = p;
            p = p.getParent();
        }

        return validP;
    }

    public void clear() {

        scoreHistogram.clear();

        for (HistogramTree node : children) {
            node.clear();
        }
    }



}
