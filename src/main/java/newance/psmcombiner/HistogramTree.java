/**
 * Copyright (C) 2019, SIB/LICR. All rights reserved
 *
 * SIB, Swiss Institute of Bioinformatics
 * Ludwig Institute for Cancer Research (LICR)
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer. Redistributions in binary
 * form must reproduce the above copyright notice, this list of conditions and
 * the following disclaimer in the documentation and/or other materials provided
 * with the distribution. Neither the name of the SIB/LICR nor the names of
 * its contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL SIB/LICR BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

package newance.psmcombiner;

import newance.psmconverter.PeptideSpectrumMatch;
import newance.util.NewAnceParams;

import java.io.File;
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

    public void writeHistogram(String outputDir, String fileprefix) {

        String filename = outputDir+File.separatorChar+fileprefix+"_"+id+".txt";
        scoreHistogram.write(new File(filename));

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

    public void importPriorHisto() {

        if (!NewAnceParams.getInstance().getReadHistos().isEmpty()) {
            File histoFile = new File(NewAnceParams.getInstance().getReadHistos() + File.separatorChar + "prior_histo_" + id + ".txt");

            if (!histoFile.exists()) {
                System.out.println("ERROR: histogram file "+histoFile.getAbsolutePath()+" for charge "+id+" does not exist. Please check your -minZ, -maxZ options. Abort.");
                System.exit(1);
            }
            CometScoreHistogram scoreHistogram = CometScoreHistogram.read(histoFile);
            setScoreHistogram(scoreHistogram);
            scoreHistogram.setCanCalculateFDR(NewAnceParams.getInstance().getMinNrPsmsPerHisto());
        } else if (!scoreHistogram.canCalculateFDR() && !isLeaf()) {

            System.out.println("Setting histogram " +id+" to prior values since there are not enough PSMs to estimate the distribution.");
            File histoFile = new File(getClass().getResource("prior_histo_" + id + ".txt").getFile());

            if (!histoFile.exists()) {
                System.out.println("ERROR: histogram file "+histoFile.getAbsolutePath()+" for charge "+id+" does not exist. Please check your -minZ, -maxZ options. Abort.");
                System.exit(1);
            }

            CometScoreHistogram scoreHistogram = CometScoreHistogram.read(histoFile);
            setScoreHistogram(scoreHistogram);
            scoreHistogram.setCanCalculateFDR(0);
        }

        for (HistogramTree node : children) {
            node.importPriorHisto();
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

            if (validParent!=null) {
                float ratio = (pi0==0)?10000:pi1/pi0;
                scoreHistogram.calcLocalFDR(ratio, parent.getScoreHistogram());
            } else {
                System.out.println("Cannot calculate local FDR of group: "+id);
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
        while (!p.scoreHistogram.canCalculateFDR()) p = p.getParent();

        return p;
    }

}
