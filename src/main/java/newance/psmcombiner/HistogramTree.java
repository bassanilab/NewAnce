package newance.psmcombiner;

import newance.psmconverter.PeptideMatchData;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * @author Markus Muller
 */
public class HistogramTree {

    private final SmoothedScoreHistogram scoreHistogram;
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

    public void add(PeptideMatchData peptideMatchData) {

        scoreHistogram.add(peptideMatchData);
    }

    public ScoreHistogram getScoreHistogram() {
        return scoreHistogram;
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
        String txt = String.format("totT=%d, totD=%d, selT=%d, selD=%d, grpFDR=%.4f, V=%d Pi_1= %.4f, Pi_0= %.4f",tCnt,dCnt,(int)cnts[1],(int)cnts[0],2*cnts[0]/(cnts[0]+cnts[1]),
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
        float[] cnts = scoreHistogram.getTargetDecoyCounts(grpThresholdMap.get(group));
        String txt = String.format("totT=%d, totD=%d, selT=%d, selD=%d, grpFDR=%.4f, V=%d Pi_1= %.4f, Pi_0= %.4f",tCnt,dCnt,(int)cnts[1],(int)cnts[0],2*cnts[0]/(cnts[0]+cnts[1]),
                scoreHistogram.canCalculateFDR()?1:0,scoreHistogram.getPi_1(),scoreHistogram.getPi_0());

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


    public void calcClassProbs() {

        scoreHistogram.calcClassProb();

        for (HistogramTree node : children) {
            node.calcClassProbs();
        }
    }

    public void smoothHistogram(int degree) {

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

            if (pi0>0 && validParent!=null) {
                scoreHistogram.calcLocalFDR(pi1/pi0, parent.getScoreHistogram());
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
