package newance.psmcombiner;

import newance.psmconverter.PeptideMatchData;
import newance.util.NewAnceParams;
import newance.util.PsmGrouper;

import java.io.File;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Created by Markus Muller on 18.11.19.
 */
public class GroupedFDRCalculator {

    protected final HistogramTree histogramTreeRoot;
    protected final Map<String, HistogramTree> histogramMap;
    protected final PsmGrouper psmGrouper;
    protected final int[] nrBins;

    public GroupedFDRCalculator(PsmGrouper psmGrouper) {

        this.nrBins = new int[3];
        nrBins[0] = NewAnceParams.getInstance().getNrXCorrBins();
        nrBins[1] = NewAnceParams.getInstance().getNrDeltaCnBins();
        nrBins[2] = NewAnceParams.getInstance().getNrSpScoreBins();

        CometScoreHistogram cometScoreHistogram = new CometScoreHistogram(nrBins);

        this.histogramTreeRoot = new HistogramTree(cometScoreHistogram, "root","");

        this.histogramMap = new HashMap<>();
        this.histogramMap.put("root",this.histogramTreeRoot);
        this.psmGrouper = psmGrouper;

        buildTree();
    }

    private void buildTree() {

        Set<String> groups = psmGrouper.getGroups();
        CometScoreHistogram cometScoreHistogram;

        for (int Z=NewAnceParams.getInstance().getMinCharge();Z<=NewAnceParams.getInstance().getMaxCharge();Z++) {
            cometScoreHistogram = new CometScoreHistogram(nrBins);
            HistogramTree node = histogramTreeRoot.addHistogram(cometScoreHistogram,"Z"+String.valueOf(Z),"");
            String label = "Z"+String.valueOf(Z);
            this.histogramMap.put(label,node);

            for (String group : groups) {
                cometScoreHistogram = new CometScoreHistogram(nrBins);
                label = "Z"+String.valueOf(Z)+"_"+group;
                HistogramTree leaf = node.addHistogram(cometScoreHistogram,label, group);
                this.histogramMap.put(label,leaf);
            }
        }
    }

    public void addAll(ConcurrentHashMap<String,List<PeptideMatchData>> psms) {

        for (String specID : psms.keySet()) {
            for (PeptideMatchData psm : psms.get(specID)) {
                add(psm);
            }
        }
    }

    public void add(PeptideMatchData peptideMatchData) {
        String id = getNodeID(peptideMatchData);

        HistogramTree histogramTree = histogramMap.get(id);

        while (histogramTree != null) {
            histogramTree.add(peptideMatchData);
            histogramTree = histogramTree.getParent();
        }
    }

    public String getNodeID(PeptideMatchData peptideMatchData) {

        return "Z"+peptideMatchData.getCharge()+"_"+psmGrouper.apply("",peptideMatchData);
    }

    public String printTree(float lFDRThreshold) {

        return histogramTreeRoot.print(lFDRThreshold);
    }

    public String printTree(Map<String, Float> grpThresholdMap) {

        return histogramTreeRoot.print(grpThresholdMap);
    }


    public void smoothHistogram(int degree) {

        histogramTreeRoot.smoothHistogram(degree);
    }

    public void writeHistograms(String outputDir, String filePrefix) {
        try {
            File dir = new File(outputDir);
            if (!dir.exists()) dir.mkdirs();
        } catch (SecurityException e) {
            System.out.println("Cannot create directory "+outputDir);
        }

        histogramTreeRoot.writeHistogram(outputDir, filePrefix);
    }

    public void setCanCalculateFDR(int minNrPsms) {

        histogramTreeRoot.setCanCalculateFDR(minNrPsms);
    }


    public void calcLocalFDR() {

        for (String id : histogramMap.keySet()) {

            HistogramTree node = histogramMap.get(id);
            if (node.isLeaf()) {
                node.calcLocalFDR();
            }
        }
    }

    public void calcClassProbs() {
        histogramTreeRoot.calcClassProbs();
    }

    public float[] getTargetDecoyCounts(float lFDR) {

        return getTargetDecoyCounts(lFDR, "");
    }

    public float[] getTargetDecoyCounts(float lFDR, String group) {

        float decoySum = 0;
        float targetSum = 0;

        for (String id : histogramMap.keySet()) {

            HistogramTree node = histogramMap.get(id);
            if (node.isLeaf() && (group.isEmpty()) || node.getGroup().equals(group)) {
                float[] values = node.getTargetDecoyCounts(lFDR);
                decoySum += values[0];
                targetSum += values[1];
            }
        }

        return new float[] {decoySum, targetSum};
    }

    public float calcGlobalFDR(float lFDR, String group) {

        float decoySum = 0;
        float targetSum = 0;

        for (String id : histogramMap.keySet()) {

            HistogramTree node = histogramMap.get(id);
            if (node.isLeaf() && (group.isEmpty()) || node.getGroup().equals(group)) {
                float[] values = node.getTargetDecoyCounts(lFDR);
                decoySum += values[0];
                targetSum += values[1];
            }
        }

        return (decoySum+targetSum>0)?2*decoySum/(decoySum+targetSum):0f;
    }

    public float calcGlobalFDR(float lFDR) {

        return calcGlobalFDR(lFDR, "");
    }

    public float calcLocalFDRThreshold(float pFDR) {
       return calcLocalFDRThreshold(pFDR, "");
    }

    public float calcLocalFDRThreshold(float pFDR, String group) {

        if (pFDR==0) return 0;

        float pfdr = 0;
        float lfdr = 0;
        float pfdrP = 0;
        float lfdrP = 0;

        for (int i=0;i<10;i++) {
            lfdr = i*0.1f;
            pfdr = calcGlobalFDR(lfdr, group);
            if (pfdr>=pFDR) break;

            pfdrP = pfdr;
            lfdrP = lfdr;
        }

        float eps = Math.abs(pfdr-pFDR);
        if (eps<0.0001) return lfdr;


        return calcLocalFDRThreshold(lfdrP, pfdrP, lfdr, pfdr, pFDR, eps, group);
    }

    private float calcLocalFDRThreshold(float lFDRLeft, float pFDRLeft, float lFDRRight, float pFDRRight, float pFDR, float eps, String group) {

        float f = (pFDR-pFDRLeft)/(pFDRRight-pFDRLeft);
        float lfdr = lFDRLeft + f*(lFDRRight-lFDRLeft);

        float pfdr = calcGlobalFDR(lfdr, group);
        float epsN = Math.abs(pfdr-pFDR);

        if (epsN<0.0001 || epsN >= eps) return lfdr;

        float lfdrN = 0;
        if (pfdr>pFDR)
            lfdrN = calcLocalFDRThreshold(lfdr, pfdr, lFDRRight, pFDRRight, pFDR, epsN, group);
        else
            lfdrN = calcLocalFDRThreshold(lFDRLeft, pFDRLeft, lfdr, pfdr, pFDR, epsN, group);

        return lfdrN;
    }

    public ConcurrentHashMap<String, List<PeptideMatchData>> filterPsms(ConcurrentHashMap<String, List<PeptideMatchData>>  psms, float lFDRThreshold, String group) {

        ConcurrentHashMap<String, List<PeptideMatchData>> filteredPsms = new ConcurrentHashMap<>();

        // do parallel later
        for (String specID : psms.keySet()) {
            for (PeptideMatchData psm : psms.get(specID)) {

                String grp = psmGrouper.apply(specID,psm);

                if (!grp.equals(group)) continue;

                String label = getNodeID(psm);
                HistogramTree node = histogramMap.get(label);
                float lFDR = node.getScoreHistogram().getLocalFDR(psm);

                if (lFDR <= lFDRThreshold) {
                    filteredPsms.putIfAbsent(specID, Collections.synchronizedList(new ArrayList<>()));
                    filteredPsms.get(specID).add(psm);
                }
            }
        }

        return filteredPsms;
    }


    public ConcurrentHashMap<String, List<PeptideMatchData>> filterPsms(ConcurrentHashMap<String, List<PeptideMatchData>>  psms, float lFDRThreshold) {

        ConcurrentHashMap<String, List<PeptideMatchData>> filteredPsms = new ConcurrentHashMap<>();

        // do parallel later
        for (String specID : psms.keySet()) {
            for (PeptideMatchData psm : psms.get(specID)) {

                String label = getNodeID(psm);
                HistogramTree node = histogramMap.get(label);
                float lFDR = node.getScoreHistogram().getLocalFDR(psm);

                if (lFDR <= lFDRThreshold) {
                    filteredPsms.putIfAbsent(specID, Collections.synchronizedList(new ArrayList<>()));
                    filteredPsms.get(specID).add(psm);
                }
            }
        }

        return filteredPsms;
    }

    public Set<String> getGroups() {
        return psmGrouper.getGroups();
    }

    public Map<String,Float> calcGroupLocalFDRThreshold(float pFDR) {

        Map<String,Float> grplFDRMap = new HashMap<>();

        for (String group : getGroups()) {
            grplFDRMap.put(group, calcLocalFDRThreshold(pFDR, group));
        }

        return grplFDRMap;
    }

}
