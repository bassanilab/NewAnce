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

import newance.proteinmatch.UniProtDB;
import newance.psmconverter.AddUniProtIds2Psm;
import newance.psmconverter.PeptideSpectrumMatch;
import newance.util.NewAnceParams;
import newance.util.PsmGrouper;

import java.io.File;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

/**
 * @author Markus Müller
 */

public class GroupedFDRCalculator {

    protected final HistogramTree histogramTreeRoot;
    protected final Map<String, HistogramTree> histogramMap;
    protected final PsmGrouper psmGrouper;
    protected final int[] nrBins;
    protected final AddUniProtIds2Psm addUniProtIds2Psm;

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
        this.addUniProtIds2Psm = null;

        buildTree();
    }

    public GroupedFDRCalculator(PsmGrouper psmGrouper, UniProtDB uniProtDB) {

        this.nrBins = new int[3];
        nrBins[0] = NewAnceParams.getInstance().getNrXCorrBins();
        nrBins[1] = NewAnceParams.getInstance().getNrDeltaCnBins();
        nrBins[2] = NewAnceParams.getInstance().getNrSpScoreBins();

        CometScoreHistogram cometScoreHistogram = new CometScoreHistogram(nrBins);

        this.histogramTreeRoot = new HistogramTree(cometScoreHistogram, "root","");

        this.histogramMap = new HashMap<>();
        this.histogramMap.put("root",this.histogramTreeRoot);
        this.psmGrouper = psmGrouper;
        if (uniProtDB!=null)
            this.addUniProtIds2Psm = new AddUniProtIds2Psm(uniProtDB);
        else
            this.addUniProtIds2Psm = null;

        buildTree();
    }

    public GroupedFDRCalculator() {

        this.nrBins = new int[3];
        nrBins[0] = NewAnceParams.getInstance().getNrXCorrBins();
        nrBins[1] = NewAnceParams.getInstance().getNrDeltaCnBins();
        nrBins[2] = NewAnceParams.getInstance().getNrSpScoreBins();

        CometScoreHistogram cometScoreHistogram = new CometScoreHistogram(nrBins);

        this.histogramTreeRoot = new HistogramTree(cometScoreHistogram, "root","");

        this.histogramMap = new HashMap<>();
        this.histogramMap.put("root",this.histogramTreeRoot);
        this.psmGrouper = null;
        this.addUniProtIds2Psm = null;

        buildTree();
    }


    private void buildTree() {

        Set<String> groups = new HashSet<>();

        if (psmGrouper!=null) groups = psmGrouper.getGroups();
        CometScoreHistogram cometScoreHistogram;

        for (int Z=NewAnceParams.getInstance().getMinCharge();Z<=NewAnceParams.getInstance().getMaxCharge();Z++) {
            cometScoreHistogram = new CometScoreHistogram(nrBins);
            String label = "Z"+String.valueOf(Z);
            HistogramTree node = histogramTreeRoot.addHistogram(cometScoreHistogram,label,"");
            this.histogramMap.put(label,node);

            for (String group : groups) {
                cometScoreHistogram = new CometScoreHistogram(nrBins);
                label = "Z"+String.valueOf(Z)+"_"+group;
                HistogramTree leaf = node.addHistogram(cometScoreHistogram,label, group);
                this.histogramMap.put(label,leaf);
            }
        }
    }


    public void readPriorHistograms(String outputDir, String fileprefix) {

        for (int Z=NewAnceParams.getInstance().getMinCharge();Z<=NewAnceParams.getInstance().getMaxCharge();Z++) {
            String label = "Z"+String.valueOf(Z);

            HistogramTree node = histogramMap.get(label);

            if (!node.getScoreHistogram().canCalculateFDR()) {

                String filename = outputDir+File.separatorChar+fileprefix+"_"+label+".txt";
                node.setScoreHistogram(CometScoreHistogram.read(new File(filename)));

            }
        }
    }


    public void addAll(ConcurrentHashMap<String,List<PeptideSpectrumMatch>> psms) {

        for (String specID : psms.keySet()) {
            if (addUniProtIds2Psm!=null) addUniProtIds2Psm.accept(specID, psms.get(specID));

            for (PeptideSpectrumMatch psm : psms.get(specID)) {
                add(psm);
            }
        }
    }


    protected void add(PeptideSpectrumMatch peptideSpectrumMatch) {

        String id = getNodeID(peptideSpectrumMatch);

        HistogramTree histogramTree = histogramMap.get(id);

        while (histogramTree != null) {
            histogramTree.add(peptideSpectrumMatch);
            histogramTree = histogramTree.getParent();
        }
    }


    public String getNodeID(PeptideSpectrumMatch peptideSpectrumMatch) {

        if (psmGrouper==null)
            return "Z"+ peptideSpectrumMatch.getCharge();
        else
            return "Z"+ peptideSpectrumMatch.getCharge()+"_"+psmGrouper.apply("", peptideSpectrumMatch);
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

    public void writeHistograms(String outputDir, String filePrefix, int level) {
        try {
            File dir = new File(outputDir);
            if (!dir.exists()) dir.mkdirs();
        } catch (SecurityException e) {
            System.out.println("Cannot create directory "+outputDir);
        }

        histogramTreeRoot.writeHistogram(outputDir, filePrefix, level);
    }

    public void setCanCalculateFDR(int minNrPsms) {

        histogramTreeRoot.setCanCalculateFDR(minNrPsms);
    }

    public void importPriorHistos() {

        histogramTreeRoot.importPriorCoreHistos();
    }


    public void calcLocalFDR() {

        for (HistogramTree node : histogramMap.values()) {

//            if (node.isLeaf()) {
                node.calcLocalFDR();
//            }
        }
    }

    public float getLocalFDR(PeptideSpectrumMatch peptideSpectrumMatch) {

        String id = getNodeID(peptideSpectrumMatch);

        HistogramTree histogramTree = histogramMap.get(id);

        return histogramTree.getScoreHistogram().getLocalFDR(peptideSpectrumMatch);
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

        for (HistogramTree node : histogramMap.values()) {

            if (node.isLeaf() && (group.isEmpty() || node.getGroup().equals(group))) {
                float[] values = node.getTargetDecoyCounts(lFDR);

                System.out.println(node.getId()+","+node.getGroup()+","+node.isLeaf()+". Decoy: "+values[0]+". Target: "+values[1]);
                decoySum += values[0];
                targetSum += values[1];
            }
        }

        return new float[] {decoySum, targetSum};
    }

    public float calcGlobalFDR(float lFDR, String group) {

        float decoySum = 0;
        float targetSum = 0;

        for (HistogramTree node : histogramMap.values()) {

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

        // find lfdr interval where pfdr becomes larger than pFDR
        for (int i=0;i<=10;i++) {
            lfdr = i*0.1f;
            pfdr = calcGlobalFDR(lfdr, group);
            if (pfdr>=pFDR) break;

            pfdrP = pfdr;
            lfdrP = lfdr;
        }

        float eps = Math.abs(pfdr-pFDR);
        if (eps<0.0001) return lfdr;

        if (pfdr<pFDR) return 1f;

        // refine search
        return calcLocalFDRThreshold(lfdrP, pfdrP, lfdr, pfdr, pFDR, eps, group);
    }

    private float calcLocalFDRThreshold(float lFDRLeft, float pFDRLeft, float lFDRRight, float pFDRRight, float pFDR, float eps, String group) {

        float f = (pFDR-pFDRLeft)/(pFDRRight-pFDRLeft);
        float lfdr = lFDRLeft + f*(lFDRRight-lFDRLeft);

        float pfdr = calcGlobalFDR(lfdr, group);
        float epsN = Math.abs(pfdr-pFDR);

        if (epsN<0.0001 || epsN >= eps) return lfdr;

        float lfdrN = 0;
        if (pfdr<pFDR)
            lfdrN = calcLocalFDRThreshold(lfdr, pfdr, lFDRRight, pFDRRight, pFDR, epsN, group);
        else
            lfdrN = calcLocalFDRThreshold(lFDRLeft, pFDRLeft, lfdr, pfdr, pFDR, epsN, group);

        return lfdrN;
    }

    public ConcurrentHashMap<String, List<PeptideSpectrumMatch>> filterPsms(ConcurrentHashMap<String, List<PeptideSpectrumMatch>>  psms, float lFDRThreshold, String group) {

        ConcurrentHashMap<String, List<PeptideSpectrumMatch>> filteredPsms = new ConcurrentHashMap<>();

        // do parallel later
        for (String specID : psms.keySet()) {
            for (PeptideSpectrumMatch psm : psms.get(specID)) {

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


    public ConcurrentHashMap<String, List<PeptideSpectrumMatch>> filterPsms(ConcurrentHashMap<String, List<PeptideSpectrumMatch>> psms, float lFDRThreshold) {

        ConcurrentHashMap<String, List<PeptideSpectrumMatch>> filteredPsms = new ConcurrentHashMap<>();

        // do parallel later
        for (String specID : psms.keySet()) {
            for (PeptideSpectrumMatch psm : psms.get(specID)) {

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
        if (psmGrouper!=null) return psmGrouper.getGroups();
        else return new HashSet<>();
    }

    public Map<String,Float> calcGroupLocalFDRThreshold(float pFDR) {

        Map<String,Float> grplFDRMap = new HashMap<>();

        for (String group : getGroups()) {
            grplFDRMap.put(group, calcLocalFDRThreshold(pFDR, group));
        }

        return grplFDRMap;
    }

}
