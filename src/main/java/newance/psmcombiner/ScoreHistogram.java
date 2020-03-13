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

import java.io.File;
import java.io.Serializable;
import java.util.*;

/**
 * @author Markus Müller
 */

public abstract class ScoreHistogram implements Serializable {

    protected final List<Float> gamma;
    protected final List<Float> lFDR;
    protected List<Integer> sortedIndexes;
    protected final List<Float> pFDR;
    protected List<Float> targetCnts;
    protected List<Float> decoyCnts;
    protected final List<Integer> indexMap;
    protected List<Integer> psmBins;
    protected float totTargetCnt;
    protected float totDecoyCnt;
    protected final int totNrBins;
    protected final int[] nrBins;
    protected final int dimension;
    protected int currIndex;
    protected double pi_0;
    protected double pi_1;
    protected boolean canCalculateFDR;

    public ScoreHistogram(int[] nrBins) {

        this.totTargetCnt = 0;
        this.totDecoyCnt = 0;
        this.nrBins = nrBins;
        this.dimension = nrBins.length;
        int cnt = 1;
        for (int c : nrBins) cnt *= c;
        totNrBins = cnt;

        targetCnts = new ArrayList<>();
        decoyCnts = new ArrayList<>();
        psmBins = new ArrayList<>();
        gamma = new ArrayList<>();
        lFDR = new ArrayList<>();
        pFDR = new ArrayList<>();
        sortedIndexes = null;
        indexMap = new ArrayList<>(totNrBins);
        for (int i=0;i<totNrBins;i++) indexMap.add(-1);
        currIndex = 0;
        canCalculateFDR = false;
        pi_0 = -1.0;
        pi_1 = -1.0;
    }

    public ScoreHistogram(ScoreHistogram scoreHistogram) {

        this.totTargetCnt = scoreHistogram.totTargetCnt;
        this.totDecoyCnt = scoreHistogram.totDecoyCnt;
        this.nrBins = scoreHistogram.nrBins;
        this.dimension = scoreHistogram.dimension;
        totNrBins = scoreHistogram.totNrBins;
        targetCnts = new ArrayList<>(scoreHistogram.targetCnts);
        decoyCnts = new ArrayList<>(scoreHistogram.decoyCnts);
        psmBins = new ArrayList<>(scoreHistogram.psmBins);
        gamma = new ArrayList<>(scoreHistogram.gamma);
        lFDR = new ArrayList<>(scoreHistogram.lFDR);
        pFDR = new ArrayList<>(scoreHistogram.pFDR);
        if (scoreHistogram.sortedIndexes!=null)
            sortedIndexes = new ArrayList<>(scoreHistogram.sortedIndexes);
        else
            sortedIndexes = null;

        indexMap = new ArrayList<>(scoreHistogram.indexMap);
        currIndex = scoreHistogram.currIndex;
        canCalculateFDR = scoreHistogram.canCalculateFDR;
        pi_0 = scoreHistogram.pi_0;
        pi_1 = scoreHistogram.pi_1;
    }

    public void add(List<PeptideSpectrumMatch> peptideSpectrumMatchList) {

        for (PeptideSpectrumMatch peptideSpectrumMatch : peptideSpectrumMatchList) {
            add(peptideSpectrumMatch);
        }
    }

    public void add(PeptideSpectrumMatch peptideSpectrumMatch) {

        int bin = index(peptideSpectrumMatch);

        int idx = indexMap.get(bin);
        if (idx<0) {
            indexMap.set(bin,currIndex);
            psmBins.add(bin);
            currIndex++;
        }

        float isDecoy =  peptideSpectrumMatch.isDecoy()?1:0;

        if (idx < 0) { // new bin
            targetCnts.add(1-isDecoy);
            decoyCnts.add(isDecoy);
        } else { // already recorded bin
            if (isDecoy==0) targetCnts.set(idx,targetCnts.get(idx)+1);
            else decoyCnts.set(idx,decoyCnts.get(idx)+1);
        }

        if (isDecoy==0) totTargetCnt++;
        else totDecoyCnt++;
    }

    protected void calcClassProb() {

        List<Float> pvalues = calcPValues();

        int n = pvalues.size();

        if (n==0) {
            pi_0 = pi_1 = 0.5;
            return;
        }

        float th = 0.05f/n;
        int i;
        for (i=0;i<n;i++) {
            if (pvalues.get(i)>th) break;
        }

        pi_1 = 1.05f*i/n;
        if (pi_1>=1.0) {
            pi_1 = totTargetCnt/(totTargetCnt+totDecoyCnt);
        }

//        if (totTargetCnt>0) {
//            pi_0 = (2.0 * totDecoyCnt) / (totTargetCnt+totDecoyCnt);
//            pi_0 = (pi_0>1.0)?1.0:pi_0;
//        } else if (totDecoyCnt>0){
//            pi_0 = 1.0;
//        }  else {
//            pi_0 = 0.5;
//        }

        pi_0 = 1.0 - pi_1;
    }

    protected List<Float> calcPValues() {

        List<Float> d = new ArrayList<>();
        List<Float> t = new ArrayList<>();

        for (Float cnt: decoyCnts) {
            t.add(cnt);
            if (cnt>0) d.add(cnt);
        }

        List<Integer> sortedIdxD = sortIndexes(d,false);
        List<Integer> sortedIdxT = sortIndexes(t,false);
        List<Float> pvalues = new ArrayList<>();

        int n = sortedIdxD.size();
        int m = sortedIdxT.size();

        for (int i=0,j=0;j<n;j++) {

            while (i<m && t.get(sortedIdxT.get(i))<d.get(sortedIdxD.get(j))) {
                i++;
                pvalues.add(1.0f*j/n);
            }
        }

        return pvalues;
    }

    protected void calcGamma() {

        gamma.clear();
        float fact = totDecoyCnt/totTargetCnt;
        for (int i=0;i<targetCnts.size();i++) {

            float ratio = -1;
            if (decoyCnts.get(i) > 0) ratio = targetCnts.get(i)/decoyCnts.get(i);

            gamma.add(ratio*fact);
        }
    }

    public float getGamma(int bin) {

        if (gamma.isEmpty()) calcGamma();

        int idx = indexMap.get(bin);

        if (idx>=0)
            return gamma.get(idx);
        else
            return -1;
    }

    public float getLocalFDR(PeptideSpectrumMatch psm) {

        if (lFDR.isEmpty()) calcLocalFDR();

        int bin = index(psm);

        return getLocalFDR(bin);
    }

    public float getLocalFDR(int bin) {

        if (lFDR.isEmpty()) calcLocalFDR();

        int idx = indexMap.get(bin);

        if (idx>=0)
            return lFDR.get(idx);
        else
            return -1;
    }

    public void calcLocalFDR() {

        if (pi_0<0 || pi_1<0) calcClassProb();

        float p1_p0_ratio = (float) (pi_1/pi_0);
        if (gamma.isEmpty()) calcGamma();

        lFDR.clear();
        for (int i=0;i<targetCnts.size();i++) {

            float lfdr = 0;
            if (gamma.get(i)>=0) {
                lfdr = 1f/(1+p1_p0_ratio*gamma.get(i));
            }

            lFDR.add((lfdr>1)?1f:lfdr);
        }

        sortedIndexes = sortIndexes(lFDR,false);
    }

    public void calcLocalFDR(float p1_p0_ratio, ScoreHistogram parent) {

        if (p1_p0_ratio<0 || parent==null) return;

        gamma.clear();
        lFDR.clear();
        for (Integer bin : psmBins) {

            float g = parent.getGamma(bin);

            gamma.add(g);
            float lfdr = 0;
            if (g>=0) {
                lfdr = 1f/(1+p1_p0_ratio*g);
            }
            lFDR.add((lfdr>1)?1f:lfdr);
        }

        sortedIndexes = sortIndexes(lFDR,false);
    }

    public float[] getTargetDecoyCounts(float maxlFDR) {

        if (lFDR.isEmpty()) calcLocalFDR();

        float decoySum = 0;
        float targetSum = 0;
        float decoySumP = 0;
        float targetSumP = 0;
        float fdrP = 0;

        for (int i = 0; i < sortedIndexes.size(); i++) {
            int idx = sortedIndexes.get(i);

            float fdr = lFDR.get(idx);
            if (fdr>maxlFDR) {
                float f = (fdr-maxlFDR)/(fdr-fdrP);
                decoySum = decoySum + f*(decoySumP-decoySum);
                targetSum = targetSum + f*(targetSumP-targetSum);
                break;
            }

            decoySum += decoyCnts.get(idx);
            targetSum += targetCnts.get(idx);

            fdrP = fdr;
            decoySumP = decoySum;
            targetSumP = targetSum;
        }

        return new float[] {decoySum, targetSum};
    }

    protected static List<Integer> sortIndexes(List<Float> array, boolean reverse) {

        List<Integer> sortedIndexes = new ArrayList<>();

        for (int i=0; i < array.size(); i++) {
            sortedIndexes.add(i);
        }

        Comparator<Integer> comparator = new Comparator<Integer>() {
            @Override
            public int compare(Integer i1, Integer i2) {
                if (reverse)
                    return array.get(i2).compareTo(array.get(i1));
                else
                    return array.get(i1).compareTo(array.get(i2));
            }
        };

        sortedIndexes.sort(comparator);

        return sortedIndexes;
    }

    public List<Float> calcMids(List<Float> breaks) {
        List<Float> mids = new ArrayList<>();

        for (int i=0;i<breaks.size()-1;i++) {
            mids.add((breaks.get(i)+breaks.get(i+1))/2);
        }

        return mids;
    }


    public List<Float> calcBreaks(float minScore, float maxScore, int nrIntervals) {
        List<Float> breaks = new ArrayList<>();

        float delta = (maxScore-minScore)/nrIntervals;

        float tmpBreak = minScore;
        breaks.add(tmpBreak);
        for (int i=0;i<nrIntervals;i++) {
            tmpBreak += delta;
            breaks.add(tmpBreak);
        }

        return breaks;
    }

    protected int get1DIndex(double score, double minScore, double deltaScore, int maxIndex) {
        int index = (int) Math.floor((score-minScore)/deltaScore);
        if (index<0) index = 0;
        if (index>maxIndex) index = maxIndex;

        return index;
    }

    public boolean isEmpty() {
        return (totDecoyCnt==0f && totTargetCnt==0f);
    }

    public void setCanCalculateFDR(int minNrPsms) {
        canCalculateFDR = (totTargetCnt+totDecoyCnt)>=minNrPsms;
    }

    protected abstract int index(PeptideSpectrumMatch peptideSpectrumMatch);
    protected abstract Set<Integer> getNeighbourIndex(int bin);
    protected abstract void write(File outputFile);
    public abstract void write(String outputDir, String fileTag, String id);
    protected abstract List<Float> getMids(int bin);

    public float getTotTargetCnt() {
        return totTargetCnt;
    }

    public float getTotDecoyCnt() {
        return totDecoyCnt;
    }

    public double getPi_0() {
        return pi_0;
    }

    public double getPi_1() {
        return pi_1;
    }

    public boolean canCalculateFDR() {
        return canCalculateFDR;
    }

    public int getTotNrBins() {
        return totNrBins;
    }

    public int[] getNrBins() {
        return nrBins;
    }

    public int getDimension() {
        return dimension;
    }

    public List<Float> getGamma() {
        return gamma;
    }

    public List<Float> getlFDR() {
        return lFDR;
    }

    public List<Integer> getSortedIndexes() {
        return sortedIndexes;
    }

    public List<Float> getpFDR() {
        return pFDR;
    }

    public List<Float> getTargetCnts() {
        return targetCnts;
    }

    public List<Float> getDecoyCnts() {
        return decoyCnts;
    }

    public List<Integer> getIndexMap() {
        return indexMap;
    }

    public List<Integer> getPsmBins() {
        return psmBins;
    }
}
