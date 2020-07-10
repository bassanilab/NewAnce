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

import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import newance.mzjava.mol.Peptide;
import newance.psmconverter.PeptideSpectrumMatch;
import newance.util.NewAnceParams;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * @author Markus Müller
 */

public class ScoreHistogram3DTest {

    @Test
    public void test_getIndex(){

        ScoreHistogram3D cometScoreHistogram = buildScoreHisto();

        TObjectDoubleMap<String> scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("xcorr",NewAnceParams.getInstance().getMinXCorr());
        scoreMap.put("deltacn",NewAnceParams.getInstance().getMinDeltaCn());
        scoreMap.put("spscore",NewAnceParams.getInstance().getMinSpScore());

        List<String> prots = new ArrayList<>();
        prots.add("protein1");
        prots.add("protein2");

        PeptideSpectrumMatch peptideSpectrumMatch = new PeptideSpectrumMatch("spectrumFile",Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, false, false, null);

        Assert.assertEquals(0, cometScoreHistogram.index(peptideSpectrumMatch));

        scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("xcorr",NewAnceParams.getInstance().getMinXCorr()-1.0);
        scoreMap.put("deltacn",NewAnceParams.getInstance().getMinDeltaCn()-1.0);
        scoreMap.put("spscore",NewAnceParams.getInstance().getMinSpScore()-1.0);

        peptideSpectrumMatch = new PeptideSpectrumMatch("spectrumFile", Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, false, false, null);

        Assert.assertEquals(0, cometScoreHistogram.index(peptideSpectrumMatch));

        scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("xcorr",NewAnceParams.getInstance().getMaxXCorr()+10.0);
        scoreMap.put("deltacn",NewAnceParams.getInstance().getMaxDeltaCn()+10.0);
        scoreMap.put("spscore",NewAnceParams.getInstance().getMaxSpScore()+10.0);
        scoreMap.put("neg_log10_p",100);

        PeptideSpectrumMatch peptideSpectrumMatch1 = new PeptideSpectrumMatch("spectrumFile",Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, false, false, null);

        scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("xcorr",NewAnceParams.getInstance().getMaxXCorr()+20.0);
        scoreMap.put("deltacn",NewAnceParams.getInstance().getMaxDeltaCn()+20.0);
        scoreMap.put("spscore",NewAnceParams.getInstance().getMaxSpScore()+40.0);
        scoreMap.put("neg_log10_p",50);
        PeptideSpectrumMatch peptideSpectrumMatch2 = new PeptideSpectrumMatch("spectrumFile",Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, false, false, null);

        int maxIdx = NewAnceParams.getInstance().getNrXCorrBins()*NewAnceParams.getInstance().getNrDeltaCnBins()*
                NewAnceParams.getInstance().getNrSpScoreBins()-1;
        Assert.assertEquals(cometScoreHistogram.index(peptideSpectrumMatch1), cometScoreHistogram.index(peptideSpectrumMatch2));
        Assert.assertEquals(maxIdx, cometScoreHistogram.index(peptideSpectrumMatch1));


        double dXCorr = (NewAnceParams.getInstance().getMaxXCorr()-NewAnceParams.getInstance().getMinXCorr())/NewAnceParams.getInstance().getNrXCorrBins();

        scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("xcorr",NewAnceParams.getInstance().getMinXCorr()+1.4*dXCorr);
        scoreMap.put("deltacn",NewAnceParams.getInstance().getMinDeltaCn());
        scoreMap.put("spscore",NewAnceParams.getInstance().getMinSpScore());

        peptideSpectrumMatch1 = new PeptideSpectrumMatch("spectrumFile",Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, false, false, null);
        Assert.assertEquals(1, cometScoreHistogram.index(peptideSpectrumMatch1));

        scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("xcorr",NewAnceParams.getInstance().getMinXCorr()+5.1*dXCorr);
        scoreMap.put("deltacn",NewAnceParams.getInstance().getMinDeltaCn());
        scoreMap.put("spscore",NewAnceParams.getInstance().getMinSpScore());
        peptideSpectrumMatch2 = new PeptideSpectrumMatch("spectrumFile",Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, false, false, null);

        Assert.assertEquals(5, cometScoreHistogram.index(peptideSpectrumMatch2));

        double dDeltaCn = (NewAnceParams.getInstance().getMaxDeltaCn()-NewAnceParams.getInstance().getMinDeltaCn())/NewAnceParams.getInstance().getNrDeltaCnBins();

        scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("xcorr",NewAnceParams.getInstance().getMinXCorr());
        scoreMap.put("deltacn",NewAnceParams.getInstance().getMinDeltaCn()+1.4*dDeltaCn);
        scoreMap.put("spscore",NewAnceParams.getInstance().getMinSpScore());

        peptideSpectrumMatch1 = new PeptideSpectrumMatch("spectrumFile",Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, false, false, null);
        Assert.assertEquals(NewAnceParams.getInstance().getNrXCorrBins(), cometScoreHistogram.index(peptideSpectrumMatch1));

        scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("xcorr",NewAnceParams.getInstance().getMinXCorr());
        scoreMap.put("deltacn",NewAnceParams.getInstance().getMinDeltaCn()+5.5*dDeltaCn);
        scoreMap.put("spscore",NewAnceParams.getInstance().getMinSpScore());
        peptideSpectrumMatch2 = new PeptideSpectrumMatch("spectrumFile",Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, false, false, null);

        Assert.assertEquals(NewAnceParams.getInstance().getNrXCorrBins()*5, cometScoreHistogram.index(peptideSpectrumMatch2));

        scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("xcorr",NewAnceParams.getInstance().getMinXCorr()+27.5*dXCorr);
        scoreMap.put("deltacn",NewAnceParams.getInstance().getMinDeltaCn()+5.5*dDeltaCn);
        scoreMap.put("spscore",NewAnceParams.getInstance().getMinSpScore());
        peptideSpectrumMatch2 = new PeptideSpectrumMatch("spectrumFile",Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, false, false, null);

        Assert.assertEquals(NewAnceParams.getInstance().getNrXCorrBins()*5+27, cometScoreHistogram.index(peptideSpectrumMatch2));


        double dSpScore = (NewAnceParams.getInstance().getMaxSpScore()-NewAnceParams.getInstance().getMinSpScore())/NewAnceParams.getInstance().getNrSpScoreBins();

        scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("xcorr",NewAnceParams.getInstance().getMinXCorr());
        scoreMap.put("deltacn",NewAnceParams.getInstance().getMinDeltaCn());
        scoreMap.put("spscore",NewAnceParams.getInstance().getMinSpScore()+3.5*dSpScore);
        peptideSpectrumMatch2 = new PeptideSpectrumMatch("spectrumFile",Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, false, false, null);

        Assert.assertEquals(NewAnceParams.getInstance().getNrXCorrBins()*NewAnceParams.getInstance().getNrDeltaCnBins()*3, cometScoreHistogram.index(peptideSpectrumMatch2));
    }

    @Test
    public void test_getNeighbourIndex() {

        ScoreHistogram3D cometScoreHistogram = buildScoreHisto();

        TObjectDoubleMap<String> scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("xcorr", NewAnceParams.getInstance().getMinXCorr());
        scoreMap.put("deltacn", NewAnceParams.getInstance().getMinDeltaCn());
        scoreMap.put("spscore", NewAnceParams.getInstance().getMinSpScore());

        List<String> prots = new ArrayList<>();
        prots.add("protein1");
        prots.add("protein2");

        PeptideSpectrumMatch peptideSpectrumMatch = new PeptideSpectrumMatch("spectrumFile",Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, false, false, null);

        Set<Integer> nn = cometScoreHistogram.getNeighbourIndex(cometScoreHistogram.index(peptideSpectrumMatch));
        Assert.assertEquals(3, nn.size());
        Assert.assertTrue(nn.contains(1));
        Assert.assertTrue(nn.contains(40));
        Assert.assertTrue(nn.contains(1600));

        nn = cometScoreHistogram.getNeighbourIndex(41);
        Assert.assertEquals(5, nn.size());
        Assert.assertTrue(nn.contains(42));
        Assert.assertTrue(nn.contains(40));
        Assert.assertTrue(nn.contains(81));
        Assert.assertTrue(nn.contains(1));
        Assert.assertTrue(nn.contains(1641));

        nn = cometScoreHistogram.getNeighbourIndex(39);
        Assert.assertEquals(3, nn.size());
        Assert.assertTrue(nn.contains(38));
        Assert.assertTrue(nn.contains(79));
        Assert.assertTrue(nn.contains(1639));

        nn = cometScoreHistogram.getNeighbourIndex(63999);
        Assert.assertEquals(3, nn.size());
        Assert.assertTrue(nn.contains(63998));
        Assert.assertTrue(nn.contains(63959));
        Assert.assertTrue(nn.contains(62399));

        nn = cometScoreHistogram.getNeighbourIndex(32000);
        Assert.assertEquals(4, nn.size());
        Assert.assertTrue(nn.contains(32001));
        Assert.assertTrue(nn.contains(32040));
        Assert.assertTrue(nn.contains(33600));
        Assert.assertTrue(nn.contains(30400));


        nn = cometScoreHistogram.getNeighbourIndex(11151);
        Assert.assertEquals(6, nn.size());
        Assert.assertTrue(nn.contains(11152));
        Assert.assertTrue(nn.contains(11150));
        Assert.assertTrue(nn.contains(11191));
        Assert.assertTrue(nn.contains(11111));
        Assert.assertTrue(nn.contains(12751));
        Assert.assertTrue(nn.contains(9551));
    }

    @Test
    public void test_getMids() {

        boolean isDecoy = true;
        ScoreHistogram3D cometScoreHistogram = buildScoreHisto();

        List<String> prots = new ArrayList<>();
        prots.add("protein1");

        NewAnceParams params = NewAnceParams.getInstance();

        List<Float> xcorrMids = cometScoreHistogram.calcMids(cometScoreHistogram.calcBreaks((float)params.getMinXCorr(),(float)params.getMaxXCorr(),params.getNrXCorrBins()));
        List<Float> deltacnMids = cometScoreHistogram.calcMids(cometScoreHistogram.calcBreaks((float)params.getMinDeltaCn(),(float)params.getMaxDeltaCn(),params.getNrDeltaCnBins()));
        List<Float> spscoreMids = cometScoreHistogram.calcMids(cometScoreHistogram.calcBreaks((float)params.getMinSpScore(),(float)params.getMaxSpScore(),params.getNrSpScoreBins()));

        for (int i=0;i<params.getNrXCorrBins();i++) {

            for (int j=0;j<params.getNrDeltaCnBins();j++) {

                for (int k=0;k<params.getNrSpScoreBins();k++) {

                    TObjectDoubleMap<String> scoreMap = new TObjectDoubleHashMap<>();
                    scoreMap.put("xcorr",xcorrMids.get(i));
                    scoreMap.put("deltacn",deltacnMids.get(j));
                    scoreMap.put("spscore",spscoreMids.get(k));

                    PeptideSpectrumMatch peptideSpectrumMatch = new PeptideSpectrumMatch("spectrumFile",Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                            100, 101, 1001.1, isDecoy, false, null);

                    int idx = cometScoreHistogram.index(peptideSpectrumMatch);

                    List<Float> mids = cometScoreHistogram.getMids(idx);

                    Assert.assertEquals(xcorrMids.get(i),mids.get(0),0.0001);
                    Assert.assertEquals(deltacnMids.get(j),mids.get(1),0.0001);
                    Assert.assertEquals(spscoreMids.get(k),mids.get(2),0.0001);
                }
            }
        }
     }

    @Test
    public void addTest() {

        ScoreHistogram3D cometScoreHistogram = buildScoreHisto();

        Assert.assertTrue(cometScoreHistogram.isEmpty());

        List<String> prots = new ArrayList<>();
        prots.add("protein1");
        prots.add("protein2");

        List<String> decoys = new ArrayList<>();
        prots.add("DECOY_p1");
        prots.add("protein2");

        double xcorr = (NewAnceParams.getInstance().getMaxXCorr()-NewAnceParams.getInstance().getMinXCorr())/2.0;
        double deltacn = (NewAnceParams.getInstance().getMaxDeltaCn()-NewAnceParams.getInstance().getMinDeltaCn())/2.0;
        double spscore = (NewAnceParams.getInstance().getMaxSpScore()-NewAnceParams.getInstance().getMinSpScore())/2.0;

        TObjectDoubleMap<String> scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("xcorr",xcorr);
        scoreMap.put("deltacn",deltacn);
        scoreMap.put("spscore",spscore);

        PeptideSpectrumMatch peptideSpectrumMatch1 = new PeptideSpectrumMatch("spectrumFile",Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, false, false, null);

        cometScoreHistogram.add(peptideSpectrumMatch1);
        cometScoreHistogram.calcClassProb();

        Assert.assertFalse(cometScoreHistogram.isEmpty());
        Assert.assertEquals(0.5,cometScoreHistogram.getPi_0(),0.00001);
        Assert.assertEquals(0.5,cometScoreHistogram.getPi_1(),0.00001);
        Assert.assertEquals(0f,cometScoreHistogram.getTotDecoyCnt(),0.0001);
        Assert.assertEquals(1f,cometScoreHistogram.getTotTargetCnt(),0.0001);

        PeptideSpectrumMatch peptideSpectrumMatch2 = new PeptideSpectrumMatch("spectrumFile",Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, true, false, null);

        cometScoreHistogram.add(peptideSpectrumMatch2);
        cometScoreHistogram.calcClassProb();

        Assert.assertFalse(cometScoreHistogram.isEmpty());
        Assert.assertEquals(0.5,cometScoreHistogram.getPi_0(),0.00001);
        Assert.assertEquals(0.5,cometScoreHistogram.getPi_1(),0.00001);
        Assert.assertEquals(1f,cometScoreHistogram.getTotDecoyCnt(),0.0001);
        Assert.assertEquals(1f,cometScoreHistogram.getTotTargetCnt(),0.0001);

        for (int i=0;i<8;i++) cometScoreHistogram.add(peptideSpectrumMatch1);
        cometScoreHistogram.calcClassProb();

        Assert.assertEquals(0.5,cometScoreHistogram.getPi_0(),0.00001);
        Assert.assertEquals(0.5,cometScoreHistogram.getPi_1(),0.00001);
        Assert.assertEquals(1f,cometScoreHistogram.getTotDecoyCnt(),0.0001);
        Assert.assertEquals(9f,cometScoreHistogram.getTotTargetCnt(),0.0001);
    }

    @Test
    public void test_Smooth() {
        smooth(true);
        smooth(false);
    }

    public void smooth(boolean isDecoy) {

        ScoreHistogram3D cometScoreHistogram = buildScoreHisto();

        List<String> prots = new ArrayList<>();
        prots.add("protein1");
        prots.add("protein2");

        NewAnceParams params = NewAnceParams.getInstance();

        List<Float> xcorrMids = cometScoreHistogram.calcMids(cometScoreHistogram.calcBreaks((float)params.getMinXCorr(),(float)params.getMaxXCorr(),params.getNrXCorrBins()));
        List<Float> deltacnMids = cometScoreHistogram.calcMids(cometScoreHistogram.calcBreaks((float)params.getMinDeltaCn(),(float)params.getMaxDeltaCn(),params.getNrDeltaCnBins()));
        List<Float> spscoreMids = cometScoreHistogram.calcMids(cometScoreHistogram.calcBreaks((float)params.getMinSpScore(),(float)params.getMaxSpScore(),params.getNrSpScoreBins()));
        int xcorrIdx = xcorrMids.size()/2;
        int deltacnIdx = deltacnMids.size()/2;
        int spscoreIdx = spscoreMids.size()/2;

        List<Integer> indexes = addPsms(xcorrIdx, deltacnIdx, spscoreIdx, cometScoreHistogram, isDecoy);

        if (isDecoy) {

            Assert.assertEquals(1.0,cometScoreHistogram.getPi_0(),0.00001);
            Assert.assertEquals(0.0,cometScoreHistogram.getPi_1(),0.00001);
            Assert.assertEquals(14f,cometScoreHistogram.getTotDecoyCnt(),0.0001);
            Assert.assertEquals(0f,cometScoreHistogram.getTotTargetCnt(),0.0001);
        } else {

            Assert.assertEquals(0.5,cometScoreHistogram.getPi_0(),0.00001);
            Assert.assertEquals(0.5,cometScoreHistogram.getPi_1(),0.00001);
            Assert.assertEquals(0f,cometScoreHistogram.getTotDecoyCnt(),0.0001);
            Assert.assertEquals(14f,cometScoreHistogram.getTotTargetCnt(),0.0001);
        }

        cometScoreHistogram.smoothHistogram(true);
        SmoothedScoreHistogram smoothed = cometScoreHistogram.smoothedHistogram;

        if (isDecoy) {
            Assert.assertEquals(2.0,smoothed.decoyCnts.get(smoothed.indexMap.get(indexes.get(0))),0.0001);
            Assert.assertEquals(1.285714,smoothed.decoyCnts.get(smoothed.indexMap.get(indexes.get(1))),0.0001);
            Assert.assertEquals(1.285714,smoothed.decoyCnts.get(smoothed.indexMap.get(indexes.get(2))),0.0001);
            Assert.assertEquals(1.285714,smoothed.decoyCnts.get(smoothed.indexMap.get(indexes.get(3))),0.0001);
            Assert.assertEquals(1.285714,smoothed.decoyCnts.get(smoothed.indexMap.get(indexes.get(4))),0.0001);
            Assert.assertEquals(1.285714,smoothed.decoyCnts.get(smoothed.indexMap.get(indexes.get(5))),0.0001);
            Assert.assertEquals(1.285714,smoothed.decoyCnts.get(smoothed.indexMap.get(indexes.get(6))),0.0001);

        } else {

            Assert.assertEquals(2.0,smoothed.targetCnts.get(smoothed.indexMap.get(indexes.get(0))),0.0001);
            Assert.assertEquals(1.285714,smoothed.targetCnts.get(smoothed.indexMap.get(indexes.get(1))),0.0001);
            Assert.assertEquals(1.285714,smoothed.targetCnts.get(smoothed.indexMap.get(indexes.get(2))),0.0001);
            Assert.assertEquals(1.285714,smoothed.targetCnts.get(smoothed.indexMap.get(indexes.get(3))),0.0001);
            Assert.assertEquals(1.285714,smoothed.targetCnts.get(smoothed.indexMap.get(indexes.get(4))),0.0001);
            Assert.assertEquals(1.285714,smoothed.targetCnts.get(smoothed.indexMap.get(indexes.get(5))),0.0001);
            Assert.assertEquals(1.285714,smoothed.targetCnts.get(smoothed.indexMap.get(indexes.get(6))),0.0001);
        }
    }

    @Test
    public void test_Smooth2() {
        smooth2(true);
        smooth2(false);
    }

    public void smooth2(boolean isDecoy) {

        ScoreHistogram3D cometScoreHistogram = buildScoreHisto();

        int xcorrIdx = 0;
        int deltacnIdx = 0;
        int spscoreIdx = 0;

        List<Integer> indexes = addPsms(xcorrIdx, deltacnIdx, spscoreIdx, cometScoreHistogram, isDecoy);

        cometScoreHistogram.smoothHistogram(false);
        cometScoreHistogram.calcClassProb();

        SmoothedScoreHistogram smoothed = cometScoreHistogram.smoothedHistogram;

        if (isDecoy) {

            Assert.assertEquals(1.0,smoothed.getPi_0(),0.00001);
            Assert.assertEquals(0.0,cometScoreHistogram.getPi_1(),0.00001);
            Assert.assertEquals(11f,cometScoreHistogram.getTotDecoyCnt(),0.0001);
            Assert.assertEquals(0f,cometScoreHistogram.getTotTargetCnt(),0.0001);

            Assert.assertEquals(11f/4,smoothed.decoyCnts.get(smoothed.indexMap.get(indexes.get(0))),0.0001);
            Assert.assertEquals(9f/5,smoothed.decoyCnts.get(smoothed.indexMap.get(indexes.get(1))),0.0001);
            Assert.assertEquals(9f/5,smoothed.decoyCnts.get(smoothed.indexMap.get(indexes.get(2))),0.0001);
            Assert.assertEquals(9f/5,smoothed.decoyCnts.get(smoothed.indexMap.get(indexes.get(3))),0.0001);

            Assert.assertEquals(10,smoothed.decoyCnts.size());
        } else {

            Assert.assertEquals(0.5,cometScoreHistogram.getPi_0(),0.00001);
            Assert.assertEquals(0.5,cometScoreHistogram.getPi_1(),0.00001);
            Assert.assertEquals(0f,cometScoreHistogram.getTotDecoyCnt(),0.0001);
            Assert.assertEquals(11f,cometScoreHistogram.getTotTargetCnt(),0.0001);

            Assert.assertEquals(11f/4,smoothed.targetCnts.get(smoothed.indexMap.get(indexes.get(0))),0.0001);
            Assert.assertEquals(9f/5,smoothed.targetCnts.get(smoothed.indexMap.get(indexes.get(1))),0.0001);
            Assert.assertEquals(9f/5,smoothed.targetCnts.get(smoothed.indexMap.get(indexes.get(2))),0.0001);
            Assert.assertEquals(9f/5,smoothed.targetCnts.get(smoothed.indexMap.get(indexes.get(3))),0.0001);

            Assert.assertEquals(10,smoothed.targetCnts.size());
        }

    }

    @Test
    public void test_Smooth3() {
        smooth3(true);
        smooth3(false);
    }

    public void smooth3(boolean isDecoy) {

        ScoreHistogram3D cometScoreHistogram = buildScoreHisto();

        int xcorrIdx = 0;
        int deltacnIdx = 0;
        int spscoreIdx = 0;

        List<Integer> indexes = addPsms(xcorrIdx, deltacnIdx, spscoreIdx, cometScoreHistogram, isDecoy);

        cometScoreHistogram.smoothHistogram(true);
        cometScoreHistogram.calcClassProb();
        SmoothedScoreHistogram smoothed = cometScoreHistogram.smoothedHistogram;

        if (isDecoy) {

            Assert.assertEquals(1.0,cometScoreHistogram.getPi_0(),0.00001);
            Assert.assertEquals(0.0,cometScoreHistogram.getPi_1(),0.00001);
            Assert.assertEquals(11f,cometScoreHistogram.getTotDecoyCnt(),0.0001);
            Assert.assertEquals(0f,cometScoreHistogram.getTotTargetCnt(),0.0001);

            Assert.assertEquals(3.102564,smoothed.decoyCnts.get(smoothed.indexMap.get(indexes.get(0))),0.0001);
            Assert.assertEquals(2.0307693,smoothed.decoyCnts.get(smoothed.indexMap.get(indexes.get(1))),0.0001);
            Assert.assertEquals(2.0307693,smoothed.decoyCnts.get(smoothed.indexMap.get(indexes.get(2))),0.0001);
            Assert.assertEquals(2.0307693,smoothed.decoyCnts.get(smoothed.indexMap.get(indexes.get(3))),0.0001);

            Assert.assertEquals(10,smoothed.decoyCnts.size());
        } else {

            Assert.assertEquals(0.5,cometScoreHistogram.getPi_0(),0.00001);
            Assert.assertEquals(0.5,cometScoreHistogram.getPi_1(),0.00001);
            Assert.assertEquals(0f,cometScoreHistogram.getTotDecoyCnt(),0.0001);
            Assert.assertEquals(11f,cometScoreHistogram.getTotTargetCnt(),0.0001);

            Assert.assertEquals(3.102564,smoothed.targetCnts.get(smoothed.indexMap.get(indexes.get(0))),0.0001);
            Assert.assertEquals(2.0307693,smoothed.targetCnts.get(smoothed.indexMap.get(indexes.get(1))),0.0001);
            Assert.assertEquals(2.0307693,smoothed.targetCnts.get(smoothed.indexMap.get(indexes.get(2))),0.0001);
            Assert.assertEquals(2.0307693,smoothed.targetCnts.get(smoothed.indexMap.get(indexes.get(3))),0.0001);

            Assert.assertEquals(10,smoothed.targetCnts.size());
        }

    }

    @Test
    public void test_removeSpikeNoise() {

        ScoreHistogram3D cometScoreHistogram = buildScoreHisto();

        List<String> prots = new ArrayList<>();
        prots.add("protein1");
        prots.add("protein2");

        NewAnceParams params = NewAnceParams.getInstance();

        List<Float> xcorrMids = cometScoreHistogram.calcMids(cometScoreHistogram.calcBreaks((float)params.getMinXCorr(),(float)params.getMaxXCorr(),params.getNrXCorrBins()));
        List<Float> deltacnMids = cometScoreHistogram.calcMids(cometScoreHistogram.calcBreaks((float)params.getMinDeltaCn(),(float)params.getMaxDeltaCn(),params.getNrDeltaCnBins()));
        List<Float> spscoreMids = cometScoreHistogram.calcMids(cometScoreHistogram.calcBreaks((float)params.getMinSpScore(),(float)params.getMaxSpScore(),params.getNrSpScoreBins()));
        int xcorrIdx = xcorrMids.size()/2;
        int deltacnIdx = deltacnMids.size()/2;
        int spscoreIdx = spscoreMids.size()/2;

        TObjectDoubleMap<String> scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("xcorr",xcorrMids.get(xcorrIdx));
        scoreMap.put("deltacn",deltacnMids.get(deltacnIdx));
        scoreMap.put("spscore",spscoreMids.get(spscoreIdx));

        PeptideSpectrumMatch peptideSpectrumMatch = new PeptideSpectrumMatch("spectrumFile",Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, false, false, null);

        for (int i=0;i<8;i++) cometScoreHistogram.add(peptideSpectrumMatch);
        int idx1 = cometScoreHistogram.index(peptideSpectrumMatch);

        Assert.assertEquals(8f,cometScoreHistogram.targetCnts.get(cometScoreHistogram.indexMap.get(idx1)),0.0001);
        Assert.assertTrue(cometScoreHistogram.getTotTargetCnt()==8);

        cometScoreHistogram.removeSpikeNoiseHistogram(true);

        Assert.assertEquals(8f,cometScoreHistogram.targetCnts.get(cometScoreHistogram.indexMap.get(idx1)),0.0001);
        Assert.assertTrue(cometScoreHistogram.getTotTargetCnt()==8);

        SmoothedScoreHistogram smoothed = cometScoreHistogram.smoothedHistogram;
        smoothed.adjustTotalCounts();

        Assert.assertTrue(smoothed.indexMap.get(idx1)==-1);
        Assert.assertTrue(smoothed.getTotTargetCnt()==0);
    }

    @Test
    public void test_Smooth4() {

        boolean isDecoy = true;

        ScoreHistogram3D cometScoreHistogram = buildScoreHisto();

        NewAnceParams params = NewAnceParams.getInstance();

        List<Float> xcorrMids = cometScoreHistogram.calcMids(cometScoreHistogram.calcBreaks((float)params.getMinXCorr(),(float)params.getMaxXCorr(),params.getNrXCorrBins()));
        List<Float> deltacnMids = cometScoreHistogram.calcMids(cometScoreHistogram.calcBreaks((float)params.getMinDeltaCn(),(float)params.getMaxDeltaCn(),params.getNrDeltaCnBins()));
        List<Float> spscoreMids = cometScoreHistogram.calcMids(cometScoreHistogram.calcBreaks((float)params.getMinSpScore(),(float)params.getMaxSpScore(),params.getNrSpScoreBins()));

        int xcorrIdx = xcorrMids.size()/2;
        int deltacnIdx = deltacnMids.size()/2;
        int spscoreIdx = spscoreMids.size()/2;

        List<Integer> indexes = addPsms(xcorrIdx, deltacnIdx, spscoreIdx, cometScoreHistogram, isDecoy);

        float value1 = cometScoreHistogram.decoyCnts.get(cometScoreHistogram.indexMap.get(indexes.get(0)));
        float value2 = cometScoreHistogram.decoyCnts.get(cometScoreHistogram.indexMap.get(indexes.get(1)));
        float value3 = cometScoreHistogram.decoyCnts.get(cometScoreHistogram.indexMap.get(indexes.get(2)));
        float value4 = cometScoreHistogram.decoyCnts.get(cometScoreHistogram.indexMap.get(indexes.get(3)));
        float value5 = cometScoreHistogram.decoyCnts.get(cometScoreHistogram.indexMap.get(indexes.get(4)));
        float value6 = cometScoreHistogram.decoyCnts.get(cometScoreHistogram.indexMap.get(indexes.get(5)));
        float value7 = cometScoreHistogram.decoyCnts.get(cometScoreHistogram.indexMap.get(indexes.get(6)));

        cometScoreHistogram.smoothHistogram(true);
        SmoothedScoreHistogram smoothed = cometScoreHistogram.smoothedHistogram;

        float value1_1 = smoothed.decoyCnts.get(smoothed.indexMap.get(indexes.get(0)));
        float value2_1 = smoothed.decoyCnts.get(smoothed.indexMap.get(indexes.get(1)));
        float value3_1 = smoothed.decoyCnts.get(smoothed.indexMap.get(indexes.get(2)));
        float value4_1 = smoothed.decoyCnts.get(smoothed.indexMap.get(indexes.get(3)));
        float value5_1 = smoothed.decoyCnts.get(smoothed.indexMap.get(indexes.get(4)));
        float value6_1 = smoothed.decoyCnts.get(smoothed.indexMap.get(indexes.get(5)));
        float value7_1 = smoothed.decoyCnts.get(smoothed.indexMap.get(indexes.get(6)));

        Assert.assertTrue(value1>value1_1);
        Assert.assertTrue(value2<value2_1);
        Assert.assertTrue(value3<value3_1);
        Assert.assertTrue(value4<value4_1);
        Assert.assertTrue(value5<value5_1);
        Assert.assertTrue(value6<value6_1);
        Assert.assertTrue(value7<value7_1);

        List<String> prots = new ArrayList<>();
        prots.add("protein1");
        prots.add("protein2");

        TObjectDoubleMap<String> scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("xcorr",xcorrMids.get(xcorrIdx+2));
        scoreMap.put("deltacn",deltacnMids.get(deltacnIdx));
        scoreMap.put("spscore",spscoreMids.get(spscoreIdx));

        PeptideSpectrumMatch peptideSpectrumMatch = new PeptideSpectrumMatch("spectrumFile",Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, isDecoy, false, null);
        int idx8 = cometScoreHistogram.index(peptideSpectrumMatch);

        scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("xcorr",xcorrMids.get(xcorrIdx));
        scoreMap.put("deltacn",deltacnMids.get(deltacnIdx-2));
        scoreMap.put("spscore",spscoreMids.get(spscoreIdx));

        peptideSpectrumMatch = new PeptideSpectrumMatch("spectrumFile",Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, isDecoy, false, null);
        int idx9 = cometScoreHistogram.index(peptideSpectrumMatch);

        Assert.assertTrue(smoothed.indexMap.get(idx8)>0);
        Assert.assertTrue(smoothed.indexMap.get(idx9)>0);

        float value8_1 = smoothed.decoyCnts.get(smoothed.indexMap.get(idx8));
        float value9_1 = smoothed.decoyCnts.get(smoothed.indexMap.get(idx9));

        Assert.assertTrue(value8_1>0);
        Assert.assertTrue(value9_1>0);

        cometScoreHistogram.smoothHistogram(true);
        smoothed = cometScoreHistogram.smoothedHistogram;

        float value1_2 = smoothed.decoyCnts.get(smoothed.indexMap.get(indexes.get(0)));
        float value2_2 = smoothed.decoyCnts.get(smoothed.indexMap.get(indexes.get(1)));
        float value3_2 = smoothed.decoyCnts.get(smoothed.indexMap.get(indexes.get(2)));
        float value4_2 = smoothed.decoyCnts.get(smoothed.indexMap.get(indexes.get(3)));
        float value5_2 = smoothed.decoyCnts.get(smoothed.indexMap.get(indexes.get(4)));
        float value6_2 = smoothed.decoyCnts.get(smoothed.indexMap.get(indexes.get(5)));
        float value7_2 = smoothed.decoyCnts.get(smoothed.indexMap.get(indexes.get(6)));

        Assert.assertTrue(value1_2<value1_1);
        Assert.assertTrue(value2_2<value2_1);
        Assert.assertTrue(value3_2<value3_1);
        Assert.assertTrue(value4_2<value4_1);
        Assert.assertTrue(value5_2<value5_1);
        Assert.assertTrue(value6_2<value6_1);
        Assert.assertTrue(value7_2<value7_1);

        float value8_2 = smoothed.decoyCnts.get(smoothed.indexMap.get(idx8));
        float value9_2 = smoothed.decoyCnts.get(smoothed.indexMap.get(idx9));

        Assert.assertTrue(value8_2>0);
        Assert.assertTrue(value9_2>0);
    }


    @Test
    public void test_getTargetDecoyCounts() {

        ScoreHistogram3D cometScoreHistogram = buildScoreHisto();

        List<String> prots = new ArrayList<>();
        prots.add("protein1");
        prots.add("protein2");

        NewAnceParams params = NewAnceParams.getInstance();

        List<Float> xcorrMids = cometScoreHistogram.calcMids(cometScoreHistogram.calcBreaks((float)params.getMinXCorr(),(float)params.getMaxXCorr(),params.getNrXCorrBins()));
        List<Float> deltacnMids = cometScoreHistogram.calcMids(cometScoreHistogram.calcBreaks((float)params.getMinDeltaCn(),(float)params.getMaxDeltaCn(),params.getNrDeltaCnBins()));
        List<Float> spscoreMids = cometScoreHistogram.calcMids(cometScoreHistogram.calcBreaks((float)params.getMinSpScore(),(float)params.getMaxSpScore(),params.getNrSpScoreBins()));

        int xcorrIdx = xcorrMids.size()/2;
        int deltacnIdx = deltacnMids.size()/2;
        int spscoreIdx = spscoreMids.size()/2;

        List<PeptideSpectrumMatch> psms = new ArrayList<>();

        addPsms2(xcorrIdx, deltacnIdx, spscoreIdx, cometScoreHistogram, prots, psms,  10, false, 2);
        addPsms2(xcorrIdx, deltacnIdx, spscoreIdx, cometScoreHistogram, prots, psms,  1, true,3);
        addPsms2(xcorrIdx+5,deltacnIdx+5,spscoreIdx+5, cometScoreHistogram, prots, psms, 1, false, 2);

        float[] counts = cometScoreHistogram.getTargetDecoyCounts(0.2f);

        int tCnt = 0;
        int dCnt = 0;
        for (PeptideSpectrumMatch psm : psms) {

            if (cometScoreHistogram.getLocalFDR(psm) > 0.2) continue;

            if (psm.isDecoy()) dCnt++;
            else tCnt++;
        }

        Assert.assertEquals(dCnt,counts[0],0.00001);
        Assert.assertEquals(tCnt,counts[1],0.00001);

        cometScoreHistogram.smoothHistogram(true);
        counts = cometScoreHistogram.getTargetDecoyCounts(0.2f);

        Assert.assertEquals(dCnt,counts[0],0.00001);
        Assert.assertEquals(tCnt,counts[1],0.00001);
    }


    @Test
    public void test_calcGamma() {

        ScoreHistogram3D cometScoreHistogram = buildScoreHisto();

        NewAnceParams params = NewAnceParams.getInstance();

        List<Float> xcorrMids = cometScoreHistogram.calcMids(cometScoreHistogram.calcBreaks((float)params.getMinXCorr(),(float)params.getMaxXCorr(),params.getNrXCorrBins()));
        List<Float> deltacnMids = cometScoreHistogram.calcMids(cometScoreHistogram.calcBreaks((float)params.getMinDeltaCn(),(float)params.getMaxDeltaCn(),params.getNrDeltaCnBins()));
        List<Float> spscoreMids = cometScoreHistogram.calcMids(cometScoreHistogram.calcBreaks((float)params.getMinSpScore(),(float)params.getMaxSpScore(),params.getNrSpScoreBins()));

        int xcorrIdx = xcorrMids.size()/2;
        int deltacnIdx = deltacnMids.size()/2;
        int spscoreIdx = spscoreMids.size()/2;

        List<Integer> indexes = addPsms(xcorrIdx, deltacnIdx, spscoreIdx, cometScoreHistogram, false);
        addPsms(xcorrIdx, deltacnIdx, spscoreIdx, cometScoreHistogram, true);

        cometScoreHistogram.calcGamma();

        Set<Integer> indexSet = new HashSet<>(indexes);
        for (int i=0;i<cometScoreHistogram.getTotNrBins();i++) {
            if (indexSet.contains(i))
                Assert.assertEquals(1f,cometScoreHistogram.getGamma(i),0.0001);
            else
                Assert.assertEquals(-1f,cometScoreHistogram.getGamma(0),0.0001);
        }
    }

    @Test
    public void test_psmCounts() {

        ScoreHistogram3D scoreHistogram3D = buildScoreHisto();

        NewAnceParams params = NewAnceParams.getInstance();

        List<Float> xcorrMids = scoreHistogram3D.calcMids(scoreHistogram3D.calcBreaks((float)params.getMinXCorr(),(float)params.getMaxXCorr(),params.getNrXCorrBins()));
        List<Float> deltacnMids = scoreHistogram3D.calcMids(scoreHistogram3D.calcBreaks((float)params.getMinDeltaCn(),(float)params.getMaxDeltaCn(),params.getNrDeltaCnBins()));
        List<Float> spscoreMids = scoreHistogram3D.calcMids(scoreHistogram3D.calcBreaks((float)params.getMinSpScore(),(float)params.getMaxSpScore(),params.getNrSpScoreBins()));

        int xcorrIdx = xcorrMids.size()/2;
        int deltacnIdx = deltacnMids.size()/2;
        int spscoreIdx = spscoreMids.size()/2;

        addPsms(xcorrIdx, deltacnIdx, spscoreIdx, scoreHistogram3D, false);
        addPsms(xcorrIdx, deltacnIdx, spscoreIdx, scoreHistogram3D, true);

        scoreHistogram3D.calcGamma();
        float[] psmCounts = scoreHistogram3D.getTargetDecoyCounts(0f);

        Assert.assertEquals(0,psmCounts[0],0.00001);
        Assert.assertEquals(0,psmCounts[1],0.00001);

        psmCounts = scoreHistogram3D.getTargetDecoyCounts(1.1f);

        Assert.assertEquals(14,psmCounts[0],0.00001);
        Assert.assertEquals(14,psmCounts[1],0.00001);

        scoreHistogram3D.smoothHistogram(true);

        psmCounts = scoreHistogram3D.getTargetDecoyCounts(1.1f);

        Assert.assertEquals(14,psmCounts[0],0.00001);
        Assert.assertEquals(14,psmCounts[1],0.00001);
    }

    public static List<Integer> addPsms(int xcorrIdx, int deltacnIdx, int spscoreIdx,
                                        ScoreHistogram3D scoreHistogram3D, boolean isDecoy)
    {
        List<String> prots = new ArrayList<>();
        prots.add("protein1");
        prots.add("protein2");

        NewAnceParams params = NewAnceParams.getInstance();

        List<Float> xcorrMids = scoreHistogram3D.calcMids(scoreHistogram3D.calcBreaks((float)params.getMinXCorr(),(float)params.getMaxXCorr(),params.getNrXCorrBins()));
        List<Float> deltacnMids = scoreHistogram3D.calcMids(scoreHistogram3D.calcBreaks((float)params.getMinDeltaCn(),(float)params.getMaxDeltaCn(),params.getNrDeltaCnBins()));
        List<Float> spscoreMids = scoreHistogram3D.calcMids(scoreHistogram3D.calcBreaks((float)params.getMinSpScore(),(float)params.getMaxSpScore(),params.getNrSpScoreBins()));

        if (xcorrIdx<0 || xcorrIdx>=xcorrMids.size()) return null;
        if (deltacnIdx<0 || deltacnIdx>=deltacnMids.size()) return null;
        if (spscoreIdx<0 || spscoreIdx>=spscoreMids.size()) return null;

        List<Integer> indexes = new ArrayList<>();

        TObjectDoubleMap<String> scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("xcorr",xcorrMids.get(xcorrIdx));
        scoreMap.put("deltacn",deltacnMids.get(deltacnIdx));
        scoreMap.put("spscore",spscoreMids.get(spscoreIdx));

        PeptideSpectrumMatch peptideSpectrumMatch = new PeptideSpectrumMatch("spectrumFile",Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, isDecoy, false, null);

        for (int i=0;i<8;i++) scoreHistogram3D.add(peptideSpectrumMatch);
        int idx1 = scoreHistogram3D.index(peptideSpectrumMatch);
        indexes.add(idx1);

        if (xcorrIdx+1<xcorrMids.size()) {
            scoreMap = new TObjectDoubleHashMap<>();
            scoreMap.put("xcorr", xcorrMids.get(xcorrIdx + 1));
            scoreMap.put("deltacn", deltacnMids.get(deltacnIdx));
            scoreMap.put("spscore", spscoreMids.get(spscoreIdx));

            peptideSpectrumMatch = new PeptideSpectrumMatch("spectrumFile",Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                    100, 101, 1001.1, isDecoy, false, null);
            scoreHistogram3D.add(peptideSpectrumMatch);
            int idx2 = scoreHistogram3D.index(peptideSpectrumMatch);
            indexes.add(idx2);
        }

        if (xcorrIdx-1>=0) {
            scoreMap = new TObjectDoubleHashMap<>();
            scoreMap.put("xcorr", xcorrMids.get(xcorrIdx - 1));
            scoreMap.put("deltacn", deltacnMids.get(deltacnIdx));
            scoreMap.put("spscore", spscoreMids.get(spscoreIdx));

            peptideSpectrumMatch = new PeptideSpectrumMatch("spectrumFile",Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                    100, 101, 1001.1, isDecoy, false, null);
            scoreHistogram3D.add(peptideSpectrumMatch);
            int idx3 = scoreHistogram3D.index(peptideSpectrumMatch);
            indexes.add(idx3);
        }

        if (deltacnIdx+1<deltacnMids.size()) {
            scoreMap = new TObjectDoubleHashMap<>();
            scoreMap.put("xcorr",xcorrMids.get(xcorrIdx));
            scoreMap.put("deltacn",deltacnMids.get(deltacnIdx+1));
            scoreMap.put("spscore",spscoreMids.get(spscoreIdx));

            peptideSpectrumMatch = new PeptideSpectrumMatch("spectrumFile",Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                    100, 101, 1001.1, isDecoy, false, null);
            scoreHistogram3D.add(peptideSpectrumMatch);
            int idx4 = scoreHistogram3D.index(peptideSpectrumMatch);
            indexes.add(idx4);
        }

        if (deltacnIdx-1>=0) {
            scoreMap = new TObjectDoubleHashMap<>();
            scoreMap.put("xcorr", xcorrMids.get(xcorrIdx));
            scoreMap.put("deltacn", deltacnMids.get(deltacnIdx - 1));
            scoreMap.put("spscore", spscoreMids.get(spscoreIdx));

            peptideSpectrumMatch = new PeptideSpectrumMatch("spectrumFile",Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                    100, 101, 1001.1, isDecoy, false, null);
            scoreHistogram3D.add(peptideSpectrumMatch);
            int idx5 = scoreHistogram3D.index(peptideSpectrumMatch);
            indexes.add(idx5);
        }

        if (spscoreIdx+1<spscoreMids.size()) {
            scoreMap = new TObjectDoubleHashMap<>();
            scoreMap.put("xcorr", xcorrMids.get(xcorrIdx));
            scoreMap.put("deltacn", deltacnMids.get(deltacnIdx));
            scoreMap.put("spscore", spscoreMids.get(spscoreIdx + 1));

            peptideSpectrumMatch = new PeptideSpectrumMatch("spectrumFile",Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                    100, 101, 1001.1, isDecoy, false, null);
            scoreHistogram3D.add(peptideSpectrumMatch);
            int idx6 = scoreHistogram3D.index(peptideSpectrumMatch);
            indexes.add(idx6);
        }

        if (spscoreIdx-1>=0) {
            scoreMap = new TObjectDoubleHashMap<>();
            scoreMap.put("xcorr", xcorrMids.get(xcorrIdx));
            scoreMap.put("deltacn", deltacnMids.get(deltacnIdx));
            scoreMap.put("spscore", spscoreMids.get(spscoreIdx - 1));

            peptideSpectrumMatch = new PeptideSpectrumMatch("spectrumFile",Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                    100, 101, 1001.1, isDecoy, false, null);
            scoreHistogram3D.add(peptideSpectrumMatch);
            int idx7 = scoreHistogram3D.index(peptideSpectrumMatch);
            indexes.add(idx7);
        }

        scoreHistogram3D.calcClassProb();

        return indexes;
    }

    public static void addPsms2(int xcorrIdx, int deltacnIdx, int spscoreIdx,
                                ScoreHistogram3D cometScoreHistogram, List<String> prots,
                                List<PeptideSpectrumMatch> psmList,
                                int freq, boolean isDecoy, int degree)
    {
        NewAnceParams params = NewAnceParams.getInstance();

        List<Float> xcorrMids = cometScoreHistogram.calcMids(cometScoreHistogram.calcBreaks((float)params.getMinXCorr(),(float)params.getMaxXCorr(),params.getNrXCorrBins()));
        List<Float> deltacnMids = cometScoreHistogram.calcMids(cometScoreHistogram.calcBreaks((float)params.getMinDeltaCn(),(float)params.getMaxDeltaCn(),params.getNrDeltaCnBins()));
        List<Float> spscoreMids = cometScoreHistogram.calcMids(cometScoreHistogram.calcBreaks((float)params.getMinSpScore(),(float)params.getMaxSpScore(),params.getNrSpScoreBins()));

        TObjectDoubleMap<String> scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("xcorr",xcorrMids.get(xcorrIdx));
        scoreMap.put("deltacn",deltacnMids.get(deltacnIdx));
        scoreMap.put("spscore",spscoreMids.get(spscoreIdx));


        for (int i=-1;i<=1;i++) {
            for (int j=-1;j<=1;j++) {
                for (int k=-1;k<=1;k++) {

                    if (xcorrIdx+i<0 || xcorrIdx+i>=xcorrMids.size()) continue;
                    if (deltacnIdx+j<0 || deltacnIdx+j>=deltacnMids.size()) continue;
                    if (spscoreIdx+k<0 || spscoreIdx+k>=spscoreMids.size()) continue;

                    scoreMap = new TObjectDoubleHashMap<>();
                    scoreMap.put("xcorr", xcorrMids.get(xcorrIdx + i));
                    scoreMap.put("deltacn", deltacnMids.get(deltacnIdx + j));
                    scoreMap.put("spscore", spscoreMids.get(spscoreIdx + k));

                    int n = (i<0)?-i:i;
                    n += (j<0)?-j:j;
                    n += (k<0)?-k:k;

                    int m = 1;
                    for (int r=0;r<degree;r++) m *= (3-n);
                    m *= freq;
                    for (int r=0; r<m; r++) {
                        PeptideSpectrumMatch psm = new PeptideSpectrumMatch("spectrumFile",Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                                100, 101, 1001.1, isDecoy, false, null);
                        cometScoreHistogram.add(psm);
                        psmList.add(psm);
                    }
                }
            }
        }

        cometScoreHistogram.calcClassProb();
    }


    public static ScoreHistogram3D buildScoreHisto() {

        return NewAnceParams.getInstance().getScoreHistogram3D(NewAnceParams.SearchTool.COMET);
    }
}
