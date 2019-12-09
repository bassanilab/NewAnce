package newance.psmcombiner;

import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import org.expasy.mzjava.proteomics.mol.Peptide;
import newance.psmconverter.PeptideMatchData;
import newance.util.NewAnceParams;
import org.junit.Assert;
import org.junit.Test;

import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;

/**
 * Created by markusmueller on 20.11.19.
 */
public class GroupedFDRCalculatorTest {
    @Test
    public void test_print() {

        RegExpProteinGrouper psmGrouper = new RegExpProteinGrouper(Pattern.compile("^sp"),"canonical","cryptic");
        GroupedFDRCalculator groupedFDRCalculator = new GroupedFDRCalculator(psmGrouper);

        System.out.print(groupedFDRCalculator.printTree(1));
    }

    @Test
    public void test_add() {

        RegExpProteinGrouper psmGrouper = new RegExpProteinGrouper(Pattern.compile("^sp"),"canonical","cryptic");
        GroupedFDRCalculator groupedFDRCalculator = new GroupedFDRCalculator(psmGrouper);

        TObjectDoubleMap<String> scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("xcorr", NewAnceParams.getInstance().getMinXCorr());
        scoreMap.put("deltacn", NewAnceParams.getInstance().getMinDeltaCn());
        scoreMap.put("spscore", NewAnceParams.getInstance().getMinSpScore());

        Set<String> prots = new HashSet<>();
        prots.add("sp|protein1");
        prots.add("protein2");

        PeptideMatchData peptideMatchData = new PeptideMatchData(Peptide.parse("PEPTIDE"), prots, scoreMap, 1);

        groupedFDRCalculator.add(peptideMatchData);

        for (String label : groupedFDRCalculator.histogramMap.keySet()) {
            if (label.equals("Z1_canonical") || label.equals("Z1") || label.equals("root"))
                Assert.assertFalse(groupedFDRCalculator.histogramMap.get(label).getScoreHistogram().isEmpty());
            else
                Assert.assertTrue(groupedFDRCalculator.histogramMap.get(label).getScoreHistogram().isEmpty());
        }
    }

    @Test
    public void test_add2() {

        RegExpProteinGrouper psmGrouper = new RegExpProteinGrouper(Pattern.compile("^sp"),"canonical","cryptic");
        GroupedFDRCalculator groupedFDRCalculator = new GroupedFDRCalculator(psmGrouper);

        NewAnceParams params = NewAnceParams.getInstance();
        CometScoreHistogram cometScoreHistogram = CometScoreHistogramTest.buildScoreHisto();

        List<Float> xcorrMids = cometScoreHistogram.calcMids(cometScoreHistogram.calcBreaks((float)params.getMinXCorr(),(float)params.getMaxXCorr(),params.getNrXCorrBins()));
        List<Float> deltacnMids = cometScoreHistogram.calcMids(cometScoreHistogram.calcBreaks((float)params.getMinDeltaCn(),(float)params.getMaxDeltaCn(),params.getNrDeltaCnBins()));
        List<Float> spscoreMids = cometScoreHistogram.calcMids(cometScoreHistogram.calcBreaks((float)params.getMinSpScore(),(float)params.getMaxSpScore(),params.getNrSpScoreBins()));

        int xcorrIdx = xcorrMids.size()/2;
        int deltacnIdx = deltacnMids.size()/2;
        int spscoreIdx = spscoreMids.size()/2;

        Set<String> prots = new HashSet<>();
        prots.add("sp|protein1");
        prots.add("protein2");
        addPsms(xcorrIdx,deltacnIdx,spscoreIdx,groupedFDRCalculator,prots, 10, false);
        addPsms(xcorrIdx,deltacnIdx,spscoreIdx,groupedFDRCalculator,prots, 1, true);
        addPsms(xcorrIdx+5,deltacnIdx+5,spscoreIdx+5,groupedFDRCalculator,prots, 1, false);

        prots = new HashSet<>();
        prots.add("protein3");
        prots.add("protein4");
        addPsms(xcorrIdx,deltacnIdx,spscoreIdx,groupedFDRCalculator,prots, 2, false);
        addPsms(xcorrIdx,deltacnIdx,spscoreIdx,groupedFDRCalculator,prots, 2, true);

        groupedFDRCalculator.setCanCalculateFDR(1);


        for (String label : groupedFDRCalculator.histogramMap.keySet()) {

            ScoreHistogram h = groupedFDRCalculator.histogramMap.get(label).getScoreHistogram();
            Assert.assertEquals(1.0, h.getPi_0()+h.getPi_1(),0.00001);

            if (label.equals("Z1_canonical")) {
                Assert.assertEquals(495f, groupedFDRCalculator.histogramMap.get(label).getScoreHistogram().getTotTargetCnt(),0.00001);
                Assert.assertEquals(45f, groupedFDRCalculator.histogramMap.get(label).getScoreHistogram().getTotDecoyCnt(),0.00001);
            } else if (label.equals("Z1_cryptic")) {
                Assert.assertEquals(90f, groupedFDRCalculator.histogramMap.get(label).getScoreHistogram().getTotTargetCnt(),0.00001);
                Assert.assertEquals(90f, groupedFDRCalculator.histogramMap.get(label).getScoreHistogram().getTotDecoyCnt(),0.00001);
            } else if (label.equals("Z1")) {
                Assert.assertEquals(585f, groupedFDRCalculator.histogramMap.get(label).getScoreHistogram().getTotTargetCnt(),0.00001);
                Assert.assertEquals(135f, groupedFDRCalculator.histogramMap.get(label).getScoreHistogram().getTotDecoyCnt(),0.00001);
            } else if (label.equals("root")) {
                Assert.assertEquals(585f, groupedFDRCalculator.histogramMap.get(label).getScoreHistogram().getTotTargetCnt(),0.00001);
                Assert.assertEquals(135f, groupedFDRCalculator.histogramMap.get(label).getScoreHistogram().getTotDecoyCnt(),0.00001);
            } else {
                Assert.assertEquals(0f, groupedFDRCalculator.histogramMap.get(label).getScoreHistogram().getTotTargetCnt(),0.00001);
                Assert.assertEquals(0f, groupedFDRCalculator.histogramMap.get(label).getScoreHistogram().getTotDecoyCnt(),0.00001);
                Assert.assertEquals(0.5, h.getPi_0(),0.00001);
            }
        }
    }

    @Test
    public void test_calcLocalFDR() {

        RegExpProteinGrouper psmGrouper = new RegExpProteinGrouper(Pattern.compile("^sp"),"canonical","cryptic");
        GroupedFDRCalculator groupedFDRCalculator = new GroupedFDRCalculator(psmGrouper);

        NewAnceParams params = NewAnceParams.getInstance();
        CometScoreHistogram cometScoreHistogram = CometScoreHistogramTest.buildScoreHisto();

        List<Float> xcorrMids = cometScoreHistogram.calcMids(cometScoreHistogram.calcBreaks((float)params.getMinXCorr(),(float)params.getMaxXCorr(),params.getNrXCorrBins()));
        List<Float> deltacnMids = cometScoreHistogram.calcMids(cometScoreHistogram.calcBreaks((float)params.getMinDeltaCn(),(float)params.getMaxDeltaCn(),params.getNrDeltaCnBins()));
        List<Float> spscoreMids = cometScoreHistogram.calcMids(cometScoreHistogram.calcBreaks((float)params.getMinSpScore(),(float)params.getMaxSpScore(),params.getNrSpScoreBins()));

        int xcorrIdx = xcorrMids.size()/2;
        int deltacnIdx = deltacnMids.size()/2;
        int spscoreIdx = spscoreMids.size()/2;

        Set<String> prots = new HashSet<>();
        prots.add("sp|protein1");
        prots.add("protein2");
        addPsms(xcorrIdx, deltacnIdx, spscoreIdx, groupedFDRCalculator, prots, 10, false);
        addPsms(xcorrIdx, deltacnIdx, spscoreIdx, groupedFDRCalculator, prots, 1, true);
        addPsms(xcorrIdx+5,deltacnIdx+5,spscoreIdx+5, groupedFDRCalculator, prots, 1, false);

        prots = new HashSet<>();
        prots.add("protein3");
        prots.add("protein4");
        addPsms(xcorrIdx, deltacnIdx, spscoreIdx, groupedFDRCalculator, prots, 2, false);
        addPsms(xcorrIdx, deltacnIdx, spscoreIdx, groupedFDRCalculator, prots, 2, true);

        groupedFDRCalculator.setCanCalculateFDR(1);

        for (String label : groupedFDRCalculator.histogramMap.keySet()) {

            float[] cnts = groupedFDRCalculator.histogramMap.get(label).getScoreHistogram().getTargetDecoyCounts(1f);
            if (label.equals("Z1_canonical")) {
                Assert.assertEquals(45f, cnts[0],0.00001);
                Assert.assertEquals(495f, cnts[1],0.00001);
            } else if (label.equals("Z1_cryptic")) {
                Assert.assertEquals(90f, cnts[0],0.00001);
                Assert.assertEquals(90f, cnts[1],0.00001);
            } else if (label.equals("Z1")) {
                Assert.assertEquals(135f, cnts[0],0.00001);
                Assert.assertEquals(585f, cnts[1],0.00001);
            } else if (label.equals("root")) {
                Assert.assertEquals(135f, cnts[0],0.00001);
                Assert.assertEquals(585f, cnts[1],0.00001);
            } else {
                Assert.assertEquals(0f, cnts[0],0.00001);
                Assert.assertEquals(0f, cnts[1],0.00001);
            }
        }

        for (String label : groupedFDRCalculator.histogramMap.keySet()) {

            float[] cnts = groupedFDRCalculator.histogramMap.get(label).getScoreHistogram().getTargetDecoyCounts(0.1f);
            if (label.equals("Z1_canonical")) {
                Assert.assertEquals(0f, cnts[0],0.00001);
                Assert.assertEquals(45f, cnts[1],0.00001);
            } else if (label.equals("Z1_cryptic")) {
                Assert.assertEquals(0f, cnts[0],0.00001);
                Assert.assertEquals(0f, cnts[1],0.00001);
            } else if (label.equals("Z1")) {
                Assert.assertEquals(0f, cnts[0],0.00001);
                Assert.assertEquals(45f, cnts[1],0.00001);
            } else if (label.equals("root")) {
                Assert.assertEquals(0f, cnts[0],0.00001);
                Assert.assertEquals(45f, cnts[1],0.00001);
            } else {
                Assert.assertEquals(0f, cnts[0],0.00001);
                Assert.assertEquals(0f, cnts[1],0.00001);
            }
        }


        for (String label : groupedFDRCalculator.histogramMap.keySet()) {

            float[] cnts = groupedFDRCalculator.histogramMap.get(label).getScoreHistogram().getTargetDecoyCounts(0.2f);
            if (label.equals("Z1_canonical")) {
                Assert.assertEquals(0f, cnts[0],0.00001);
                Assert.assertEquals(45f, cnts[1],0.00001);
            } else if (label.equals("Z1_cryptic")) {
                Assert.assertEquals(0f, cnts[0],0.00001);
                Assert.assertEquals(0f, cnts[1],0.00001);
            } else if (label.equals("Z1")) {
                Assert.assertEquals(0f, cnts[0],0.00001);
                Assert.assertEquals(45f, cnts[1],0.00001);
            } else if (label.equals("root")) {
                Assert.assertEquals(0f, cnts[0],0.00001);
                Assert.assertEquals(45f, cnts[1],0.00001);
            } else {
                Assert.assertEquals(0f, cnts[0],0.00001);
                Assert.assertEquals(0f, cnts[1],0.00001);
            }
        }


        for (String label : groupedFDRCalculator.histogramMap.keySet()) {

            float[] cnts = groupedFDRCalculator.histogramMap.get(label).getScoreHistogram().getTargetDecoyCounts(0.5f);
            if (label.equals("Z1_canonical")) {
                Assert.assertEquals(45f, cnts[0],0.00001);
                Assert.assertEquals(495f, cnts[1],0.00001);
            } else if (label.equals("Z1_cryptic")) {
                Assert.assertEquals(0f, cnts[0],0.00001);
                Assert.assertEquals(0f, cnts[1],0.00001);
            } else if (label.equals("Z1")) {
                Assert.assertEquals(135f, cnts[0],0.00001);
                Assert.assertEquals(585f, cnts[1],0.00001);
            } else if (label.equals("root")) {
                Assert.assertEquals(135f, cnts[0],0.00001);
                Assert.assertEquals(585f, cnts[1],0.00001);
            } else {
                Assert.assertEquals(0f, cnts[0],0.00001);
                Assert.assertEquals(0f, cnts[1],0.00001);
            }
        }
    }


    @Test
    public void test_calcGlobalFDR() {

        RegExpProteinGrouper psmGrouper = new RegExpProteinGrouper(Pattern.compile("^sp"),"canonical","cryptic");
        GroupedFDRCalculator groupedFDRCalculator = new GroupedFDRCalculator(psmGrouper);

        groupedFDRCalculator.setCanCalculateFDR(1);

        NewAnceParams params = NewAnceParams.getInstance();
        CometScoreHistogram cometScoreHistogram = CometScoreHistogramTest.buildScoreHisto();

        List<Float> xcorrMids = cometScoreHistogram.calcMids(cometScoreHistogram.calcBreaks((float)params.getMinXCorr(),(float)params.getMaxXCorr(),params.getNrXCorrBins()));
        List<Float> deltacnMids = cometScoreHistogram.calcMids(cometScoreHistogram.calcBreaks((float)params.getMinDeltaCn(),(float)params.getMaxDeltaCn(),params.getNrDeltaCnBins()));
        List<Float> spscoreMids = cometScoreHistogram.calcMids(cometScoreHistogram.calcBreaks((float)params.getMinSpScore(),(float)params.getMaxSpScore(),params.getNrSpScoreBins()));

        int xcorrIdx = xcorrMids.size()/2;
        int deltacnIdx = deltacnMids.size()/2;
        int spscoreIdx = spscoreMids.size()/2;

        Set<String> prots = new HashSet<>();
        prots.add("sp|protein1");
        prots.add("protein2");
        addPsms(xcorrIdx, deltacnIdx, spscoreIdx, groupedFDRCalculator, prots, 10, false);
        addPsms(xcorrIdx, deltacnIdx, spscoreIdx, groupedFDRCalculator, prots, 1, true);
        addPsms(xcorrIdx+5,deltacnIdx+5,spscoreIdx+5, groupedFDRCalculator, prots, 1, false);

        prots = new HashSet<>();
        prots.add("protein3");
        prots.add("protein4");
        addPsms(xcorrIdx, deltacnIdx, spscoreIdx, groupedFDRCalculator, prots, 2, false);
        addPsms(xcorrIdx, deltacnIdx, spscoreIdx, groupedFDRCalculator, prots, 2, true);

        Assert.assertEquals(0f, groupedFDRCalculator.calcGlobalFDR(0f),0.000001);
        Assert.assertEquals(0f, groupedFDRCalculator.calcGlobalFDR(0.1f),0.000001);
        Assert.assertEquals(0f, groupedFDRCalculator.calcGlobalFDR(0.2f),0.000001);
        Assert.assertEquals(7f/22, groupedFDRCalculator.calcGlobalFDR(0.5f),0.000001);
        Assert.assertEquals(270.0f/720, groupedFDRCalculator.calcGlobalFDR(1.0f),0.000001);
    }


    public static void addPsms(int xcorrIdx, int deltacnIdx, int spscoreIdx, GroupedFDRCalculator groupedFDRCalculator, Set<String> prots, int freq, boolean isDecoy)
    {
        NewAnceParams params = NewAnceParams.getInstance();
        CometScoreHistogram cometScoreHistogram = CometScoreHistogramTest.buildScoreHisto();

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

                    n = freq*(3-n)*(3-n);
                    for (int m=0; m<n; m++) groupedFDRCalculator.add(new PeptideMatchData(Peptide.parse("PEPTIDE"), prots, scoreMap, 1, isDecoy));
                }
            }
        }

        groupedFDRCalculator.calcClassProbs();

    }

}
