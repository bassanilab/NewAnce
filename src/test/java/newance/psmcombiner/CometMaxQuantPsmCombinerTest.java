package newance.psmcombiner;

import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import org.expasy.mzjava.proteomics.mol.Peptide;
import newance.psmconverter.PeptideMatchData;
import newance.util.NewAnceParams;
import org.junit.Assert;
import org.junit.Test;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Created by markusmueller on 02.12.19.
 */
public class CometMaxQuantPsmCombinerTest {

    @Test
    public void test_merge() {

        ConcurrentHashMap<String, List<PeptideMatchData>> cometPsms = getCometPsms();
        ConcurrentHashMap<String, List<PeptideMatchData>> mqPsms = getMQPsms();
        ConcurrentHashMap<String, List<PeptideMatchData>> combined = new ConcurrentHashMap<>();

        CometMaxQuantPsmCombiner combiner = new CometMaxQuantPsmCombiner(mqPsms,combined);
        cometPsms.forEach(combiner);

        Assert.assertTrue(combined.size()==1);
        Assert.assertTrue(combined.get("spec1").size()==2);
        Assert.assertEquals("PEPTIDE",combined.get("spec1").get(0).getPeptide().toString());
        Assert.assertEquals("PEPTIDER",combined.get("spec1").get(1).getPeptide().toString());

    }

    public static ConcurrentHashMap<String, List<PeptideMatchData>> getCometPsms()
    {

        ConcurrentHashMap<String, List<PeptideMatchData>> psms = new ConcurrentHashMap<>();

        Set<String> prots = new HashSet<>();
        prots.add("protein1");
        prots.add("protein2");

        NewAnceParams params = NewAnceParams.getInstance();


        TObjectDoubleMap<String> scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("xcorr",1.0);
        scoreMap.put("deltacn",0.1);
        scoreMap.put("spscore",100.0);

        PeptideMatchData peptideMatchData1 = new PeptideMatchData(Peptide.parse("PEPTIDE"), prots, scoreMap, 1, false);

        scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("xcorr",2.0);
        scoreMap.put("deltacn",0.2);
        scoreMap.put("spscore",200.0);

        PeptideMatchData peptideMatchData2 = new PeptideMatchData(Peptide.parse("PEPTIDER"), prots, scoreMap, 1, false);

        psms.put("spec1",Collections.synchronizedList(new ArrayList<>()));
        psms.get("spec1").add(peptideMatchData1);
        psms.get("spec1").add(peptideMatchData2);

        scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("xcorr",3.0);
        scoreMap.put("deltacn",0.3);
        scoreMap.put("spscore",300.0);

        peptideMatchData1 = new PeptideMatchData(Peptide.parse("PEPTIDE"), prots, scoreMap, 1, false);
        psms.put("spec2",Collections.synchronizedList(new ArrayList<>()));
        psms.get("spec2").add(peptideMatchData1);

        peptideMatchData1 = new PeptideMatchData(Peptide.parse("PEPTIDE"), prots, scoreMap, 1, false);
        psms.put("spec5",Collections.synchronizedList(new ArrayList<>()));
        psms.get("spec5").add(peptideMatchData1);

        return psms;
    }

    public static ConcurrentHashMap<String, List<PeptideMatchData>> getMQPsms()
    {

        ConcurrentHashMap<String, List<PeptideMatchData>> psms = new ConcurrentHashMap<>();

        Set<String> prots = new HashSet<>();
        prots.add("protein1");
        prots.add("protein2");

        NewAnceParams params = NewAnceParams.getInstance();


        TObjectDoubleMap<String> scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("score",1.0);

        PeptideMatchData peptideMatchData1 = new PeptideMatchData(Peptide.parse("PEPTIDE"), prots, scoreMap, 1, false);

        scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("score",2.0);

        PeptideMatchData peptideMatchData2 = new PeptideMatchData(Peptide.parse("PEPTIDER"), prots, scoreMap, 1, false);

        psms.put("spec1",Collections.synchronizedList(new ArrayList<>()));
        psms.get("spec1").add(peptideMatchData1);
        psms.get("spec1").add(peptideMatchData2);

        scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("score",3.0);

        peptideMatchData1 = new PeptideMatchData(Peptide.parse("PEPDIDE"), prots, scoreMap, 1, false);
        psms.put("spec2",Collections.synchronizedList(new ArrayList<>()));
        psms.get("spec2").add(peptideMatchData1);

        scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("score",4.0);

        peptideMatchData1 = new PeptideMatchData(Peptide.parse("PEPTIDE"), prots, scoreMap, 1, false);
        psms.put("spec3",Collections.synchronizedList(new ArrayList<>()));
        psms.get("spec3").add(peptideMatchData1);

        return psms;
    }
}
