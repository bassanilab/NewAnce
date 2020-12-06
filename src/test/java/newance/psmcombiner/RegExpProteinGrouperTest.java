package newance.psmcombiner;

import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import newance.mzjava.mol.Peptide;
import newance.psmconverter.PeptideSpectrumMatch;
import newance.util.NewAnceParams;
import org.junit.Assert;
import org.junit.Test;

import java.util.*;
import java.util.regex.Pattern;

/**
 * Created by markusmueller on 20.04.20.
 */
public class RegExpProteinGrouperTest {

    @Test
    public void testApply() {
        List<Pattern> regex = new ArrayList<>();
        regex.add(Pattern.compile("^sp"));
        List<String> groups = new ArrayList<>();
        groups.add("canonical");
        groups.add("cryptic");

        ProteinModifGrouper psmGrouper = new ProteinModifGrouper(regex,groups);

        TObjectDoubleMap<String> scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("xcorr", NewAnceParams.getInstance().getMinXCorr());

        List<String> prots = new ArrayList<>();
        prots.add("sp|protein1");
        prots.add("protein2");

        PeptideSpectrumMatch peptideSpectrumMatch = new PeptideSpectrumMatch("spectrumFile", Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, false, false, null);

        Assert.assertEquals("canonical",psmGrouper.apply("",peptideSpectrumMatch));

        prots = new ArrayList<>();
        prots.add("protein2");

        peptideSpectrumMatch = new PeptideSpectrumMatch("spectrumFile", Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, false, false, null);


        Assert.assertEquals("cryptic",psmGrouper.apply("",peptideSpectrumMatch));
    }

    @Test
    public void testApply2() {
        List<Pattern> regex = new ArrayList<>();
        regex.add(Pattern.compile("^sp"));
        regex.add(Pattern.compile("ENST"));
        List<String> groups = new ArrayList<>();
        groups.add("canonical");
        groups.add("lnc");
        groups.add("ere");

        ProteinModifGrouper psmGrouper = new ProteinModifGrouper(regex,groups);

        TObjectDoubleMap<String> scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("xcorr", NewAnceParams.getInstance().getMinXCorr());

        List<String> prots = new ArrayList<>();
        prots.add("ENSTprotein1");
        prots.add("ENSTprotein2");

        PeptideSpectrumMatch peptideSpectrumMatch = new PeptideSpectrumMatch("spectrumFile", Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, false, false, null);

        Assert.assertEquals("lnc",psmGrouper.apply("",peptideSpectrumMatch));

        prots = new ArrayList<>();
        prots.add("protein2");

        peptideSpectrumMatch = new PeptideSpectrumMatch("spectrumFile", Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, false, false, null);

        Assert.assertEquals("ere",psmGrouper.apply("",peptideSpectrumMatch));

        prots = new ArrayList<>();
        prots.add("ENSTprotein1");
        prots.add("sp|protein1");

        peptideSpectrumMatch = new PeptideSpectrumMatch("spectrumFile", Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, false, false, null);

        Assert.assertEquals("canonical",psmGrouper.apply("",peptideSpectrumMatch));
    }

    @Test
    public void testApply3() {
        List<Pattern> regex = new ArrayList<>();
        regex.add(Pattern.compile("^sp"));
        regex.add(Pattern.compile("ENST"));
        List<String> groups = new ArrayList<>();
        groups.add("canonical");
        groups.add("lnc");
        groups.add("ere");

        Map<String,Set<String>> exclude = new HashMap<>();
        exclude.put("ere",new HashSet<>());
        exclude.get("ere").add("sp|protein1");

        RegExpProteinGrouper psmGrouper = new RegExpProteinGrouper(regex,groups,exclude);

        TObjectDoubleMap<String> scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("xcorr", NewAnceParams.getInstance().getMinXCorr());

        List<String> prots = new ArrayList<>();
        prots.add("sp|protein1");
        prots.add("ENSTprotein2");

        PeptideSpectrumMatch peptideSpectrumMatch = new PeptideSpectrumMatch("spectrumFile",
                Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,100, 101,
                1001.1, false, false, null);

        Assert.assertEquals("ere",psmGrouper.apply("",peptideSpectrumMatch));

        exclude.get("ere").add("ENSTprotein2");

        prots = new ArrayList<>();
        prots.add("ENSTprotein2");

        peptideSpectrumMatch = new PeptideSpectrumMatch("spectrumFile", Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, false, false, null);


        Assert.assertEquals("ere",psmGrouper.apply("",peptideSpectrumMatch));

        prots = new ArrayList<>();
        prots.add("protein2");

        peptideSpectrumMatch = new PeptideSpectrumMatch("spectrumFile", Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, false, false, null);


        Assert.assertEquals("ere",psmGrouper.apply("",peptideSpectrumMatch));

        prots = new ArrayList<>();
        prots.add("sp|protein2");

        peptideSpectrumMatch = new PeptideSpectrumMatch("spectrumFile", Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, false, false, null);


        Assert.assertEquals("canonical",psmGrouper.apply("",peptideSpectrumMatch));

        prots = new ArrayList<>();
        prots.add("ENSTprotein3");
        prots.add("protein3");

        peptideSpectrumMatch = new PeptideSpectrumMatch("spectrumFile", Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, false, false, null);


        Assert.assertEquals("lnc",psmGrouper.apply("",peptideSpectrumMatch));
    }
}
