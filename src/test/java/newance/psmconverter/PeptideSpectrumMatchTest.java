package newance.psmconverter;

import gnu.trove.map.hash.TObjectDoubleHashMap;
import newance.mzjava.mol.Peptide;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

/**
 * Created by markusmueller on 14.05.20.
 */
public class PeptideSpectrumMatchTest {

    @Test
    public void test_getWTSequence() {

        Peptide peptide = Peptide.parse("QSEDGSHTIQIMY");

        PeptideSpectrumMatch psm = new PeptideSpectrumMatch("Spectrum_file", peptide, new HashSet<>(),
                new TObjectDoubleHashMap<>(), 2, 1, 100, 100, 1507.661308,
                false, false, null, null);

        Assert.assertEquals("QSEDGSHTIQIMY",psm.getWTSequence());

        List<Integer> varPos = new ArrayList<>();
        List<Character> wtAAs = new ArrayList<>();

        varPos.add(3);
        wtAAs.add('A');

        psm = new PeptideSpectrumMatch("Spectrum_file", peptide, new HashSet<>(),
                new TObjectDoubleHashMap<>(), 2, 1, 100, 100, 1507.661308,
                false, true, varPos, wtAAs);

        Assert.assertEquals("QSEAGSHTIQIMY",psm.getWTSequence());

        varPos.add(6);
        wtAAs.add('D');
        psm = new PeptideSpectrumMatch("Spectrum_file", peptide, new HashSet<>(),
                new TObjectDoubleHashMap<>(), 2, 1, 100, 100, 1507.661308,
                false, true, varPos, wtAAs);

        Assert.assertEquals("QSEAGSDTIQIMY",psm.getWTSequence());
    }
}
