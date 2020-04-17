package newance.psmconverter;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import newance.mzjava.mol.modification.ModAttachment;
import newance.mzjava.mol.modification.Modification;
import org.junit.Assert;
import org.junit.Test;

import java.sql.SQLException;

/**
 * Created by markusmueller on 27.03.20.
 */
public class MaxQuantReaderTest {

    @Test
    public void testParseModif() {

        MaxQuantPsmReader maxQuantPsmReader = new MaxQuantPsmReader();

        String modif = "_AVDWWS(ph)LGALM(ox)Y_";

        ListMultimap<Object, Modification> modMatchMap = ArrayListMultimap.create();

        int endIdx = maxQuantPsmReader.parseModification(7,5, modif.toCharArray(), modMatchMap);

        Assert.assertEquals(11,endIdx);

        Assert.assertTrue(modMatchMap.containsKey(5));
        Assert.assertEquals("Phospho:HO3P",modMatchMap.get(5).get(0).toString());

        endIdx = maxQuantPsmReader.parseModification(16,10, modif.toCharArray(), modMatchMap);

        Assert.assertEquals(20,endIdx);

        Assert.assertTrue(modMatchMap.containsKey(10));
        Assert.assertEquals("Oxidation:O",modMatchMap.get(10).get(0).toString());

        modif = "_LGALM(Oxidation (M))Y_";

        modMatchMap = ArrayListMultimap.create();

        endIdx = maxQuantPsmReader.parseModification(6,4, modif.toCharArray(), modMatchMap);

        Assert.assertEquals(21,endIdx);

        Assert.assertTrue(modMatchMap.containsKey(4));
        Assert.assertEquals("Oxidation:O",modMatchMap.get(4).get(0).toString());

        modif = "_LGALS(Phospho (STY))Y_";

        modMatchMap = ArrayListMultimap.create();

        endIdx = maxQuantPsmReader.parseModification(6,4, modif.toCharArray(), modMatchMap);

        Assert.assertEquals(21,endIdx);

        Assert.assertTrue(modMatchMap.containsKey(4));
        Assert.assertEquals("Phospho:HO3P",modMatchMap.get(4).get(0).toString());

        modif = "_(Acetyl (Protein N-term))ADISLDELIRKRGAAA_";

        modMatchMap = ArrayListMultimap.create();

        endIdx = maxQuantPsmReader.parseModification(1,-1, modif.toCharArray(), modMatchMap);

        Assert.assertEquals(26,endIdx);

        Assert.assertTrue(modMatchMap.containsKey(ModAttachment.N_TERM));
        Assert.assertEquals("Acetyl:C2H2O",modMatchMap.get(ModAttachment.N_TERM).get(0).toString());
    }

    @Test
    public void testMakeModifiedPeptideMatch() throws SQLException{

        MaxQuantPsmReader maxQuantPsmReader = new MaxQuantPsmReader();

        String modif = "_AVDWWS(ph)LGALM(ox)Y_";

        PeptideMatchDataWrapper peptideMatchDataWrapper = maxQuantPsmReader.makeModifiedPeptideMatch(modif);

        Assert.assertEquals(2,peptideMatchDataWrapper.getModificationCount());
        Assert.assertEquals("AVDWWS(Phospho)LGALM(Oxidation)Y",peptideMatchDataWrapper.toPeptide().toString());

        modif = "_AVDWWS(Phospho (STY))LGALM(Oxidation (M))Y_";

        peptideMatchDataWrapper = maxQuantPsmReader.makeModifiedPeptideMatch(modif);

        Assert.assertEquals(2,peptideMatchDataWrapper.getModificationCount());
        Assert.assertEquals("AVDWWS(Phospho)LGALM(Oxidation)Y",peptideMatchDataWrapper.toPeptide().toString());
    }
}
