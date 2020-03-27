package newance.psmconverter;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
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

        int endIdx = maxQuantPsmReader.parseModification(7,6, modif.toCharArray(), modMatchMap);

        Assert.assertEquals(11,endIdx);

        Assert.assertTrue(modMatchMap.containsKey(6));
        Assert.assertEquals("Phospho:HO3P",modMatchMap.get(6).get(0).toString());

        endIdx = maxQuantPsmReader.parseModification(16,11, modif.toCharArray(), modMatchMap);

        Assert.assertEquals(20,endIdx);

        Assert.assertTrue(modMatchMap.containsKey(11));
        Assert.assertEquals("Oxidation:O",modMatchMap.get(11).get(0).toString());

        modif = "_LGALM(Oxidation (M))Y_";

        modMatchMap = ArrayListMultimap.create();

        endIdx = maxQuantPsmReader.parseModification(6,5, modif.toCharArray(), modMatchMap);

        Assert.assertEquals(21,endIdx);

        Assert.assertTrue(modMatchMap.containsKey(5));
        Assert.assertEquals("Oxidation:O",modMatchMap.get(5).get(0).toString());

        modif = "_LGALS(Phospho (STY))Y_";

        modMatchMap = ArrayListMultimap.create();

        endIdx = maxQuantPsmReader.parseModification(6,5, modif.toCharArray(), modMatchMap);

        Assert.assertEquals(21,endIdx);

        Assert.assertTrue(modMatchMap.containsKey(5));
        Assert.assertEquals("Phospho:HO3P",modMatchMap.get(5).get(0).toString());
    }

    @Test
    public void testMakeModifiedPeptideMatch() throws SQLException{

        MaxQuantPsmReader maxQuantPsmReader = new MaxQuantPsmReader();

        String modif = "_AVDWWS(ph)LGALM(ox)Y_";

        PeptideMatchDataWrapper peptideMatchDataWrapper = maxQuantPsmReader.makeModifiedPeptideMatch(modif);

        Assert.assertEquals(2,peptideMatchDataWrapper.getModificationCount());
        Assert.assertEquals("AVDWWSL(Phospho)GALMY(Oxidation)",peptideMatchDataWrapper.toPeptide().toString());

        modif = "_AVDWWS(Phospho (STY))LGALM(Oxidation (M))Y_";

        peptideMatchDataWrapper = maxQuantPsmReader.makeModifiedPeptideMatch(modif);

        Assert.assertEquals(2,peptideMatchDataWrapper.getModificationCount());
        Assert.assertEquals("AVDWWSL(Phospho)GALMY(Oxidation)",peptideMatchDataWrapper.toPeptide().toString());
    }
}
