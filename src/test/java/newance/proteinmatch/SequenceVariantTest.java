package newance.proteinmatch;

import org.junit.Assert;
import org.junit.Test;

/**
 * Created by markusmueller on 18.06.20.
 */
public class SequenceVariantTest {

    @Test
    public void test_VariantSimple() {

        String mutation = "(141|A|rs75062661_0)";
        SequenceVariant sequenceVariant = SequenceVariant.parseSimpleVariantString("ENSP00000334393.3", 305, mutation);

        Assert.assertEquals(141,sequenceVariant.getStartWT());
        Assert.assertEquals(141,sequenceVariant.getEndWT());
        Assert.assertEquals("A",sequenceVariant.getMutatedSequence());
        Assert.assertEquals(1,sequenceVariant.getLength());
        Assert.assertEquals("rs75062661_0",sequenceVariant.getInfo());
        Assert.assertEquals(SequenceVariant.VariantType.AA_SUBSTITUTION,sequenceVariant.getType());

        mutation = "(24|*|rs6671527_0)";
        sequenceVariant = SequenceVariant.parseSimpleVariantString("ENSP00000271139.8", 268, mutation);

        Assert.assertEquals(24,sequenceVariant.getStartWT());
        Assert.assertEquals(24,sequenceVariant.getEndWT());
        Assert.assertEquals("",sequenceVariant.getMutatedSequence());
        Assert.assertEquals(0,sequenceVariant.getLength());
        Assert.assertEquals("rs6671527_0",sequenceVariant.getInfo());
        Assert.assertEquals(SequenceVariant.VariantType.STOP,sequenceVariant.getType());

        mutation = "(205|R|rs4970490_0)";
        sequenceVariant = SequenceVariant.parseSimpleVariantString("ENSP00000363278.1",204, mutation);

        Assert.assertEquals(205,sequenceVariant.getStartWT());
        Assert.assertEquals(205,sequenceVariant.getEndWT());
        Assert.assertEquals("R",sequenceVariant.getMutatedSequence());
        Assert.assertEquals(1,sequenceVariant.getLength());
        Assert.assertEquals("rs4970490_0",sequenceVariant.getInfo());
        Assert.assertEquals(SequenceVariant.VariantType.INSERTION,sequenceVariant.getType());
    }

    @Test
    public void test_VariantComplex() {

        String mutation = "(117|119|GL|rs10578519_3)";
        SequenceVariant sequenceVariant = SequenceVariant.parseComplexVariantString("ENSP00000483131.1", 387, mutation);

        Assert.assertEquals(117,sequenceVariant.getStartWT());
        Assert.assertEquals(119,sequenceVariant.getEndWT());
        Assert.assertEquals("GL",sequenceVariant.getMutatedSequence());
        Assert.assertEquals(2,sequenceVariant.getLength());
        Assert.assertEquals("rs10578519_3",sequenceVariant.getInfo());
        Assert.assertEquals(SequenceVariant.VariantType.DELETION,sequenceVariant.getType());

        mutation = "(117|119||rs10578519_3)";
        sequenceVariant = SequenceVariant.parseComplexVariantString("ENSP00000483131.1", 387, mutation);

        Assert.assertEquals(117,sequenceVariant.getStartWT());
        Assert.assertEquals(119,sequenceVariant.getEndWT());
        Assert.assertEquals("",sequenceVariant.getMutatedSequence());
        Assert.assertEquals(0,sequenceVariant.getLength());
        Assert.assertEquals("rs10578519_3",sequenceVariant.getInfo());
        Assert.assertEquals(SequenceVariant.VariantType.DELETION,sequenceVariant.getType());

        mutation = "(117|119|GLA|rs10578519_3)";
        sequenceVariant = SequenceVariant.parseComplexVariantString("ENSP00000483131.1", 387, mutation);

        Assert.assertEquals(117,sequenceVariant.getStartWT());
        Assert.assertEquals(119,sequenceVariant.getEndWT());
        Assert.assertEquals("GLA",sequenceVariant.getMutatedSequence());
        Assert.assertEquals(3,sequenceVariant.getLength());
        Assert.assertEquals("rs10578519_3",sequenceVariant.getInfo());
        Assert.assertEquals(SequenceVariant.VariantType.SUBSTITUTION,sequenceVariant.getType());

        mutation = "(4|4|KRK|rs144636354_0)";
        sequenceVariant = SequenceVariant.parseComplexVariantString("ENST00000598846.1.3103-3153", 17, mutation);

        Assert.assertEquals(4,sequenceVariant.getStartWT());
        Assert.assertEquals(4,sequenceVariant.getEndWT());
        Assert.assertEquals("KRK",sequenceVariant.getMutatedSequence());
        Assert.assertEquals(3,sequenceVariant.getLength());
        Assert.assertEquals("rs144636354_0",sequenceVariant.getInfo());
        Assert.assertEquals(SequenceVariant.VariantType.INSERTION,sequenceVariant.getType());

        mutation = "(342|342|QWHLSMRTPTAWAGGSRMRWPSATPRLPTRQRP*|rs10709483_0)";
        sequenceVariant = SequenceVariant.parseComplexVariantString("ENSP00000483125.1", 694, mutation);

        Assert.assertEquals(342,sequenceVariant.getStartWT());
        Assert.assertEquals(342,sequenceVariant.getEndWT());
        Assert.assertEquals("QWHLSMRTPTAWAGGSRMRWPSATPRLPTRQRP",sequenceVariant.getMutatedSequence());
        Assert.assertEquals("QWHLSMRTPTAWAGGSRMRWPSATPRLPTRQRP".length(),sequenceVariant.getLength());
        Assert.assertEquals("rs10709483_0",sequenceVariant.getInfo());
        Assert.assertEquals(SequenceVariant.VariantType.STOP,sequenceVariant.getType());

    }

    @Test
    public void test_match() {

        String mutation = "(141|A|rs75062661_0)";
        SequenceVariant sequenceVariant = SequenceVariant.parseSimpleVariantString("ENSP00000334393.3", 305, mutation);

        Assert.assertTrue(sequenceVariant.match("ACDCACDC", 4));
        Assert.assertEquals(5, sequenceVariant.getPosAfterVariant());
        Assert.assertTrue(sequenceVariant.match("ACDCACDC", 0));
        Assert.assertEquals(1, sequenceVariant.getPosAfterVariant());
        Assert.assertFalse(sequenceVariant.match("ACDCCDC", 4));
        Assert.assertEquals(-1, sequenceVariant.getPosAfterVariant());

        mutation = "(117|119|GL|rs10578519_3)";
        sequenceVariant = SequenceVariant.parseComplexVariantString("ENSP00000483131.1", 387, mutation);

        Assert.assertTrue(sequenceVariant.match("ACDCGLACDC", 4));
        Assert.assertEquals(6, sequenceVariant.getPosAfterVariant());
        Assert.assertTrue(sequenceVariant.match("LACDC", 0));
        Assert.assertEquals(1, sequenceVariant.getPosAfterVariant());
        Assert.assertFalse(sequenceVariant.match("ACDCGCDC", 4));
        Assert.assertEquals(-1, sequenceVariant.getPosAfterVariant());
        Assert.assertFalse(sequenceVariant.match("ACDCLCDC", 4));
        Assert.assertEquals(-1, sequenceVariant.getPosAfterVariant());

        mutation = "(342|342|QWHLSMRTPTAWAGGSRMRWPSATPRLPTRQRP*|rs10709483_0)";
        sequenceVariant = SequenceVariant.parseComplexVariantString("ENSP00000483125.1", 694, mutation);

        Assert.assertTrue(sequenceVariant.match("QWHLSMRTPTAWAGGSRMRWPSATPRLPTRQRP", 0));
        Assert.assertEquals(33, sequenceVariant.getPosAfterVariant());
        Assert.assertTrue(sequenceVariant.match("QWHLSMRTPTAWAGGSRMRWPSATPRL", 0));
        Assert.assertEquals(27, sequenceVariant.getPosAfterVariant());
        Assert.assertTrue(sequenceVariant.match("SMRTPTAWAGGSRMRWPSAT", 0));
        Assert.assertEquals(20, sequenceVariant.getPosAfterVariant());
        Assert.assertTrue(sequenceVariant.match("ABCDQWHLSMRTPTAWAGGSRMRWPSAT", 4));
        Assert.assertEquals(28, sequenceVariant.getPosAfterVariant());
        Assert.assertFalse(sequenceVariant.match("ABCDQWHLSMRTPTAWAGGSRMRWPSAT", 2));
        Assert.assertEquals(-1, sequenceVariant.getPosAfterVariant());
        Assert.assertFalse(sequenceVariant.match("ABCDQWHLSMRTPTAWAGGSRMRWPSAT", 0));
        Assert.assertEquals(-1, sequenceVariant.getPosAfterVariant());
        Assert.assertFalse(sequenceVariant.match("QWHLSMRTPTAWAGGSRMRWPSATPRLPTRQRPABC", 0));
        Assert.assertEquals(-1, sequenceVariant.getPosAfterVariant());

    }
}
