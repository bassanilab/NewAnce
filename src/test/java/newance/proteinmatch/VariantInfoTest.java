package newance.proteinmatch;

import org.junit.Assert;
import org.junit.Test;

/**
 * Created by markusmueller on 13.05.20.
 */
public class VariantInfoTest {

    @Test
    public void testConstructor() {
        String header = ">ENSP00000343864.2 \\Length=676 \\VariantSimple=(579|Q|rs13302979_0) (410|Q|rs34903276_0) (406|E|rs13303368_0) (398|G|rs13302983_0) (663|A|rs7417106_0)";

        VariantInfo variantInfo = new VariantInfo(header);

        Assert.assertEquals("ENSP00000343864.2",variantInfo.getProteinID());

        Assert.assertEquals(5,variantInfo.getVariants().size());
        Assert.assertEquals(579,variantInfo.getVariants().get(0).getPosition());
        Assert.assertEquals('Q',variantInfo.getVariants().get(0).getMutatedAA());
        Assert.assertEquals("rs13302979_0",variantInfo.getVariants().get(0).getInfo());

        Assert.assertEquals(410,variantInfo.getVariants().get(1).getPosition());
        Assert.assertEquals('Q',variantInfo.getVariants().get(1).getMutatedAA());
        Assert.assertEquals("rs34903276_0",variantInfo.getVariants().get(1).getInfo());

        Assert.assertEquals(663,variantInfo.getVariants().get(4).getPosition());
        Assert.assertEquals('A',variantInfo.getVariants().get(4).getMutatedAA());
        Assert.assertEquals("rs7417106_0",variantInfo.getVariants().get(4).getInfo());
        Assert.assertEquals("ENSP00000343864.2 \\VariantSimple=(579|Q|rs13302979_0) (410|Q|rs34903276_0) " +
                "(406|E|rs13303368_0) (398|G|rs13302983_0) (663|A|rs7417106_0)",variantInfo.toString());

        header = ">ENSP00000343864.2 \\Length=676";

        variantInfo = new VariantInfo(header);

        Assert.assertEquals("ENSP00000343864.2",variantInfo.getProteinID());
        Assert.assertEquals(0,variantInfo.getVariants().size());
        Assert.assertEquals("ENSP00000343864.2",variantInfo.toString());
    }

    @Test
    public void test_getVariantInfo() {

        String header = ">ENSP00000343864.2 \\Length=676 \\VariantSimple=(579|Q|rs13302979_0) (410|Q|rs34903276_0) (406|E|rs13303368_0) (398|G|rs13302983_0) (663|A|rs7417106_0)";

        VariantInfo variantInfo = new VariantInfo(header);

        Assert.assertEquals("rs13302979_0", variantInfo.getVariantInfo(570, 'Q', 8));
        Assert.assertEquals("rs34903276_0", variantInfo.getVariantInfo(405, 'Q', 4));
        Assert.assertEquals("rs7417106_0", variantInfo.getVariantInfo(660, 'A', 2));
        Assert.assertEquals("", variantInfo.getVariantInfo(570, 'P', 8));
        Assert.assertEquals("", variantInfo.getVariantInfo(570, 'Q', 9));
        Assert.assertEquals("", variantInfo.getVariantInfo(570, 'Q', 7));
    }
}
