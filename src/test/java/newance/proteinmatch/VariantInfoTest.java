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

        VariantInfo variantInfo = VariantInfo.parseFastaHeader(header);

        Assert.assertEquals("ENSP00000343864.2",variantInfo.getProteinID());

        Assert.assertEquals(5,variantInfo.getVariants().size());
        Assert.assertEquals(579,variantInfo.getVariants().get(0).getStartWT());
        Assert.assertEquals("Q",variantInfo.getVariants().get(0).getMutatedSequence());
        Assert.assertEquals("rs13302979_0",variantInfo.getVariants().get(0).getInfo());

        Assert.assertEquals(410,variantInfo.getVariants().get(1).getStartWT());
        Assert.assertEquals("Q",variantInfo.getVariants().get(1).getMutatedSequence());
        Assert.assertEquals("rs34903276_0",variantInfo.getVariants().get(1).getInfo());

        Assert.assertEquals(663,variantInfo.getVariants().get(4).getStartWT());
        Assert.assertEquals("A",variantInfo.getVariants().get(4).getMutatedSequence());
        Assert.assertEquals("rs7417106_0",variantInfo.getVariants().get(4).getInfo());
        Assert.assertEquals("ENSP00000343864.2 AA_SUBSTITUTION:579,579,Q,rs13302979_0 " +
                "AA_SUBSTITUTION:410,410,Q,rs34903276_0 " +
                "AA_SUBSTITUTION:406,406,E,rs13303368_0 " +
                "AA_SUBSTITUTION:398,398,G,rs13302983_0 " +
                "AA_SUBSTITUTION:663,663,A,rs7417106_0",variantInfo.toString());

        header = ">ENSP00000343864.2 \\Length=676";

        variantInfo = VariantInfo.parseFastaHeader(header);

        Assert.assertEquals(0, variantInfo.size());

        header = ">ENSP00000400001.1 \\Length=75  \\VariantComplex=(50|50|CVWWCGSLEAVRESKGGNYQLCTRKKI|rr189_0) (50|50|CVWWCGSLEAVRESKGGNYQLCTRKKI|rr189_0)";
        variantInfo = VariantInfo.parseFastaHeader(header);

        Assert.assertEquals("ENSP00000400001.1",variantInfo.getProteinID());

        Assert.assertEquals(2,variantInfo.getVariants().size());
        Assert.assertEquals(50,variantInfo.getVariants().get(0).getStartWT());
        Assert.assertEquals(50,variantInfo.getVariants().get(0).getEndWT());
        Assert.assertEquals("CVWWCGSLEAVRESKGGNYQLCTRKKI",variantInfo.getVariants().get(0).getMutatedSequence());
        Assert.assertEquals("rr189_0",variantInfo.getVariants().get(0).getInfo());

        Assert.assertEquals(50,variantInfo.getVariants().get(1).getStartWT());
        Assert.assertEquals(50,variantInfo.getVariants().get(1).getEndWT());
        Assert.assertEquals("CVWWCGSLEAVRESKGGNYQLCTRKKI",variantInfo.getVariants().get(1).getMutatedSequence());
        Assert.assertEquals("rr189_0",variantInfo.getVariants().get(1).getInfo());

        header = ">ENSP00000387671.1 \\Length=420 \\VariantSimple=(12|N|rs12026825_0) (283|V|rs2275253_0) (298|S|rs2275254_0) (321|L|rs36011905_0) (376|G|rs2256721_0)  \\VariantComplex=(2|2|TSPRLCTNRWPTLEPRLLLSSLVVNPP*|rs34698010_0)";

        variantInfo = VariantInfo.parseFastaHeader(header);

        Assert.assertEquals("ENSP00000387671.1",variantInfo.getProteinID());

        Assert.assertEquals(6,variantInfo.getVariants().size());

        Assert.assertEquals(12,variantInfo.getVariants().get(0).getStartWT());
        Assert.assertEquals(12,variantInfo.getVariants().get(0).getEndWT());
        Assert.assertEquals("N",variantInfo.getVariants().get(0).getMutatedSequence());
        Assert.assertEquals("rs12026825_0",variantInfo.getVariants().get(0).getInfo());
        Assert.assertEquals(SequenceVariant.VariantType.AA_SUBSTITUTION, variantInfo.getVariants().get(0).getType());

        Assert.assertEquals(376,variantInfo.getVariants().get(4).getStartWT());
        Assert.assertEquals(376,variantInfo.getVariants().get(4).getEndWT());
        Assert.assertEquals("G",variantInfo.getVariants().get(4).getMutatedSequence());
        Assert.assertEquals("rs2256721_0",variantInfo.getVariants().get(4).getInfo());
        Assert.assertEquals(SequenceVariant.VariantType.AA_SUBSTITUTION, variantInfo.getVariants().get(4).getType());

        Assert.assertEquals(2,variantInfo.getVariants().get(5).getStartWT());
        Assert.assertEquals(2,variantInfo.getVariants().get(5).getEndWT());
        Assert.assertEquals("TSPRLCTNRWPTLEPRLLLSSLVVNPP",variantInfo.getVariants().get(5).getMutatedSequence());
        Assert.assertEquals("rs34698010_0",variantInfo.getVariants().get(5).getInfo());
        Assert.assertEquals(SequenceVariant.VariantType.STOP, variantInfo.getVariants().get(5).getType());

        header = ">ENST00000387461.2.16-66";

        variantInfo = VariantInfo.parseFastaHeader(header);

        Assert.assertEquals("ENST00000387461.2.16-66",variantInfo.getProteinID());

        Assert.assertEquals(0,variantInfo.getVariants().size());
    }


}
