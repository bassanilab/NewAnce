package newance.proteinmatch;

import org.junit.Assert;
import org.junit.Test;

import java.util.List;

/**
 * Created by markusmueller on 03.07.20.
 */
public class VariantProteinTest {

    private String fastaStr = ">ENSP00000368717.2 \\Length=576 \\VariantSimple=(452|P|rs3829740_2)\n" +
            "MGNSHCVPQAPRRLRASFSRKPSLKGNREDSARMSAGLPGPEAARSGDAAANKLFHYIPG\n" +
            "TDILDLENQRENLEQPFLSVFKKGRRRVPVRNLGKVVHYAKVQLRFQHSQDVSDCYLELF\n" +
            "PAHLYFQAHGSEGLTFQGLLPLTELSVCPLEGSREHAFQITGPLPAPLLVLCPSRAELDR\n" +
            "WLYHLEKQTALLGGPRRCHSAPPQGSCGDELPWTLQRRLTRLRTASGHEPGGSAVCASRV\n" +
            "KLQHLPAQEQWDRLLVLYPTSLAIFSEELDGLCFKGELPLRAVHINLEEKEKQIRSFLIE\n" +
            "GPLINTIRVVCASYEDYGHWLLCLRAVTHREGAPPLPGAESFPGSQVMGSGRGSLSSGGQ\n" +
            "TSWDSGCLAPPSTRTSHSLPESSVPSTVGCSSQHTPLHRLSLESSPDAPDHTSETSHSPL\n" +
            "YADPYTPPATSHRRVTDVRGLEEFLSAMQSARGPTPSSPLPSVPVSVPASDPRSCSSGPA\n" +
            "GPYLLSKKGALQSRAAQRHRGSAKDGGPQPPDAPQLVSSAREGSPEPWLPLTDGRSPRRS\n" +
            "RDPGYDHLWDETLSSSHQKCPQLGGPEASGGLVQWI\n" +
            ">ENSP00000462558.1 \\Length=176 \\VariantSimple=(52|P|rs3829740_3)\n" +
            "SLESSPDAPDHTSETSHSPLYADPYTPPATSHRRVTDVRGLEEFLSAMQSARGPTPSSPL\n" +
            "PSVPVSVPASDPRSCSSGPAGPYLLSKKGALQSRAAQRHRGSAKDGGPQPPDAPQLVSSA\n" +
            "REGSPEPWLPLTDGRSPRRSRDPGYDHLWDETLSSSHQKCPQLGGPEASGGLVQWI\n" +
            ">ENSP00000343864.2 \\Length=676 \\VariantSimple=(579|Q|rs13302979_0) (410|Q|rs34903276_0) (406|E|rs13303368_0) (398|G|rs13302983_0) (663|A|rs7417106_0)\n" +
            "MEPRGGGSSQFSSCPGPASSGDQMQRLLQGPAPRPPGEPPGSPKSPGHSTGSQRPPDSPG\n" +
            "APPRSPSRKKRRAVGAKGGGHTGASASAQTGSPLLPAASPETAKLMAKAGQEELGPGPAG\n" +
            "APEPGPRSPVQEDRPGPGLGLSTPVPVTEQGTDQIRTPRRAKLHTVSTTVWEALPDVSRA\n" +
            "KSDMAVSTPASEPQPDRDMAVSTPASEPQSDRDMAVSTPASEPQPDTDMAVSTPASEPQP\n" +
            "DRDMAVSIPASKPQSDTAVSTPASEPQSSVALSTPISKPQLDTDVAVSTPASKHGLDVAL\n" +
            "PTAGPVAKLEVASSPPVSEAVPRMTESSGLVSTPVPRADAAGLAWPPTRRAGPDVVEMEA\n" +
            "VVSEPSAGAPGCCSGAPALGLTQVPRKKKVRFSVAGPSPNKPGSGQASARPSAPQTATGA\n" +
            "HGGPGAWEAVAVGPRPHQPRILKHLPRPPPSAVTRVGPGSSFAVTLPEAYEFFFCDTIEE\n" +
            "NEEAEAAAAGQDPAGVQWPDMCEFFFPDVGAQRSRRRGSPEPLPRADPVPAPIPGDPVPI\n" +
            "SIPEVYEHFFFGEDRLEGVLGPAVPLPLQALEPPRSASEGAGPGTPLKPAVVERLHLALR\n" +
            "RAGELRGPVPSFAFSQNDMCLVFVAFATWAVRTSDPHTPDAWKTALLANVGTISAIRYFR\n" +
            "RQVGQGRRSHSPSPSS\n" +
            ">ENSP00000414022.3 \\Length=790 \\VariantSimple=(693|Q|rs13302979_1) (524|Q|rs34903276_1) (520|E|rs13303368_1) (512|G|rs13302983_1) (777|A|rs7417106_1)\n" +
            "MENFQYSVQLSDQDWAEFSATADECGLLQAGLASGDELLSSDIDQGDSSGSSPPRAPPLP\n" +
            "TGQLAAGGRSRRGCEEEDVATQQPVSRSQGEPVLALGTGQQTPSTSARAEAPPSLGPGAS\n" +
            "PPSQFSSCPGPASSGDQMQRLLQGPAPRPPGEPPGSPKSPGHSTGSQRPPDSPGAPPRSP\n" +
            "SRKKRRAVGAKGGGHTGASASAQTGSPLLPAASPETAKLMAKAGQEELGPGPAGAPEPGP\n" +
            "RSPVQEDRPGPGLGLSTPVPVTEQGTDQIRTPRRAKLHTVSTTVWEALPDVSRAKSDMAV\n" +
            "STPASEPQPDRDMAVSTPASEPQSDRDMAVSTPASEPQPDTDMAVSTPASEPQPDRDMAV\n" +
            "SIPASKPQSDTAVSTPASEPQSSVALSTPISKPQLDTDVAVSTPASKHGLDVALPTAGPV\n" +
            "AKLEVASSPPVSEAVPRMTESSGLVSTPVPRADAAGLAWPPTRRAGPDVVEMEAVVSEPS\n" +
            "AGAPGCCSGAPALGLTQVPRKKKVRFSVAGPSPNKPGSGQASARPSAPQTATGAHGGPGA\n" +
            "WEAVAVGPRPHQPRILKHLPRPPPSAVTRVGPGSSFAVTLPEAYEFFFCDTIEENEEAEA\n" +
            "AAAGQDPAGVQWPDMCEFFFPDVGAQRSRRRGSPEPLPRADPVPAPIPGDPVPISIPEVY\n" +
            "EHFFFGEDRLEGVLGPAVPLPLQALEPPRSASEGAGPGTPLKPAVVERLHLALRRAGELR\n" +
            "GPVPSFAFSQNDMCLVFVAFATWAVRTSDPHTPDAWKTALLANVGTISAIRYFRRQVGQG\n" +
            "RRSHSPSPSS\n" +
            ">ENSP00000393198.2 \\Length=247 \\VariantSimple=(44|S|rs2298214_0)\n" +
            "MAADTPGKPSASPMAGAPASASRTPDKPRSAAEHRKVGSRPGVRGATGGREGRGTQPVPD\n" +
            "PQSSKPVMEKRRRARINESLAQLKTLILDALRKESSRHSKLEKADILEMTVRHLRSLRRV\n" +
            "QVTAALSADPAVLGKYRAGFHECLAEVNRFLAGCEGVPADVRSRLLGHLAACLRQLGPSR\n" +
            "RPASLSPAAPAEAPAPEVYAGRPLLPSLGGPFPLLAPPLLPGLTRALPAAPRAGPQGPGG\n" +
            "PWRPWLR\n" +
            ">ENSP00000304595.6 \\Length=221\n" +
            "MAADTPGKPSASPMAGAPASASRTPDKPRSAAEHRKSSKPVMEKRRRARINESLAQLKTL\n" +
            "ILDALRKESSRHSKLEKADILEMTVRHLRSLRRVQVTAALSADPAVLGKYRAGFHECLAE\n" +
            "VNRFLAGCEGVPADVRSRLLGHLAACLRQLGPSRRPASLSPAAPAEAPAPEVYAGRPLLP\n" +
            "SLGGPFPLLAPPLLPGLTRALPAAPRAGPQGPGGPWRPWLR\n"+
            ">ENSP00000379873.1 \\Length=365 \\VariantSimple=(10|V|rs1143146_0) (33|S|rs2075684_0) (68|K|rs707910_0) (86|E|rs1059455_0) (89|G|rs199474430_0) (90|K|rs199474436_0) (91|M|rs79361534_0) (94|H|rs78306866_0) (100|A|rs1071742_0) (101|N|rs1136688_0) (103|R|rs1136689_0) (114|D|rs1136692_0) (119|L|rs1071743_0) (121|M|rs1136695_0) (123|F|rs1136697_0) (129|P|rs1136700_0) (138|Q|rs3173420_0) (140|Y|rs3173419_0) (151|K|rs1059509_0) (168|Q|rs1059517_0) (174|V|rs1059535_0) (175|R|rs1059536_0) (176|A|rs9256983_0) (180|W|rs9260156_0) (182|V|rs9260157_0) (185|E|rs1059542_0) (187|P|rr1730_0) (187|R|rs3129018_0) (190|D|rs3098019_0) (191|G|rs76185201_0) (251|H|rs145046067_0) (300|P|rs1136903_0) (306|V|rs1136949_0) (307|H|rs41558618_0) (335|N|rs1137160_0) (345|S|rs2231119_0)  \\VariantComplex=(103|105|GM|rr1725_0) (105|105|LAATTTRARPVLTPSR*|rr1726_0) (107|107|APLLQPERGRFSHHPDNVWLRRGVGRALPPRVPAGRLRRQGLHRPERGPALLDRGGHGGSDHQAQVGGGP*|rr1727_0) (250|250|QTRSSWRPGLQGMEPSRSGRLWWCLLERSRDTPAMCSMRVCPSPSP*|rs45576436_0)\n" +
            "MAVMAPRTLLLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRF\n" +
            "DSDAASQRMEPRAPWIEQEGPEYWDQETRNVKAQSQTDRVDLGTLRGYYNQSEAGSHTIQ\n" +
            "IMYGCDVGSDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAAHEAEQL\n" +
            "RAYLDGTCVEWLRRYLENGKETLQRTDPPKTHMTHHPISDHEATLRCWALGFYPAEITLT\n" +
            "WQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWEL\n" +
            "SSQPTIPIVGIIAGLVLLGAVITGAVVAAVMWRRKSSDRKGGSYTQAASSDSAQGSDVSL\n" +
            "TACKV\n" +
            ">ENSP00000366002.5 \\Length=371 \\VariantSimple=(10|V|rs1143146_0) (33|S|rs2075684_0) (68|K|rs707910_0) (86|E|rs1059455_0) (89|G|rs199474430_0) (90|K|rs199474436_0) (91|M|rs79361534_0) (94|H|rs78306866_0) (100|A|rs1071742_0) (101|N|rs1136688_0) (103|R|rs1136689_0) (114|D|rs1136692_0) (119|L|rs1071743_0) (121|M|rs1136695_0) (123|F|rs1136697_0) (129|P|rs1136700_0) (138|Q|rs3173420_0) (140|Y|rs3173419_0) (151|K|rs1059509_0) (168|Q|rs1059517_0) (174|V|rs1059535_0) (175|R|rs1059536_0) (176|A|rs9256983_0) (180|W|rs9260156_0) (182|V|rs9260157_0) (185|E|rs1059542_0) (187|P|rr1730_0) (187|R|rs3129018_0) (190|D|rs3098019_0) (191|G|rs76185201_0) (251|H|rs145046067_0) (300|P|rs1136903_0) (306|V|rs1136949_0) (307|H|rs41558618_0) (335|N|rs1137160_0) (351|S|rs2231119_1)  \\VariantComplex=(103|105|GM|rr1725_0) (105|105|LAATTTRARPVLTPSR*|rr1726_0) (107|107|APLLQPERGRFSHHPDNVWLRRGVGRALPPRVPAGRLRRQGLHRPERGPALLDRGGHGGSDHQAQVGGGP*|rr1727_0) (250|250|QTRSSWRPGLQGMEPSRSGRLWWCLLERSRDTPAMCSMRVCPSPSP*|rs45576436_0)\n" +
            "MAVMAPRTLLLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRF\n" +
            "DSDAASQRMEPRAPWIEQEGPEYWDQETRNVKAQSQTDRVDLGTLRGYYNQSEAGSHTIQ\n" +
            "IMYGCDVGSDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAAHEAEQL\n" +
            "RAYLDGTCVEWLRRYLENGKETLQRTDPPKTHMTHHPISDHEATLRCWALGFYPAEITLT\n" +
            "WQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWEL\n" +
            "SSQPTIPIVGIIAGLVLLGAVITGAVVAAVMWRRKSSGGEGVKDRKGGSYTQAASSDSAQ\n" +
            "GSDVSLTACKV\n" +
            ">ENSP00000366005.5 \\Length=365 \\VariantSimple=(10|V|rs1143146_0) (33|S|rs2075684_0) (68|K|rs707910_0) (86|E|rs1059455_0) (89|G|rs199474430_0) (90|K|rs199474436_0) (91|M|rs79361534_0) (94|H|rs78306866_0) (100|A|rs1071742_0) (101|N|rs1136688_0) (103|R|rs1136689_0) (114|D|rs1136692_0) (119|L|rs1071743_0) (121|M|rs1136695_0) (123|F|rs1136697_0) (129|P|rs1136700_0) (138|Q|rs3173420_0) (140|Y|rs3173419_0) (151|K|rs1059509_0) (168|Q|rs1059517_0) (174|V|rs1059535_0) (175|R|rs1059536_0) (176|A|rs9256983_0) (180|W|rs9260156_0) (182|V|rs9260157_0) (185|E|rs1059542_0) (187|P|rr1730_0) (187|R|rs3129018_0) (190|D|rs3098019_0) (191|G|rs76185201_0) (251|H|rs145046067_0) (300|P|rs1136903_0) (306|V|rs1136949_0) (307|H|rs41558618_0) (335|N|rs1137160_0) (345|S|rs2231119_0)  \\VariantComplex=(103|105|GM|rr1725_0) (105|105|LAATTTRARPVLTPSR*|rr1726_0) (107|107|APLLQPERGRFSHHPDNVWLRRGVGRALPPRVPAGRLRRQGLHRPERGPALLDRGGHGGSDHQAQVGGGP*|rr1727_0) (250|250|QTRSSWRPGLQGMEPSRSGRLWWCLLERSRDTPAMCSMRVCPSPSP*|rs45576436_0)\n" +
            "MAVMAPRTLLLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRF\n" +
            "DSDAASQRMEPRAPWIEQEGPEYWDQETRNVKAQSQTDRVDLGTLRGYYNQSEAGSHTIQ\n" +
            "IMYGCDVGSDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAAHEAEQL\n" +
            "RAYLDGTCVEWLRRYLENGKETLQRTDPPKTHMTHHPISDHEATLRCWALGFYPAEITLT\n" +
            "WQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWEL\n" +
            "SSQPTIPIVGIIAGLVLLGAVITGAVVAAVMWRRKSSDRKGGSYTQAASSDSAQGSDVSL\n" +
            "TACKV\n" +
            ">ENSP00000365998.2 \\Length=299 \\VariantSimple=(10|V|rs1143146_0) (33|S|rs2075684_0) (68|K|rs707910_0) (86|E|rs1059455_0) (89|G|rs199474430_0) (90|K|rs199474436_0) (91|M|rs79361534_0) (94|H|rs78306866_0) (100|A|rs1071742_0) (101|N|rs1136688_0) (103|R|rs1136689_0) (114|D|rs1136692_0) (119|L|rs1071743_0) (121|M|rs1136695_0) (123|F|rs1136697_0) (129|P|rs1136700_0) (138|Q|rs3173420_0) (140|Y|rs3173419_0) (151|K|rs1059509_0) (168|Q|rs1059517_0) (174|V|rs1059535_0) (175|R|rs1059536_0) (176|A|rs9256983_0) (180|W|rs9260156_0) (182|V|rs9260157_0) (185|E|rs1059542_0) (187|P|rr1730_0) (187|R|rs3129018_0) (190|D|rs3098019_0) (191|G|rs76185201_0) (251|H|rs145046067_0)  \\VariantComplex=(103|105|GM|rr1725_0) (105|105|LAATTTRARPVLTPSR*|rr1726_0) (107|107|APLLQPERGRFSHHPDNVWLRRGVGRALPPRVPAGRLRRQGLHRPERGPALLDRGGHGGSDHQAQVGGGP*|rr1727_0) (250|250|QTRSSWRPGLQGMEPSRSGRLWWCLLERSRDTPAMCSMRVCPSPSP*|rs45576436_0)\n" +
            "MAVMAPRTLLLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRF\n" +
            "DSDAASQRMEPRAPWIEQEGPEYWDQETRNVKAQSQTDRVDLGTLRGYYNQSEAGSHTIQ\n" +
            "IMYGCDVGSDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAAHEAEQL\n" +
            "RAYLDGTCVEWLRRYLENGKETLQRTDPPKTHMTHHPISDHEATLRCWALGFYPAEITLT\n" +
            "WQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWV\n" +
            ">ENST00000495449.1.38-256 \\Length=73 \\VariantSimple=(61|F|rs1079_1)  \\VariantComplex=(56|56|YIYIYFFNSSDDSNVQP|rs79545619_2) (60|60|LFF*|rs200683836_2)\n" +
            "FERAVTNLSGMSNYHAAEVCEAVVLSPLHRLKPQSQPHPRPINQDLWRWGPGIYIYVYIY\n" +
            "IFLIPQMIPMYSQ\n" +
            ">ENST00000524091.1.2-271 \\Length=90 \\VariantSimple=(20|E|rs327934_1)  \\VariantComplex=(56|56|RH|rr2849_1)\n" +
            "WQQRRRRTGAPRGQIASGAGRARPSNVIYVWRLLGKLWSVCVATCTVGHVFISGWHGQNG\n" +
            "KSVQYVKLGSAERRLSRFMGEGARSPRIPD\n" +
            ">ENSP00000475352.2 \\Length=205  \\VariantComplex=(163|163|GLAGRPRRGGRGARARP*|rs78783575_0)\n" +
            "PAMNMGPGVRGPWASPSGNSIPYSSSSPGSYTGPPGGGGPPGTPIMPSPGDSTNSSENMY\n" +
            "TIMNPIGQGAGRANFPLGPGPEGPMAAMSAMEPHHVNGSLGSGDMDGLPKSSPGAVAGLS\n" +
            "NAPGTPRDDGEMAAAGTFLHPFPSESVSDCVDSPPAAASGRRGWRAGPGGAAGGPEQDRD\n" +
            "RGGPVLARDDHERVMGRQPRASLRA\n";

    @Test
    public void test_findSimpleVariant() {

        String[] fastaLines = fastaStr.split("\n");

        VariantProtDB variantProtDB = new VariantProtDB(fastaLines);

        VariantProtein variantProtein = (VariantProtein) variantProtDB.getProtein("ENSP00000462558.1");

        List<SequenceVariant> variants = variantProtein.matchWithVariant("ESSPDAPDHTSET");
        Assert.assertEquals(0, variants.size());

        variants = variantProtein.matchWithVariant("ESSTOPF");
        Assert.assertEquals(0, variants.size());

        variants = variantProtein.matchWithVariant("AMQSAPGPTPS");
        Assert.assertEquals(1, variants.size());
        Assert.assertEquals("rs3829740_3", variants.get(0).getInfo());

        variants = variantProtein.matchWithVariant("AMQSARGPTPS");
        Assert.assertEquals(0, variants.size());

        variantProtein = (VariantProtein) variantProtDB.getProtein("ENSP00000379873.1");

        variants = variantProtein.matchWithVariant("MAVMAPRTLVLLLSGAL");
        Assert.assertEquals(1, variants.size());
        Assert.assertEquals("rs1143146_0", variants.get(0).getInfo());
    }


    @Test
    public void test2_findSimpleVariant() {

        String[] fastaLines = fastaStr.split("\n");

        VariantProtDB variantProtDB = new VariantProtDB(fastaLines);

        VariantProtein variantProtein = (VariantProtein) variantProtDB.getProtein("ENSP00000379873.1");
        // (86|E|rs1059455_0) (89|G|rs199474430_0) (90|K|rs199474436_0) (91|M|rs79361534_0) (94|H|rs78306866_0)
        // DSDAASQRMEPRAPWIEQEGPEYWDQETRNVKAQSQTDRVDLGTLRGYYNQSEAGSHTIQ
        //                     PEYWDEETGKMKAHSQ
        //                    81    E  GKM  H
        List<SequenceVariant> variants = variantProtein.matchWithVariant("PEYWDEETGKMKAHSQ");
        Assert.assertEquals(5, variants.size());
        Assert.assertEquals("rs1059455_0", variants.get(4).getInfo());
        Assert.assertEquals("rs199474430_0", variants.get(3).getInfo());
        Assert.assertEquals("rs199474436_0", variants.get(2).getInfo());
        Assert.assertEquals("rs79361534_0", variants.get(1).getInfo());
        Assert.assertEquals("rs78306866_0", variants.get(0).getInfo());

        variants = variantProtein.matchWithVariant("PEYWDEETGNMKAQSQ");
        Assert.assertEquals(3, variants.size());
        Assert.assertEquals("rs1059455_0", variants.get(2).getInfo());
        Assert.assertEquals("rs199474430_0", variants.get(1).getInfo());
        Assert.assertEquals("rs79361534_0", variants.get(0).getInfo());

    }


    @Test
    public void test_findComplexVariant() {

        String[] fastaLines = fastaStr.split("\n");

        VariantProtDB variantProtDB = new VariantProtDB(fastaLines);

        //">ENSP00000379873.1 \\Length=365
        // \\VariantSimple=(10|V|rs1143146_0)
        // (33|S|rs2075684_0)
        // (68|K|rs707910_0)
        // (86|E|rs1059455_0)
        // (89|G|rs199474430_0)
        // (90|K|rs199474436_0)
        // (91|M|rs79361534_0)
        // (94|H|rs78306866_0)
        // (100|A|rs1071742_0)
        // (101|N|rs1136688_0)
        // (103|R|rs1136689_0)
        // (114|D|rs1136692_0)
        // (119|L|rs1071743_0)
        // (121|M|rs1136695_0)
        // (123|F|rs1136697_0)
        // (129|P|rs1136700_0)
        // (138|Q|rs3173420_0)
        // (140|Y|rs3173419_0)
        // (151|K|rs1059509_0)
        // (168|Q|rs1059517_0)
        // (174|V|rs1059535_0)
        // (175|R|rs1059536_0)
        // (176|A|rs9256983_0)
        // (180|W|rs9260156_0)
        // (182|V|rs9260157_0)
        // (185|E|rs1059542_0)
        // (187|P|rr1730_0)
        // (187|R|rs3129018_0)
        // (190|D|rs3098019_0)
        // (191|G|rs76185201_0)
        // (251|H|rs145046067_0)
        // (300|P|rs1136903_0)
        // (306|V|rs1136949_0)
        // (307|H|rs41558618_0)
        // (335|N|rs1137160_0)
        // (345|S|rs2231119_0)
        //\VariantComplex=(103|105|GM|rr1725_0) (105|105|LAATTTRARPVLTPSR*|rr1726_0)
        // (107|107|APLLQPERGRFSHHPDNVWLRRGVGRALPPRVPAGRLRRQGLHRPERGPALLDRGGHGGSDHQAQVGGGP*|rr1727_0)
        // (250|250|QTRSSWRPGLQGMEPSRSGRLWWCLLERSRDTPAMCSMRVCPSPSP*|rs45576436_0)
        VariantProtein variantProtein = (VariantProtein) variantProtDB.getProtein("ENSP00000379873.1");

//       61       70        80        90        100       110
//      "DSDAASQRMEPRAPWIEQEGPEYWDQETRNVKAQSQTDRVDLGTLRGYYNQSEAGSHTIQ\n" +
        //                                      VDLGM RGYYNQSETGSHTI
//                                              VDLGM RGYYNQSEAGSHTI
//                                              VDLGM RGYYNQSEDGSHTI

        List<SequenceVariant> variants = variantProtein.matchWithVariant("VDLGMRGYYNQSEAGSHTI");
        Assert.assertEquals(1, variants.size());
        Assert.assertEquals("rr1725_0", variants.get(0).getInfo());

        variants = variantProtein.matchWithVariant("VDLGMRGYYNQSEDGSHTI");
        Assert.assertEquals(2, variants.size());
        Assert.assertEquals("rr1725_0", variants.get(1).getInfo());
        Assert.assertEquals("rs1136692_0", variants.get(0).getInfo());

        variants = variantProtein.matchWithVariant("VDLGMRGYYNQSETGSHTI");
        Assert.assertEquals(0, variants.size());
    }


    @Test
    public void test2_findComplexVariant() {

        String[] fastaLines = fastaStr.split("\n");

        VariantProtDB variantProtDB = new VariantProtDB(fastaLines);

        VariantProtein variantProtein = (VariantProtein) variantProtDB.getProtein("ENST00000495449.1.38-256");

        List<SequenceVariant> variants = variantProtein.matchWithVariant("VYIYFFLI");
        Assert.assertEquals(1, variants.size());
        Assert.assertEquals("rs1079_1", variants.get(0).getInfo());
        Assert.assertEquals("GIYIYVYIYIFLIPQMIP", variantProtein.getWTSequence(variantProtein.getPeptideStart(), "VYIYFFLI".length(), 5));
    }

    //AAASGRRGL: ENSP00000475352.2 STOP:163,163,GLAGRPRRGGRGARARP,rs78783575_0

    @Test
    public void test3_findComplexVariant() {

        String[] fastaLines = fastaStr.split("\n");

        VariantProtDB variantProtDB = new VariantProtDB(fastaLines);

        VariantProtein variantProtein = (VariantProtein) variantProtDB.getProtein("ENST00000524091.1.2-271");

//        ">ENST00000524091.1.2-271 \\Length=90 \\VariantSimple=(20|E|rs327934_1)  \\VariantComplex=(56|56|RH|rr2849_1)\n" +
        //                                                                56    60
//                "WQQRRRRTGAPRGQIASGAGRARPSNVIYVWRLLGKLWSVCVATCTVGHVFISGWHGQNGKSVQYVKLGSAERRLSRFMGEGARSPRIPD\n" +
        //                                                                RGQNGKSVQY
        List<SequenceVariant> variants = variantProtein.matchWithVariant("RHGQNGKSVQY");

        Assert.assertEquals(1, variants.size());

        variants = variantProtein.matchWithVariant("HGQNGKSVQY");
        Assert.assertEquals(0, variants.size());

        variants = variantProtein.matchWithVariant("RGQNGKSVQY");
        Assert.assertEquals(0, variants.size());

        variants = variantProtein.matchWithVariant("FISGWRHGQNGKSVQY");
        Assert.assertEquals(1, variants.size());
    }

    @Test
    public void test4_findComplexVariant() {

        String[] fastaLines = fastaStr.split("\n");

        VariantProtDB variantProtDB = new VariantProtDB(fastaLines);

        VariantProtein variantProtein = (VariantProtein) variantProtDB.getProtein("ENSP00000475352.2");

//        ">ENSP00000475352.2 \\Length=205  \\VariantComplex=(163|163|GLAGRPRRGGRGARARP*|rs78783575_0)\n" +
//                "PAMNMGPGVRGPWASPSGNSIPYSSSSPGSYTGPPGGGGPPGTPIMPSPGDSTNSSENMY\n" +
//                "TIMNPIGQGAGRANFPLGPGPEGPMAAMSAMEPHHVNGSLGSGDMDGLPKSSPGAVAGLS\n" +
//                "NAPGTPRDDGEMAAAGTFLHPFPSESVSDCVDSPPAAASGRRGWRAGPGGAAGGPEQDRD\n" +
//                 ------------------------------------------163
//                 -----------------------------------AAASGRRGL-----------------
//                "RGGPVLARDDHERVMGRQPRASLRA\n";

        List<SequenceVariant> variants = variantProtein.matchWithVariant("AAASGRRGL");
        Assert.assertEquals(1, variants.size());
        Assert.assertEquals("rs78783575_0", variants.get(0).getInfo());
        Assert.assertEquals("AAASGRRGW", variantProtein.getWTSequence(variantProtein.getPeptideStart(), "AAASGRRGL".length(), 0));
        Assert.assertEquals("AAASGRRG", variantProtein.getWTSequence());
    }


}
