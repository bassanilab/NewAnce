package newance.proteinmatch;

import gnu.trove.map.hash.TObjectDoubleHashMap;
import newance.mzjava.mol.Peptide;
import newance.psmcombiner.CometPsm2StringFunction;
import newance.psmcombiner.Psm2StringFunction;
import newance.psmconverter.AddVariantIDs2Psm;
import newance.psmconverter.PeptideSpectrumMatch;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

/**
 * Created by markusmueller on 14.05.20.
 */
public class AddVariantIDs2PsmTest {

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
            "SLGGPFPLLAPPLLPGLTRALPAAPRAGPQGPGGPWRPWLR\n" +
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
            ">ENSP00000366005.5 \\Length=365 \\VariantSimple=(10|V|rs1143146_0) (33|S|rs2075684_0) (68|K|rs707910_0) (86|E|rs1059455_0) (89|G|rs199474430_0) (90|K|rs199474436_0) (91|M|rs79361534_0) (94|H|rs78306866_0) (100|A|rs1071742_0) (101|N|rs1136688_0) (103|R|rs1136689_0) (114|D|rs1136692_1) (119|L|rs1071743_0) (121|M|rs1136695_0) (123|F|rs1136697_0) (129|P|rs1136700_0) (138|Q|rs3173420_0) (140|Y|rs3173419_0) (151|K|rs1059509_0) (168|Q|rs1059517_0) (174|V|rs1059535_0) (175|R|rs1059536_0) (176|A|rs9256983_0) (180|W|rs9260156_0) (182|V|rs9260157_0) (185|E|rs1059542_0) (187|P|rr1730_0) (187|R|rs3129018_0) (190|D|rs3098019_0) (191|G|rs76185201_0) (251|H|rs145046067_0) (300|P|rs1136903_0) (306|V|rs1136949_0) (307|H|rs41558618_0) (335|N|rs1137160_0) (345|S|rs2231119_0)  \\VariantComplex=(103|105|GM|rr1725_0) (105|105|LAATTTRARPVLTPSR*|rr1726_0) (107|107|APLLQPERGRFSHHPDNVWLRRGVGRALPPRVPAGRLRRQGLHRPERGPALLDRGGHGGSDHQAQVGGGP*|rr1727_0) (250|250|QTRSSWRPGLQGMEPSRSGRLWWCLLERSRDTPAMCSMRVCPSPSP*|rs45576436_0)\n" +
            "MAVMAPRTLLLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRF\n" +
            "DSDAASQRMEPRAPWIEQEGPEYWDQETRNVKAQSQTDRVDLGTLRGYYNQSEAGSHTIQ\n" +
            "IMYGCDVGSDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAAHEAEQL\n" +
            "RAYLDGTCVEWLRRYLENGKETLQRTDPPKTHMTHHPISDHEATLRCWALGFYPAEITLT\n" +
            "WQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWEL\n" +
            "SSQPTIPIVGIIAGLVLLGAVITGAVVAAVMWRRKSSDRKGGSYTQAASSDSAQGSDVSL\n" +
            "TACKV\n" +
            ">ENSP00000365998.2 \\Length=299 \\VariantSimple=(10|V|rs1143146_0) (33|S|rs2075684_0) (68|K|rs707910_0) (86|E|rs1059455_0) (89|G|rs199474430_0) (90|K|rs199474436_0) (91|M|rs79361534_0) (94|H|rs78306866_0) (100|A|rs1071742_0) (101|N|rs1136688_0) (103|R|rs1136689_0) (114|D|rs1136692_2) (119|L|rs1071743_0) (121|M|rs1136695_0) (123|F|rs1136697_0) (129|P|rs1136700_0) (138|Q|rs3173420_0) (140|Y|rs3173419_0) (151|K|rs1059509_0) (168|Q|rs1059517_0) (174|V|rs1059535_0) (175|R|rs1059536_0) (176|A|rs9256983_0) (180|W|rs9260156_0) (182|V|rs9260157_0) (185|E|rs1059542_0) (187|P|rr1730_0) (187|R|rs3129018_0) (190|D|rs3098019_0) (191|G|rs76185201_0) (251|H|rs145046067_0)  \\VariantComplex=(103|105|GM|rr1725_0) (105|105|LAATTTRARPVLTPSR*|rr1726_0) (107|107|APLLQPERGRFSHHPDNVWLRRGVGRALPPRVPAGRLRRQGLHRPERGPALLDRGGHGGSDHQAQVGGGP*|rr1727_0) (250|250|QTRSSWRPGLQGMEPSRSGRLWWCLLERSRDTPAMCSMRVCPSPSP*|rs45576436_0)\n" +
            "MAVMAPRTLLLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRF\n" +
            "DSDAASQRMEPRAPWIEQEGPEYWDQETRNVKAQSQTDRVDLGTLRGYYNQSEAGSHTIQ\n" +
            "IMYGCDVGSDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAAHEAEQL\n" +
            "RAYLDGTCVEWLRRYLENGKETLQRTDPPKTHMTHHPISDHEATLRCWALGFYPAEITLT\n" +
            "WQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWV";


    @Test
    public void test_accept() {

        String[] fastaLines = fastaStr.split("\n");

        VariantProtDB variantProtDB = new VariantProtDB(fastaLines);

        AddVariantIDs2Psm addVariantIDs2Psm = new AddVariantIDs2Psm(variantProtDB);

        Peptide peptide = Peptide.parse("QSEDGSHTIQIMY");

        List<String> proteins = new ArrayList<>();
        proteins.add("ENSP00000379873.1");

        PeptideSpectrumMatch psm = new PeptideSpectrumMatch("Spectrum_file", peptide, proteins,
                new TObjectDoubleHashMap<>(), 2, 1, 100, 100, 1507.661308,
                false, true);


        List<PeptideSpectrumMatch> psms = new ArrayList<>();
        psms.add(psm);

        addVariantIDs2Psm.accept("spec_id", psms);

        Assert.assertEquals(1, psm.getVariants().size());
        Assert.assertEquals("rs1136692_0",psm.getVariants().get(0).getInfo());

    }

    @Test
    public void test_accept2() {

        String[] fastaLines = fastaStr.split("\n");

        VariantProtDB variantProtDB = new VariantProtDB(fastaLines);
        VariantProtein variantProtein = (VariantProtein) variantProtDB.getProtein("ENSP00000379873.1");

        AddVariantIDs2Psm addVariantIDs2Psm = new AddVariantIDs2Psm(variantProtDB);

        Peptide peptide = Peptide.parse("QSEDGSHTLQIMY");

        List<SequenceVariant> variants = new ArrayList<>();

        List<String> proteins = new ArrayList<>();
        proteins.add("ENSP00000379873.1");

        PeptideSpectrumMatch psm = new PeptideSpectrumMatch("Spectrum_file", peptide, proteins,
                new TObjectDoubleHashMap<>(), 2, 1, 100, 100, 1507.661308,
                false, true);


        List<PeptideSpectrumMatch> psms = new ArrayList<>();
        psms.add(psm);

        addVariantIDs2Psm.accept("spec_id", psms);

        Assert.assertEquals(2, psm.getVariants().size());
        Assert.assertEquals("rs1071743_0",psm.getVariants().get(0).getInfo());
        Assert.assertEquals("rs1136692_0",psm.getVariants().get(1).getInfo());

    }


    @Test
    public void test_accept3() {

        String[] fastaLines = fastaStr.split("\n");

        VariantProtDB variantProtDB = new VariantProtDB(fastaLines);

        AddVariantIDs2Psm addVariantIDs2Psm = new AddVariantIDs2Psm(variantProtDB);

        Peptide peptide = Peptide.parse("QSEDGSHTLQMMY");

        List<String> proteins = new ArrayList<>();
        proteins.add("ENSP00000379873.1");
        proteins.add("ENSP00000366002.5");
        proteins.add("ENSP00000366005.5");
        proteins.add("ENSP00000365998.2");

        PeptideSpectrumMatch psm = new PeptideSpectrumMatch("Spectrum_file", peptide, proteins,
                new TObjectDoubleHashMap<>(), 2, 1, 100, 100, 1507.661308,
                false, true);


        List<PeptideSpectrumMatch> psms = new ArrayList<>();
        psms.add(psm);

        addVariantIDs2Psm.accept("spec_id", psms);

        Assert.assertEquals(3, psm.getVariants().size());
        Assert.assertEquals("rs1136695_0",psm.getVariants().get(0).getInfo());
        Assert.assertEquals("rs1071743_0",psm.getVariants().get(1).getInfo());
        Assert.assertEquals("rs1136692_0",psm.getVariants().get(2).getInfo());

    }


    @Test
    public void test_accept4() {

        String[] fastaLines = fastaStr.split("\n");

        VariantProtDB variantProtDB = new VariantProtDB(fastaLines);

        AddVariantIDs2Psm addVariantIDs2Psm = new AddVariantIDs2Psm(variantProtDB);

        Peptide peptide = Peptide.parse("QSEDGSHTLQMMY");

        List<String> proteins = new ArrayList<>();
        proteins.add("ENSP00000379873.1");
        proteins.add("ENSP00000366002.5");
        proteins.add("ENSP00000366005.5");
        proteins.add("ENSP00000365998.2");

        PeptideSpectrumMatch psm = new PeptideSpectrumMatch("Spectrum_file", peptide, proteins,
                new TObjectDoubleHashMap<>(), 2, 1, 100, 100, 1507.661308,
                false, true);


        List<PeptideSpectrumMatch> psms = new ArrayList<>();
        psms.add(psm);

        addVariantIDs2Psm.accept("spec_id", psms);


        CometPsm2StringFunction psm2StringFunction = new CometPsm2StringFunction(null,null);
        System.out.println(psm2StringFunction.getVariantString(psm));

    }


    @Test
    public void test_accept5() {

        String[] fastaLines = fastaStr.split("\n");

        VariantProtDB variantProtDB = new VariantProtDB(fastaLines);

        AddVariantIDs2Psm addVariantIDs2Psm = new AddVariantIDs2Psm(variantProtDB);

        Peptide peptide = Peptide.parse("QSEDGSHTIQLMY");

        List<String> proteins = new ArrayList<>();
        proteins.add("ENSP00000379873.1");

        PeptideSpectrumMatch psm = new PeptideSpectrumMatch("Spectrum_file", peptide, proteins,
                new TObjectDoubleHashMap<>(), 2, 1, 100, 100, 1507.661308,
                false, true);


        List<PeptideSpectrumMatch> psms = new ArrayList<>();
        psms.add(psm);

        addVariantIDs2Psm.accept("spec_id", psms);

        Assert.assertEquals(1, psm.getVariants().size());
        Assert.assertEquals("rs1136692_0",psm.getVariants().get(0).getInfo());

        peptide = Peptide.parse("QSEDGSHTLQLMY");

        psm = new PeptideSpectrumMatch("Spectrum_file", peptide, proteins,
                new TObjectDoubleHashMap<>(), 2, 1, 100, 100, 1507.661308,
                false, true);


        psms.clear();
        psms.add(psm);

        addVariantIDs2Psm.accept("spec_id", psms);

        Assert.assertEquals(2, psm.getVariants().size());
        Assert.assertEquals("rs1071743_0",psm.getVariants().get(0).getInfo());
        Assert.assertEquals("rs1136692_0",psm.getVariants().get(1).getInfo());

    }


    @Test
    public void test_accept7() {

        String fastaEntry = ">sp|O60664|PLIN3_HUMAN Perilipin-3 OS=Homo sapiens OX=9606 GN=PLIN3 PE=1 SV=3\n"+
                "MSADGAEADGSTQVTVEEPVQQPSVVDRVASMPLISSTCDMVSAAYASTKESYPHIKTVC\n"+
                "DAAEKGVRTLTAAAVSGAQPILSKLEPQIASASEYAHRGLDKLEENLPILQQPTEKVLAD\n"+
                "TKELVSSKVSGAQEMVSSAKDTVATQLSEAVDATRGAVQSGVDKTKSVVTGGVQSVMGSR\n"+
                "LGQMVLSGVDTVLGKSEEWADNHLPLTDAELARIATSLDGFDVASVQQQRQEQSYFVRLG\n"+
                "SLSERLRQHAYEHSLGKLRATKQRAQEALLQLSQVLSLMETVKQGVDQKLVEGQEKLHQM\n"+
                "WLSWNQKQLQGPEKEPPKPEQVESRALTMFRDIAQQLQATCTSLGSSIQGLPTNVKDQVQ\n"+
                "QARRQVEDLQATFSSIHSFQDLSSSILAQSRERVASAREALDHMVEYVAQNTPVTWLVGP\n"+
                "FAPGITEKAPEEKK";

        String[] fastaLines = fastaEntry.split("\n");

        VariantProtDB variantProtDB = new VariantProtDB(fastaLines);

        AddVariantIDs2Psm addVariantIDs2Psm = new AddVariantIDs2Psm(variantProtDB);

        Peptide peptide = Peptide.parse("FRDLAQQL");

        List<String> proteins = new ArrayList<>();
        proteins.add("sp|O60664|PLIN3_HUMAN");

        PeptideSpectrumMatch psm = new PeptideSpectrumMatch("Spectrum_file", peptide, proteins,
                new TObjectDoubleHashMap<>(), 2, 1, 100, 100, 1507.661308,
                false, false);


        List<PeptideSpectrumMatch> psms = new ArrayList<>();
        psms.add(psm);

        addVariantIDs2Psm.accept("spec_id", psms);

        Assert.assertEquals("EQVESRALTMFRDIAQQLQATCTSLGSS", psm.getWtSequence());

    }


    @Test
    public void test_accept8() {

        String fastaEntry = ">sp|Q9H0D6|XRN2_HUMAN 5'-3' exoribonuclease 2 OS=Homo sapiens OX=9606 GN=XRN2 PE=1 SV=1\n"+
                "MGVPAFFRWLSRKYPSIIVNCVEEKPKECNGVKIPVDASKPNPNDVEFDNLYLDMNGIIH\n"+
                "PCTHPEDKPAPKNEDEMMVAIFEYIDRLFSIVRPRRLLYMAIDGVAPRAKMNQQRSRRFR\n"+
                "ASKEGMEAAVEKQRVREEILAKGGFLPPEEIKERFDSNCITPGTEFMDNLAKCLRYYIAD\n"+
                "RLNNDPGWKNLTVILSDASAPGEGEHKIMDYIRRQRAQPNHDPNTHHCLCGADADLIMLG\n"+
                "LATHEPNFTIIREEFKPNKPKPCGLCNQFGHEVKDCEGLPREKKGKHDELADSLPCAEGE\n"+
                "FIFLRLNVLREYLERELTMASLPFTFDVERSIDDWVFMCFFVGNDFLPHLPSLEIRENAI\n"+
                "DRLVNIYKNVVHKTGGYLTESGYVNLQRVQMIMLAVGEVEDSIFKKRKDDEDSFRRRQKE\n"+
                "KRKRMKRDQPAFTPSGILTPHALGSRNSPGSQVASNPRQAAYEMRMQNNSSPSISPNTSF\n"+
                "TSDGSPSPLGGIKRKAEDSDSEPEPEDNVRLWEAGWKQRYYKNKFDVDAADEKFRRKVVQ\n"+
                "SYVEGLCWVLRYYYQGCASWKWYYPFHYAPFASDFEGIADMPSDFEKGTKPFKPLEQLMG\n"+
                "VFPAASGNFLPPSWRKLMSDPDSSIIDFYPEDFAIDLNGKKYAWQGVALLPFVDERRLRA\n"+
                "ALEEVYPDLTPEETRRNSLGGDVLFVGKHHPLHDFILELYQTGSTEPVEVPPELCHGIQG\n"+
                "KFSLDEEAILPDQIVCSPVPMLRDLTQNTVVSINFKDPQFAEDYIFKAVMLPGARKPAAV\n"+
                "LKPSDWEKSSNGRQWKPQLGFNRDRRPVHLDQAAFRTLGHVMPRGSGTGIYSNAAPPPVT\n"+
                "YQGNLYRPLLRGQAQIPKLMSNMRPQDSWRGPPPLFQQQRFDRGVGAEPLLPWNRMLQTQ\n"+
                "NAAFQPNQYQMLAGPGGYPPRRDDRGGRQGYPREGRKYPLPPPSGRYNWN";

        String[] fastaLines = fastaEntry.split("\n");

        VariantProtDB variantProtDB = new VariantProtDB(fastaLines);

        AddVariantIDs2Psm addVariantIDs2Psm = new AddVariantIDs2Psm(variantProtDB);

        Peptide peptide = Peptide.parse("IDRIVNIY");

        List<String> proteins = new ArrayList<>();
        proteins.add("sp|Q9H0D6|XRN2_HUMAN");

        PeptideSpectrumMatch psm = new PeptideSpectrumMatch("Spectrum_file", peptide, proteins,
                new TObjectDoubleHashMap<>(), 2, 1, 100, 100, 1507.661308,
                false, false);


        List<PeptideSpectrumMatch> psms = new ArrayList<>();
        psms.add(psm);

        addVariantIDs2Psm.accept("spec_id", psms);

        Assert.assertEquals("LPSLEIRENAIDRLVNIYKNVVHKTGGY", psm.getWtSequence());

    }



}
