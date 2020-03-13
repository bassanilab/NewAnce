/*
Copyright (C) SIB - Swiss Institute of Bioinformatics, Lausanne, Switzerland
Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
*/


package newance.proteinmatch;

import org.junit.Assert;
import org.junit.Test;

import java.util.List;
import java.util.Map;

/**
 * @author Markus Müller
 */

public class UniProtDBTest {

    private String fastaStr =">sp|P61513|RL37A_HUMAN 60S ribosomal protein L37a OS=Homo sapiens OX=9606 GN=RPL37A PE=1 SV=2\n" +
            "MAKKTKKVGIVGKYGTRYGASLRKMVKKIEISQHAKYTCSFCGKTKMKRRAVGIWHCGSC\n" +
            "MKTVAGGAWTYNTTSAVTVKSAIRRLKELKDQ\n" +
            ">sp|P17812|PYRG1_HUMAN CTP synthase 1 OS=Homo sapiens OX=9606 GN=CTPS1 PE=1 SV=2\n" +
            "MKYILVTGGVISGIGKGIIASSVGTILKSCGLHVTSIKIDPYINIDAGTFSPYEHGEVFV\n" +
            "LDDGGEVDLDLGNYERFLDIRLTKDNNLTTGKIYQYVINKERKGDYLGKTVQVVPHITDA\n" +
            "IQEWVMRQALIPVDEDGLEPQVCVIELGGTVGDIESMPFIEAFRQFQFKVKRENFCNIHV\n" +
            "SLVPQPSSTGEQKTKPTQNSVRELSRPIKPSPPYVCRCSNPLDTSVKEKISMFCHVEPEQ\n" +
            "VICVHDVSSIYRVPLLLEEQGVVDYFLRRLDLPIERQPRKMLMKWKEMADRYDRLLETCS\n" +
            "IALVGKYTKFSDSYASVIKALEHSALAINHKLEIKYIDSADLEPITSQEEPVRYHEAWQK\n" +
            "LCSAHGVLVPGGFGVRGTEGKIQAIAWARNQKKPFLGVCLGMQLAVVEFSRNVLGWQDAN\n" +
            "STEFDPTTSHPVVVDMPEHNPGQMGGTMRLGKRRTLFQTKNSVMRKLYGDADYLEERHRH\n" +
            "RFEVNPVWKKCLEEQGLKFVGQDVEGERMEIVELEDHPFFVGVQYHPEFLSRPIKPSPPY\n" +
            "FGLLLASVGRLSHYLQKGCRLSPRDTYSDRSGSSSPDSEITELKFPSINHD\n" +
            ">sp|P17812-2|PYRG1_HUMAN Isoform 2 of CTP synthase 1 OS=Homo sapiens OX=9606 GN=CTPS1\n" +
            "MFCHVEPEQVICVHDVSSIYRVPLLLEEQGVVDYFLRRLDLPIERQPRKMLMKWKEMADR\n" +
            "YDRLLETCSIALVGKYTKFSDSMAKKIKALEHSALAINHKLEIKYIDSADLEPITSQEEP\n" +
            "VRYHEAWQKLCSAHGVLVPGGFGSRPIKPSPPYIAWARNQKKPFLGVCLGMQLAVVEFSR\n" +
            "NVLGWQDANSTEFDPTTSHPVVVDMPEHNPGQMGGTMRLGKRRTLFQTKNSVMRKLYGDA\n" +
            "DYLEERHRHRFEVNPVWKKCLEEQGLKFVGQDVEGERMEIVELEDHPFFVGVQYHPEFLS\n" +
            "RPIKPSPPYFGLLLASVGRLSHYLQKGCRLSPRDTYSDRSGSSSPDSEITELKFPSINHD\n" +
            ">sp|Q02878|RL6_HUMAN 60S ribosomal protein L6 OS=Homo sapiens OX=9606 GN=RPL6 PE=1 SV=3\n" +
            "MAGEKVEKPDTKEKKPEAKKVDAGGKVKKGNLKAKKPKKGKPHCSRNPVLVRGIGRYSRS\n" +
            "AMYSRKAMYKRKYSAAKSKVEKKKKEKVLATVTKPVGGDKNGGTRVVKLRKMPRYYPTED\n" +
            "VPRKLLSHGKKPFSQHVRKLRASITPGTILIILTGRHRGKRVVFLKQLASGLLLVTGPLV\n" +
            "LNRVPLRRTHQKFVIATSTKIDISNVKIPKHLTDAYFKKKKLRKPRHQEGEIFDTEKEKY\n" +
            "EITEQRKIDQKAVDSQILPKIKAIPQLQGYLRSVFALTNGIYPHKLVF\n" +
            ">sp|P61026|RAB10_HUMAN Ras-related protein Rab-10 OS=Homo sapiens OX=9606 GN=RAB10 PE=1 SV=1\n" +
            "MAKKTYDLLFKLLLIGDSGVGKTCVLFRFSDDAFNTTFISTIGIDFKIKTVELQGKKIKL\n" +
            "QIWDTAGQERFHTITTSYYRGAMSRPIKPSPPYGKSFENISKWLRNIDEHANEDVERMLL\n" +
            "GNKCDMDDKRVVPKGKGEQIAREHGIRFFETSAKANINIEKAFLTLAEDILRKTPVKEPN\n" +
            "SENVDISSGGGVTGWKSKCC\n" +
            ">sp|Q15771|RAB30_HUMAN Ras-related protein Rab-30 OS=Homo sapiens OX=9606 GN=RAB30 PE=1 SV=2\n" +
            "MSMEDYDFLFKIVLIGNAGVGKTCLVRRFTQGLFPPGQGATIGVDFMIKTVEINGEKVKL\n" +
            "QIWDTAGQERFRSITQSYYRSANALILTYDITCEESFRCLPEWLREIEQYASNKVITVLV\n" +
            "GNKIDLAERREVSQQRAEEFSEAQDMYYLETSAKESDNVEKLFLDLACRLISEARQNTLV\n" +
            "NNVSSPLPGEGKSISYLTCCNFN\n";


    @Test
    public void parseHeaderTest() {
        String header = ">tr|A0A024R3B9|A0A024R3B9_HUMAN Alpha-crystallin B chain OS=Homo sapiens GN=CRYAB PE=3 SV=1";
        Map<String,String> headerFieldMap = UniProtDB.parseHeader(header);

        Assert.assertEquals("A0A024R3B9",headerFieldMap.get("uniProtAC"));
        Assert.assertEquals("A0A024R3B9_HUMAN",headerFieldMap.get("standardName"));
        Assert.assertEquals("CRYAB",headerFieldMap.get("geneName"));
        Assert.assertEquals("Alpha-crystallin B chain OS=Homo sapiens GN=CRYAB PE=3 SV=1",headerFieldMap.get("description"));
        Assert.assertEquals("tr",headerFieldMap.get("dbFlag"));

        header = ">sp|A0A183|LCE6A_HUMAN Late cornified envelope protein 6A OS=Homo sapiens GN=LCE6A PE=2 SV=1";
        headerFieldMap = UniProtDB.parseHeader(header);

        Assert.assertEquals("A0A183",headerFieldMap.get("uniProtAC"));
        Assert.assertEquals("LCE6A_HUMAN",headerFieldMap.get("standardName"));
        Assert.assertEquals("LCE6A",headerFieldMap.get("geneName"));
        Assert.assertEquals("Late cornified envelope protein 6A OS=Homo sapiens GN=LCE6A PE=2 SV=1",headerFieldMap.get("description"));
        Assert.assertEquals("sp",headerFieldMap.get("dbFlag"));
    }

    @Test
    public void indexOfTest() {


        String[] fastaLines = fastaStr.split("\n");

        UniProtDB uniProtDB = new UniProtDB(fastaLines);

        String source = "MKYILVTGGVISGIGKGIIASSVGTILKSCGLHVTSIKIDPYINIDAGTFSPYEHGEVFV";

        String target = "MKYILVTGGVI";
        Assert.assertEquals(0,uniProtDB.indexOf(source.toCharArray(),target.toCharArray(),0));

        target = "MKYLIVTGGVL";
        Assert.assertEquals(0,uniProtDB.indexOf(source.toCharArray(),target.toCharArray(),0));

        target = "MKYLVTGGVL";
        Assert.assertEquals(-1,uniProtDB.indexOf(source.toCharArray(),target.toCharArray(),0));

        target = "ISGIGKGIIASSV";
        Assert.assertEquals(10,uniProtDB.indexOf(source.toCharArray(),target.toCharArray(),10));
    }

    @Test
    public void getHashTest() {

        char[] seq = "MKYILVTGGVISGIGKGIIASSVGTILKSCGLHVTSIKIDPYINIDAGTFSPYEHGEVFV".toCharArray();

        for (int i=0;i<seq.length-3;i++) {
            int hash1 = UniProtDB.getHash(seq,i,i+3);
            for (int j=0;j<seq.length-3;j++) {
                int hash2 = UniProtDB.getHash(seq,j,j+3);

                Assert.assertEquals(hash1==hash2,equal(seq,i,i+3,j,j+3));
            }
        }
    }

    private boolean equal(char[] seq, int start1, int end1, int start2, int end2) {

        if (end1-start1 != end2-start2) return false;

        boolean match = true;
        for (int i=start1, j = start2; i<end1 && j <end2; i++,j++) {
            match = match && UniProtDB.match(seq[i],seq[j]);
        }

        return match;
    }

    @Test
    public void testSpeedIndexOf() {


        String[] fastaLines = fastaStr.split("\n");

        UniProtDB uniProtDB = new UniProtDB(fastaLines);

        int ntIter = 1000000;

        char[] target = "GGVISGIGKGII".toCharArray();
        char[] source = "MKYILVTGGVISGIGKGIIASSVGTILKSCGLHVTSIKIDPYINIDAGTFSPYEHGEVFV".toCharArray();

        int len = source.length;

        long start = System.currentTimeMillis();

        for (int i=0;i<ntIter;i++) {

            int j = i%(len-10);
            uniProtDB.indexOf(source,target,j);
        }

        long end = System.currentTimeMillis();

        System.out.println("UniProtDB algo: "+1.0*(end-start)/1000.0);

        start = System.currentTimeMillis();

        for (int i=0;i<ntIter;i++) {

            int j = i%(len-10);
            "MKYILVTGGVISGIGKGIIASSVGTILKSCGLHVTSIKIDPYINIDAGTFSPYEHGEVFV".indexOf("GGVISGIGKGII",j);
        }

        end = System.currentTimeMillis();

        System.out.println("Java algo: "+1.0*(end-start)/1000.0);
    }

    @Test
    public void testFindPeptide() {

        String[] fastaLines = fastaStr.split("\n");

        UniProtDB uniProtDB = new UniProtDB(fastaLines);

        Assert.assertEquals(1, uniProtDB.findPeptide("RQNTLVNNVSSP").size());
        Assert.assertEquals(1, uniProtDB.findPeptide("RQNTIVNNVSSP").size());
        Assert.assertEquals(0, uniProtDB.findPeptide("RQNTXVNNVSSP").size());
        Assert.assertEquals(1, uniProtDB.findPeptide("SLTPGTLILIIT").size());
        Assert.assertEquals(0, uniProtDB.findPeptide("SLTPGLILIIT").size());
        Assert.assertEquals(1, uniProtDB.findPeptide("ISYLTCCNFN").size());
        Assert.assertEquals(3, uniProtDB.findPeptide("SRPIKPSPPY").size());

        Map<String,List<PeptideUniProtSequenceMatch>> matches = uniProtDB.findPeptide("SRPIKPSPPY");

        int cnt = 0;
        for (String ac : matches.keySet()) {
            cnt += matches.get(ac).size();
        }

        Assert.assertEquals(5, cnt);
    }

    @Test
    public void testSpeedFindPeptide() {

        String[] fastaLines = fastaStr.split("\n");

        UniProtDB uniProtDB = new UniProtDB(fastaLines);
        int ntIter = 1000000;

        String source = "MKYILVTGGVISGIGKGIIASSVGTILKSCGLHVTSIKIDPYINIDAGTFSPYEHGEVFV";

        int len = source.length();

        long start = System.currentTimeMillis();

        for (int i=0;i<ntIter;i++) {

            int j = i%(len-10);
            uniProtDB.findPeptide(source.substring(j,j+10));
        }

        long end = System.currentTimeMillis();

        System.out.println("UniProtDB algo: "+1.0*(end-start)/1000.0);
    }

}
