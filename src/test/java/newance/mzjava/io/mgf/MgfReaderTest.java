package newance.mzjava.io.mgf;

import newance.mzjava.ms.io.mgf.MgfReader;
import newance.mzjava.ms.io.mgf.MsConvertTitleParser;
import newance.mzjava.ms.spectrum.MsnSpectrum;
import org.junit.Assert;
import org.junit.Test;

import java.io.IOException;
import java.io.StringReader;
import java.util.Map;

/**
 * Created by markusmueller on 20.03.20.
 */
public class MgfReaderTest {
    @Test //review does this test any production code?
    public void testOverridingParseUnknownTag() throws IOException {

        String mgf =
                "BEGIN IONS\n" +
                        "TITLE=20170913_QEh1_LC1_CHC_SA_HLAIp_OD5P_ctrl_2_R1.2.2.1 File:\"20170913_QEh1_LC1_CHC_SA_HLAIp_OD5P_ctrl_2_R1.raw\", NativeID:\"controllerType=0 controllerNumber=1 scan=2\"\n" +
                        "PEPMASS=445.1210834492\t1904656.3750000000\n" +
                        "CHARGE=1+\n" +
                        "SCANS=2\n" +
                        "RTINSECONDS=0.44\n" +
                        "57.1284476067\t1093.4456787109\n" +
                        "59.0388165164\t4442.9296875000\n" +
                        "68.3025339467\t914.6438598633\n" +
                        "69.6124195748\t1156.7563476563\n" +
                        "73.0474787824\t61707.1523437500\n" +
                        "75.0266253817\t2264.8027343750\n" +
                        "79.0215902174\t1226.2163085938\n" +
                        "79.5794894844\t1060.9635009766\n" +
                        "81.0683999603\t1578.1314697266\n" +
                        "91.0576834120\t290889.2187500000\n" +
                        "93.0368139087\t6113.2407226563\n" +
                        "94.5515549743\t1030.1414794922\n" +
                        "103.4858220824\t1173.5506591797\n" +
                        "109.0677369017\t4103.8510742188\n" +
                        "148.4035624141\t1954.3228759766\n" +
                        "149.0446948764\t47402.3437500000\n" +
                        "167.0551652884\t178169.7187500000\n" +
                        "218.9653877967\t2278.2463378906\n" +
                        "219.3252452523\t2547.3107910156\n" +
                        "223.0632741315\t21005.4062500000\n" +
                        "225.0426992719\t50548.5195312500\n" +
                        "227.0219462036\t5190.2172851563\n" +
                        "236.9750940332\t13852.2050781250\n" +
                        "239.0942330346\t4595.9921875000\n" +
                        "254.9860583888\t6052.1269531250\n" +
                        "274.9161430507\t4707.6557617188\n" +
                        "275.0589503399\t2378.9123535156\n" +
                        "285.0094681657\t2171.7495117188\n" +
                        "292.9269359836\t17429.7304687500\n" +
                        "293.0718016238\t1294.2999267578\n" +
                        "299.0624591115\t5583.1533203125\n" +
                        "304.9626642939\t6607.0190429688\n" +
                        "310.9349965911\t6202.7338867188\n" +
                        "312.9529338041\t1516.8671875000\n" +
                        "322.9731782290\t9381.2441406250\n" +
                        "324.9865753285\t17201.9785156250\n" +
                        "330.9613263112\t10956.0966796875\n" +
                        "341.0175947120\t50165.8398437500\n" +
                        "342.9963605803\t47087.7070312500\n" +
                        "348.8603859098\t2009.3826904297\n" +
                        "348.9770688815\t8498.7402343750\n" +
                        "350.9666404642\t1440.0946044922\n" +
                        "355.0683366970\t4794.5034179688\n" +
                        "359.0274931854\t642019.4375000000\n" +
                        "368.8914363310\t2059.7834472656\n" +
                        "369.0413853210\t2023.9552001953\n" +
                        "386.9116315096\t9137.8242187500\n" +
                        "387.0598295892\t1237.5529785156\n" +
                        "398.8390830105\t1316.1827392578\n" +
                        "404.9233429021\t2045.7086181641\n" +
                        "408.8220689769\t14985.5566406250\n" +
                        "416.9592778488\t6797.2202148438\n" +
                        "426.8338829523\t5504.9379882813\n" +
                        "429.0886065844\t9346.5898437500\n" +
                        "444.8180136702\t7240.9594726563\n" +
                        "447.3454701673\t9342.7255859375\n" +
                        "448.3501667733\t1956.0030517578\n" +
                        "END IONS";

        MgfReader reader =
                new MgfReader(new StringReader(mgf), new MsConvertTitleParser());

        MsnSpectrum spectrum = reader.next();

        Assert.assertEquals("20170913_QEh1_LC1_CHC_SA_HLAIp_OD5P_ctrl_2_R1.2.2.1 File:\"20170913_QEh1_LC1_CHC_SA_HLAIp_OD5P_ctrl_2_R1.raw\", NativeID:\"controllerType=0 controllerNumber=1 scan=2\"",spectrum.getComment());
        Assert.assertEquals(2, spectrum.getScanNumbers().getFirst().getValue());
        Assert.assertEquals(1, spectrum.getPrecursor().getCharge());
    }


}
