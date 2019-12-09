package newance.psmconverter;

import newance.util.NewAnceParams;
import newance.util.PsmPredicate;

import java.io.File;
import java.util.*;
import java.util.concurrent.CountDownLatch;
import java.util.logging.Logger;

/**
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class MaxQuantMSMSEntryConverter  implements Runnable {

    private static final Logger LOGGER = Logger.getLogger(MaxQuantMSMSEntryConverter.class.getName());

    protected final File msmsFile;
    protected final NewAnceParams params;
    protected final Map<String,List<PeptideMatchData>> psms;
    protected final CountDownLatch latch;

    public MaxQuantMSMSEntryConverter(File msmsFile, Map<String,List<PeptideMatchData>> psms, CountDownLatch latch) {

        this.params = NewAnceParams.getInstance();
        this.msmsFile = msmsFile;
        this.psms = psms;
        this.latch = latch;
    }

    public void run() {

        System.out.println("Reading " + msmsFile);

        final Map<String, List<PeptideMatchData>> psmMap = new HashMap<>();

        PsmPredicate psmPredicate = new PsmPredicate(params.getMinCharge(), params.getMaxCharge(), params.getMinPeptideLength(), params.getMaxPeptideLength(), params.getMaxRank(),
                "score", 10.0f, PsmPredicate.ScoreOrder.LARGER);

        PSMReaderCallbackImpl callback = new PSMReaderCallbackImpl(params, new SpectrumKeyFunctionImpl(), psmPredicate, psmMap);

        MaxQuantPsmReader2 psmReader = new MaxQuantPsmReader2();
        psmReader.parse(msmsFile, callback);

        addPsms(psmMap);

        latch.countDown();
        System.out.println("Finished reading " + msmsFile+". Latch count: "+latch.getCount());
    }

    protected void addPsms(Map<String, List<PeptideMatchData>> psmMap) {

        for (String spectrumID : psmMap.keySet()) {

            List<PeptideMatchData> matches = psmMap.get(spectrumID);

            if (psms.containsKey(spectrumID)) {
                psms.get(spectrumID).addAll(matches);
            } else {
                psms.put(spectrumID,Collections.synchronizedList(matches));
            }
        }
    }
}
