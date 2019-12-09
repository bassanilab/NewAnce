package newance.psmconverter;

import newance.util.PsmPredicate;
import org.expasy.mzjava.core.ms.AbsoluteTolerance;
import org.expasy.mzjava.proteomics.io.ms.ident.PepXmlReader;
import org.expasy.mzjava.proteomics.mol.modification.Modification;
import org.expasy.mzjava.proteomics.ms.ident.ModListModMatchResolver;
import newance.util.NewAnceParams;

import java.io.File;
import java.util.*;
import java.util.concurrent.CountDownLatch;
import java.util.logging.Logger;
import java.util.stream.Collectors;

/**
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class CometPepXmlEntryConverter implements Runnable {

    private static final Logger LOGGER = Logger.getLogger(CometPepXmlEntryConverter.class.getName());

    protected final File pepXmlFile;
    protected final NewAnceParams params;
    protected final Map<String,List<PeptideMatchData>> psms;
    protected final CountDownLatch latch;

    public CometPepXmlEntryConverter(File pepXmlFile, Map<String,List<PeptideMatchData>> psms, CountDownLatch latch) {

        this.params = NewAnceParams.getInstance();
        this.pepXmlFile = pepXmlFile;
        this.psms = psms;
        this.latch = latch;
    }

    @Override
    public void run() {

        System.out.println("Reading " + pepXmlFile);

        final Map<String, List<PeptideMatchData>> psmMap = new HashMap<>();

        PsmPredicate psmPredicate = new PsmPredicate(params.getMinCharge(), params.getMaxCharge(), params.getMinPeptideLength(), params.getMaxPeptideLength(), params.getMaxRank(),
                "xcorr", 1.0f, PsmPredicate.ScoreOrder.LARGER);

        PSMReaderCallbackImpl callback = new PSMReaderCallbackImpl(params, new SpectrumKeyFunctionImpl(), psmPredicate, psmMap);

        Collection<Modification> modifications = Arrays.asList(params.getModificationStr().split(",")).stream().map(mod -> mod.trim()).map(Modification::parseModification).collect(Collectors.toSet());
        ModListModMatchResolver modMatchResolver = new ModListModMatchResolver(new AbsoluteTolerance(params.getModifMatchMassTol()), modifications);

        CometPepXMLReader psmReader = new CometPepXMLReader(PepXmlReader.ModMassStorage.AA_MASS_PLUS_MOD_MASS, true, modMatchResolver);
        psmReader.parse(pepXmlFile, callback);

        addPsms(psmMap);

        latch.countDown();
        System.out.println("Finished reading " + pepXmlFile+". Latch count: "+latch.getCount());
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
