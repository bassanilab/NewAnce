package newance.psmconverter;

import newance.util.NewAnceParams;
import newance.util.PsmPredicate;
import org.expasy.mzjava.core.ms.AbsoluteTolerance;
import org.expasy.mzjava.proteomics.mol.modification.Modification;
import org.expasy.mzjava.proteomics.ms.ident.ModListModMatchResolver;
import org.junit.Test;

import javax.xml.bind.JAXBException;
import javax.xml.stream.XMLStreamException;
import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by markusmueller on 27.01.20.
 */
public class PepXMLParserTest {

    @Test
    public void testPeptXMLParser() throws IOException, XMLStreamException, JAXBException {
        String psmFile = "/Users/markusmueller/Documents/data/0D5P/lncRNA/complete/Comet/20170913_QEh1_LC1_CHC_SA_HLAIp_OD5P_DAC_2_R2.pep.xml";

        final Map<String, List<PeptideSpectrumMatch>> psmMap = new HashMap<>();

        NewAnceParams params = NewAnceParams.getInstance();

        PsmPredicate psmPredicate = new PsmPredicate(params.getMinCharge(), params.getMaxCharge(), params.getMinPeptideLength(), params.getMaxPeptideLength(), params.getMaxRank(),
                params.getCometMainScore(), params.getCometMainScoreMinValue(), PsmPredicate.ScoreOrder.LARGER);

        PeptideSpectrumMatchList peptideSpectrumMatchList = new PeptideSpectrumMatchList(new SpectrumKeyFunctionImpl(), psmPredicate, psmMap);

        Collection<Modification> modifications = params.getModifications();
        ModListModMatchResolver modMatchResolver = new ModListModMatchResolver(new AbsoluteTolerance(params.getModifMatchMassTol()), modifications);

        CometPEFFPepXmlReader psmReader = new CometPEFFPepXmlReader( true, modMatchResolver);
        psmReader.parse(new File(psmFile), peptideSpectrumMatchList);
    }

    @Test
    public void testMaxQuantParser() {
        String msmsFile = "/Users/markusmueller/Documents/data/0D5P/lncRNA/complete/MaxQuant/msms.txt";

        final Map<String, List<PeptideSpectrumMatch>> psmMap = new HashMap<>();

        NewAnceParams params = NewAnceParams.getInstance();

        PsmPredicate psmPredicate = new PsmPredicate(params.getMinCharge(), params.getMaxCharge(), params.getMinPeptideLength(), params.getMaxPeptideLength(), params.getMaxRank(),
                params.getMaxQuantMainScore(), params.getMaxQuantMainScoreMinValue(), PsmPredicate.ScoreOrder.LARGER);

        PeptideSpectrumMatchList peptideSpectrumMatchList = new PeptideSpectrumMatchList(new SpectrumKeyFunctionImpl(), psmPredicate, psmMap);

        MaxQuantPsmReader2 psmReader = new MaxQuantPsmReader2();
        psmReader.parse(new File(msmsFile), peptideSpectrumMatchList);

    }
}
