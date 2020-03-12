package newance.psmconverter;

import com.google.common.base.Optional;
import newance.psmconverter.pepxml_v120.AltProteinDataType;
import newance.psmconverter.pepxml_v120.ModInfoDataType;
import newance.psmconverter.pepxml_v120.MsmsPipelineAnalysis;
import newance.psmconverter.pepxml_v120.NameValueType;
import newance.util.NewAnceParams;
import org.expasy.mzjava.core.mol.PeriodicTable;
import org.expasy.mzjava.proteomics.mol.AminoAcid;
import org.expasy.mzjava.proteomics.mol.modification.ModAttachment;
import org.expasy.mzjava.proteomics.mol.modification.Modification;
import org.expasy.mzjava.proteomics.ms.ident.*;


import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBElement;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;
import javax.xml.stream.util.StreamReaderDelegate;
import java.io.*;
import java.util.*;

import static newance.psmconverter.pepxml_v120.MsmsPipelineAnalysis.MsmsRunSummary.SpectrumQuery;
import static newance.psmconverter.pepxml_v120.MsmsPipelineAnalysis.MsmsRunSummary.SpectrumQuery.SearchResult.SearchHit;

/**
 * @author Markus Muller
 */
public class CometPEFFPepXmlReader {

    protected final boolean discardAmbiguous;
    protected final ModificationMatchResolver modMatchResolver;

    protected static final BitSet UNKNOWN_AA = new BitSet();

    public CometPEFFPepXmlReader(boolean discardAmbiguousSequences, ModificationMatchResolver modMatchResolver) {

        this.discardAmbiguous = discardAmbiguousSequences;
        this.modMatchResolver = modMatchResolver;

        UNKNOWN_AA.set('B');
        UNKNOWN_AA.set('J');
        UNKNOWN_AA.set('O');
        UNKNOWN_AA.set('U');
        UNKNOWN_AA.set('Z');
        UNKNOWN_AA.set('X');
    }

    public void parse(File file, PeptideSpectrumMatchList peptideSpectrumMatchList) {

        try {

            parse(new NamespaceRewriteDelegate(XMLInputFactory.newFactory().createXMLStreamReader(new FileInputStream(file), "UTF-8")),
                    peptideSpectrumMatchList);
        } catch (XMLStreamException | JAXBException | FileNotFoundException e) {
            throw new IllegalStateException(e);
        }
    }

    private void parse(XMLStreamReader xsr, PeptideSpectrumMatchList peptideSpectrumMatchList) throws XMLStreamException, JAXBException {

        JAXBContext jc = JAXBContext.newInstance(SpectrumQuery.class);
        Unmarshaller unmarshaller = jc.createUnmarshaller();

        do {

            xsr.next();
            if (xsr.isStartElement() && "spectrum_query".equals(xsr.getLocalName())) {

                JAXBElement<MsmsPipelineAnalysis.MsmsRunSummary.SpectrumQuery> jb =
                        unmarshaller.unmarshal(xsr, MsmsPipelineAnalysis.MsmsRunSummary.SpectrumQuery.class);

                MsmsPipelineAnalysis.MsmsRunSummary.SpectrumQuery query = jb.getValue();

                SpectrumInfo spectrumInfo = new SpectrumInfo(query.getSpectrum());
                spectrumInfo.setScanNumber((int) query.getStartScan());

                spectrumInfo.setPrecursorNeutralMass(query.getPrecursorNeutralMass().doubleValue());
                spectrumInfo.setCharge(query.getAssumedCharge().intValue());
                spectrumInfo.setIndex((int) query.getIndex());
                spectrumInfo.setRetentionTime(query.getRetentionTimeSec().doubleValue()/60.0);
                spectrumInfo.setPrecursorIntensity(0.0);

                for (MsmsPipelineAnalysis.MsmsRunSummary.SpectrumQuery.SearchResult result : query.getSearchResult()) {

                    for (MsmsPipelineAnalysis.MsmsRunSummary.SpectrumQuery.SearchResult.SearchHit searchHit : result.getSearchHit()) {

                        processSearchHit(peptideSpectrumMatchList, spectrumInfo, searchHit);
                    }
                }
            }

            if (xsr.isEndElement() && "msms_pipeline_analysis".equals(xsr.getLocalName()))
                break;
        } while (xsr.hasNext());

        xsr.close();
    }

    protected void processSearchHit(PeptideSpectrumMatchList peptideSpectrumMatchList, SpectrumInfo spectrumInfo, SearchHit searchHit) {

        if (!peptideSpectrumMatchList.isValidSpectrum(spectrumInfo)) return;

        String peptideSequence = searchHit.getPeptide();

        if(discardAmbiguous && containsUnknownAA(peptideSequence))
            return;

        Set<String> proteins = new HashSet<>();
        proteins.add(searchHit.getProtein());

        for (AltProteinDataType altProtein : searchHit.getAlternativeProtein()) {
            proteins.add(altProtein.getProtein());
        }

        if (!peptideSpectrumMatchList.isValidProtein(proteins)) return;

        PeptideMatchDataWrapper peptideMatch = new PeptideMatchDataWrapper(peptideSequence);
        peptideMatch.setRank((int) searchHit.getHitRank());
        peptideMatch.setProteins(proteins);

        boolean isDecoy = containsOnlyProteinPattern(proteins, NewAnceParams.getInstance().getCometDecoyProtPrefix());
        peptideMatch.setDecoy(isDecoy); // must be called before copyModInfo

        int numMissedCleavages = (searchHit.getNumMissedCleavages() != null) ? searchHit.getNumMissedCleavages().intValue() : null;
        peptideMatch.setNumMissedCleavages(numMissedCleavages);

        ModInfoDataType modInfo = searchHit.getModificationInfo();
        copyModInfo(peptideMatch, modInfo);

        for (NameValueType searchScore : searchHit.getSearchScore()) {

            final String name = searchScore.getName();
            final double score = parseDouble(searchScore.getValueAttribute());
            peptideMatch.addScore(name, score);
        }

        peptideMatch.addScore("mass_diff",searchHit.getMassdiff().doubleValue());
        peptideMatch.addScore("tot_num_ions", searchHit.getTotNumIons().doubleValue());
        peptideMatch.addScore("matched_num_ions", searchHit.getNumMatchedIons().doubleValue());

        peptideSpectrumMatchList.resultRead(spectrumInfo, peptideMatch);
    }

    protected boolean containsUnknownAA(String peptideSequence) {

        for(int i = 0, size = peptideSequence.length(); i < size; i++){

            if(UNKNOWN_AA.get(peptideSequence.charAt(i)))
                return true;
        }

        return false;
    }

    protected void copyModInfo(PeptideMatchDataWrapper peptideMatch, ModInfoDataType modInfo) {

        if (modInfo == null) return;

        if (modInfo.getModAminoacidMass() != null) {
            for (ModInfoDataType.ModAminoacidMass modAaMass : modInfo.getModAminoacidMass()) {

                int position = modAaMass.getPosition().intValue() - 1;
                AminoAcid residue = peptideMatch.getAminoAcid(position);
                ModificationMatch modMatch = peptideMatch.addModificationMatch(position, adjustMass(modAaMass.getMass(), residue));
                resolveMod(modMatch);
            }
        }
        if (modInfo.getModNtermMass() != null) {

            ModificationMatch modMatch = peptideMatch.addModificationMatch(ModAttachment.N_TERM, adjustMass(modInfo.getModNtermMass(), ModAttachment.N_TERM));
            resolveMod(modMatch);
        }
        if (modInfo.getModCtermMass() != null) {

            ModificationMatch modMatch = peptideMatch.addModificationMatch(ModAttachment.C_TERM, adjustMass(modInfo.getModCtermMass(), ModAttachment.C_TERM));
            resolveMod(modMatch);
        }

        if (modInfo.getAminoacidSubstitution()!=null && !modInfo.getAminoacidSubstitution().isEmpty() && !peptideMatch.isDecoy()) {

            peptideMatch.setVariant(true);
            List<Integer> positions = new ArrayList<>();
            List<Character> wtAAs = new ArrayList<>();
            for (ModInfoDataType.AminoacidSubstitution aminoacidSubstitution : modInfo.getAminoacidSubstitution()) {

                int position = 0;
                if (aminoacidSubstitution.getPosition()!=null)
                    position = aminoacidSubstitution.getPosition().intValue() - 1;
                char wtAA = aminoacidSubstitution.getOrigAa().charAt(0);
                wtAAs.add(wtAA);
                positions.add(position);
            }

            peptideMatch.setVariantPositions(positions);
            peptideMatch.setVariantWTAAs(wtAAs);
        }
    }

    protected void resolveMod(ModificationMatch modMatch) {

        Optional<Modification> modOpt = modMatchResolver.resolve(modMatch);
        if(modOpt.isPresent())
            modMatch.addPotentialModification(modOpt.get());
    }

    protected double adjustMass(Double mass, AminoAcid residue) {

        return mass - residue.getMassOfMonomer();
    }

    protected double adjustMass(Double mass, ModAttachment modAttachment) {

        switch (modAttachment) {

            case N_TERM:
                return mass - PeriodicTable.H_MASS;
            case C_TERM:
                return mass - PeriodicTable.H_MASS * 2 + PeriodicTable.O_MASS;
            case SIDE_CHAIN:
            default:
                throw new IllegalStateException("Unknown mod attachment " + modAttachment);
        }
    }

    public boolean containsOnlyProteinPattern(Set<String> proteins, String prefix) {

        if (proteins.isEmpty()) return false;

        for (String protein : proteins) {

            if (!protein.startsWith(prefix)) {
                return false;
            }
        }

        return true;
    }


    protected double parseDouble(String number) {

        return number.startsWith("+-") ? Double.parseDouble(number.substring(2)) : Double.parseDouble(number);
    }

    protected static class NamespaceRewriteDelegate extends StreamReaderDelegate {

        public NamespaceRewriteDelegate(XMLStreamReader reader) {

            super(reader);
        }

        @Override
        public String getNamespaceURI() {

            return "http://regis-web.systemsbiology.net/pepXML";
        }
    }
}
