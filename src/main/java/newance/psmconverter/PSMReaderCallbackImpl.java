package newance.psmconverter;

import newance.util.PsmPredicate;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.expasy.mzjava.proteomics.io.ms.ident.PSMReaderCallback;
import org.expasy.mzjava.proteomics.mol.Peptide;
import org.expasy.mzjava.proteomics.ms.ident.PeptideProteinMatch;
import org.expasy.mzjava.proteomics.ms.ident.SpectrumIdentifier;
import newance.util.SpectrumFilter;
import newance.util.NewAnceParams;
import newance.util.SpectrumKeyFunction;

import java.util.*;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Copyright (C) 2019
 * @author Markus Müller
 * @Institutions: SIB, Swiss Institute of Bioinformatics; Ludwig Institute for Cancer Research
 */

public class PSMReaderCallbackImpl implements PSMReaderCallback {

    private static final Logger LOGGER = Logger.getLogger(PSMReaderCallbackImpl.class.getName());

    protected final Map<String, List<PeptideMatchData>> psmMap;
    protected final NewAnceParams params;
    protected final Pattern decoyProtPattern;
    protected final Pattern excludedProtPattern;
    protected final SpectrumKeyFunction<MsnSpectrum> spectrumKeyFunction;
    protected final SpectrumFilter spectrumFilter;
    protected final PsmPredicate psmPredicate;


    public PSMReaderCallbackImpl(SpectrumKeyFunction<MsnSpectrum> spectrumKeyFunction,
                                 PsmPredicate psmPredicate,
                                 Map<String, List<PeptideMatchData>> psmMap) {

        this.psmMap = psmMap;
        this.params = NewAnceParams.getInstance();

        if (!params.getCometDecoyProtPrefix().isEmpty()) this.decoyProtPattern = Pattern.compile(params.getCometDecoyProtPrefix());
        else this.decoyProtPattern = null;

        if (params.getExcludedProtPattern()!=null) this.excludedProtPattern = params.getExcludedProtPattern();
        else this.excludedProtPattern = null;

        this.spectrumKeyFunction = spectrumKeyFunction;
        if (params.getSpectrumRegExp()!=null)
            this.spectrumFilter = new SpectrumFilter(params.getSpectrumRegExp());
        else
            this.spectrumFilter = null;

        this.psmPredicate = psmPredicate;

    }

    @Override
    public void resultRead(SpectrumIdentifier identifier, org.expasy.mzjava.proteomics.ms.ident.PeptideMatch searchResult)  {

        if (excludedProtPattern!=null && searchResult.containsOnlyProteinMatch(excludedProtPattern))
            return;

        String key = spectrumKeyFunction.apply(identifier);
        if (spectrumFilter != null && !spectrumFilter.apply(key)) return;

        if (searchResult.getRank().isPresent()) {
            searchResult.addScore("rank", searchResult.getRank().get());
        }

        if (!identifier.getRetentionTimes().isEmpty()) {
            searchResult.addScore("rt",identifier.getRetentionTimes().getFirst().getTime());
        }

        boolean isDecoy = false;
        if (decoyProtPattern != null) isDecoy = searchResult.containsOnlyProteinMatch(decoyProtPattern);

        List<PeptideProteinMatch> ppms = searchResult.getProteinMatches();
        Set<String> protACs = new HashSet<>(ppms.size());
        for (PeptideProteinMatch ppm : ppms) {
            protACs.add(ppm.getAccession());
        }

        if (!isDecoy) {
            if (decoyProtPattern != null) protACs =  removeProt(protACs,decoyProtPattern);
            if (excludedProtPattern != null) protACs =  removeProt(protACs,excludedProtPattern);
        }

        if (protACs.isEmpty()) return;


        if (!psmPredicate.check(searchResult, identifier.getAssumedCharge().get())) return;

        Peptide peptide = searchResult.toPeptide();
        PeptideMatchData psm = new PeptideMatchData(peptide, protACs,
                searchResult.getScoreMap(),identifier.getAssumedCharge().get(),isDecoy);


        if (!psmMap.containsKey(key)) {
            psmMap.put(key,new ArrayList<>());
        }

        psmMap.get(key).add(psm);
    }

    private Set<String> removeProt(Set<String> acs, Pattern proteinPattern) {

        Set<String> trueACs = new HashSet<>(acs);

        Matcher matcher;
        for (String ac : acs) {
            matcher = proteinPattern.matcher(ac);
            if (matcher.find()) trueACs.remove(ac);
        }

        return trueACs;
    }

}
