/**
 * Copyright (C) 2019, SIB/LICR. All rights reserved
 *
 * SIB, Swiss Institute of Bioinformatics
 * Ludwig Institute for Cancer Research (LICR)
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer. Redistributions in binary
 * form must reproduce the above copyright notice, this list of conditions and
 * the following disclaimer in the documentation and/or other materials provided
 * with the distribution. Neither the name of the SIB/LICR nor the names of
 * its contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL SIB/LICR BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

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
 * @author Markus Müller
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

        if (!identifier.getScanNumbers().isEmpty()) {
            searchResult.addScore("sn",identifier.getScanNumbers().getFirst().getValue());
        }

        if (identifier.getPrecursorNeutralMass().isPresent()) {
            searchResult.addScore("mass",identifier.getPrecursorNeutralMass().get());
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


        psmMap.putIfAbsent(key,new ArrayList<>());

        psmMap.get(key).add(psm);

//        System.out.println(key+"\t"+psm.getPeptide().toString());
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
