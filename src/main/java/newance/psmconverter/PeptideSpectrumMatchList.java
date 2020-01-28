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
import org.expasy.mzjava.proteomics.mol.Peptide;
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

public class PeptideSpectrumMatchList {

    private static final Logger LOGGER = Logger.getLogger(PeptideSpectrumMatchList.class.getName());

    protected final Map<String, List<PeptideSpectrumMatch>> psmMap;
    protected final NewAnceParams params;
    protected final Pattern decoyProtPattern;
    protected final Pattern excludedProtPattern;
    protected final SpectrumKeyFunction<MsnSpectrum> spectrumKeyFunction;
    protected final SpectrumFilter spectrumFilter;
    protected final PsmPredicate psmPredicate;


    public PeptideSpectrumMatchList(SpectrumKeyFunction<MsnSpectrum> spectrumKeyFunction,
                                    PsmPredicate psmPredicate,
                                    Map<String, List<PeptideSpectrumMatch>> psmMap) {

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

    public void resultRead(SpectrumInfo spectrumInfo, PeptideMatchDataWrapper searchResult)  {

        String key = spectrumKeyFunction.apply(spectrumInfo);
        if (spectrumFilter != null && !spectrumFilter.apply(key)) return;

        boolean isDecoy = searchResult.isDecoy();

        Set<String> protACs = searchResult.getProteins();

        if (!isDecoy) {
            if (decoyProtPattern != null) protACs =  removeProt(protACs,decoyProtPattern);
            if (excludedProtPattern != null) protACs =  removeProt(protACs,excludedProtPattern);
        }

        if (protACs.isEmpty()) return;

        if (!psmPredicate.check(searchResult, spectrumInfo.getCharge())) return;

        String spectrumFile = spectrumInfo.getSpectrumFile();
        float rt = (float) spectrumInfo.getRetentionTime();
        int scanNr = spectrumInfo.getScanNumber();
        double precMass = spectrumInfo.getPrecursorNeutralMass();
        int rank = searchResult.getRank();
        boolean isVariant = searchResult.isVariant();
        List<Integer> variantPositions = searchResult.getVariantPositions();
        List<Character> variantWTAAs = searchResult.getVariantWTAAs();

        Peptide peptide = searchResult.toPeptide();

        PeptideSpectrumMatch psm = new PeptideSpectrumMatch(spectrumFile, peptide, protACs,
                searchResult.getScoreMap(),spectrumInfo.getCharge(),rank, rt,
                scanNr, precMass, isDecoy, isVariant, variantPositions, variantWTAAs);


        psmMap.putIfAbsent(key,new ArrayList<>());

        psmMap.get(key).add(psm);
    }

    public boolean isValidProtein(Set<String> proteins) {

        if (proteins==null || proteins.isEmpty()) return false;
        if (excludedProtPattern==null) return true;

        for (String protein : proteins) {
            Matcher matcher = excludedProtPattern.matcher(protein);
            if (!matcher.find()) {
               return false;
            }
        }
        return true;
    }

    public boolean isValidSpectrum(SpectrumInfo spectrumInfo) {

        String key = spectrumKeyFunction.apply(spectrumInfo);
        if (spectrumFilter != null && !spectrumFilter.apply(key)) return false;

        return true;
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
