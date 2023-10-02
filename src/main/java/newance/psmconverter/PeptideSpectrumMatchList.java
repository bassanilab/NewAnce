/*
Copyright (C) SIB - Swiss Institute of Bioinformatics, Lausanne, Switzerland
Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
*/

package newance.psmconverter;

import newance.mzjava.mol.Peptide;
import newance.mzjava.mol.modification.UnresolvableModificationMatchException;
import newance.util.PsmPredicate;
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
    protected final SpectrumKeyFunction spectrumKeyFunction;
    protected final SpectrumFilter spectrumFilter;
    protected final PsmPredicate psmPredicate;


    public PeptideSpectrumMatchList(SpectrumKeyFunction spectrumKeyFunction,
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

        List<String> protACs = searchResult.getProteins();

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

        try {
            Peptide peptide = searchResult.toPeptide();

            PeptideSpectrumMatch psm = new PeptideSpectrumMatch(spectrumFile, peptide, protACs,
                    searchResult.getScoreMap(),spectrumInfo.getCharge(),rank, rt,
                    scanNr, precMass, isDecoy, isVariant);


            psmMap.putIfAbsent(key,new ArrayList<>());

            psmMap.get(key).add(psm);

        } catch(UnresolvableModificationMatchException e) {
            System.out.println(e.getMessage());
        }
    }

    public boolean isValidProtein(List<String> proteins) {

        if (proteins==null || proteins.isEmpty()) return false;
        if (excludedProtPattern==null) return true;

        for (String protein : proteins) {
            Matcher matcher = excludedProtPattern.matcher(protein);
            if (!matcher.find() && !protein.isEmpty()) return true;
        }
        return false;
    }

    public boolean isValidSpectrum(SpectrumInfo spectrumInfo) {

        String key = spectrumKeyFunction.apply(spectrumInfo);
        if (spectrumFilter != null && !spectrumFilter.apply(key)) return false;

        return true;
    }


    private List<String> removeProt(List<String> acs, Pattern proteinPattern) {

        List<String> trueACs = new ArrayList<>(acs);

        Matcher matcher;
        for (String ac : acs) {
            matcher = proteinPattern.matcher(ac);
            if (matcher.find()) trueACs.remove(ac);
        }

        return trueACs;
    }

}
