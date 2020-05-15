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

import newance.proteinmatch.PeptideUniProtSequenceMatch;
import newance.proteinmatch.VariantInfo;
import newance.proteinmatch.VariantProtDB;
import newance.proteinmatch.VariantProtein;

import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.BiConsumer;

/**
 * @author Markus Müller
 */

public class AddVariantIDs2Psm implements BiConsumer<String,List<PeptideSpectrumMatch>> {

    private final VariantProtDB variantProtDB;

    public AddVariantIDs2Psm(VariantProtDB variantProtDB) {
        this.variantProtDB = variantProtDB;
    }

    public void add(PeptideSpectrumMatch psm) {

        if (!psm.isVariant() || psm.isDecoy()) return;

        Map<String, List<PeptideUniProtSequenceMatch>> matches = variantProtDB.findPeptide(psm.getWTSequence());

        for (String ac : matches.keySet()) {
            for (PeptideUniProtSequenceMatch peptideMatch : matches.get(ac)) {
                VariantInfo variantInfo = ((VariantProtein)peptideMatch.getProtein()).getVariantInfo();
                psm.addVariantAnnot(peptideMatch.getStart(), variantInfo);
            }
        }
    }

    @Override
    public void accept(String s, List<PeptideSpectrumMatch> peptideSpectrumMatchList) {

        for (PeptideSpectrumMatch peptideSpectrumMatch : peptideSpectrumMatchList) {
            add(peptideSpectrumMatch);
        }
    }

}

