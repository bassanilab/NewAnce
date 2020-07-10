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

import newance.proteinmatch.*;
import newance.util.NewAnceParams;

import java.util.List;
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

        String seq = psm.getPeptide().toSymbolString();
        String firstProtein = psm.getFirstProteinAC();

        if (psm.isDecoy()) {
            seq = reverse(seq).substring(0,seq.length()-1);
            if (firstProtein.startsWith("DECOY_"))
                firstProtein = firstProtein.substring("DECOY_".length());
        }

        VariantProtein protein = (VariantProtein) variantProtDB.getProtein(firstProtein);

        if (protein == null)
            System.out.println("No protein in VariantDB found for peptide " + seq + " from "+
            firstProtein);
        else {
            if (psm.isVariant()) {
                List<SequenceVariant> variants = protein.matchWithVariant(seq);
                if (variants != null && !variants.isEmpty()) {
                    psm.setSequenceVariants(variants);
                    psm.setPeptideStart(protein.getPeptideStart());
                    psm.setWtSequence(protein.getWTSequence(psm.getPeptideStart(), seq.length(), 10));
                } else {
                    if (!psm.isDecoy()) System.out.println(psm.getPeptide().toSymbolString() + "(" + protein.toString() +
                            ")" + " does not contain variants.");
                    psm.setVariant(false);
                }
            } else {
                int start = protein.getSequence().indexOf(seq);
                psm.setPeptideStart(start);
                psm.setWtSequence(reverse(protein.getWTSequence(start, seq.length(),10)));
            }
        }
    }

    @Override
    public void accept(String s, List<PeptideSpectrumMatch> peptideSpectrumMatchList) {

        for (PeptideSpectrumMatch peptideSpectrumMatch : peptideSpectrumMatchList) {
            add(peptideSpectrumMatch);
        }
    }

    public VariantProtDB getVariantProtDB() {
        return variantProtDB;
    }

    protected static String reverse(String seq) {

        char[] array = seq.toCharArray();
        for (int i = 0; i<array.length/2; i++) {
            char tmp = array[i];
            array[i] = array[array.length-1-i];
            array[array.length-1-i] = tmp;
        }

        return new String(array);
    }


}

