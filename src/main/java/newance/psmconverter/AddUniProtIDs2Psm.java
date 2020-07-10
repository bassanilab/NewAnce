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
import newance.proteinmatch.UniProtDB;

import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.BiConsumer;

/**
 * @author Markus Müller
 */

public class AddUniProtIDs2Psm implements BiConsumer<String,List<PeptideSpectrumMatch>> {

    private final UniProtDB uniProtDB;

    public AddUniProtIDs2Psm(UniProtDB uniProtDB) {
        this.uniProtDB = uniProtDB;
    }

    public void add(PeptideSpectrumMatch psm) {

        // in case the peptide is mutated, we look for a uniprot entry that contains the mutated sequence
        String peptideStr = psm.toSymbolString();

        if (psm.isDecoy()) {
            peptideStr = reverse(peptideStr);
        }

        Map<String, List<PeptideUniProtSequenceMatch>> matches = uniProtDB.findPeptide(peptideStr);


        Set<String> protNames = new HashSet<>();

        for (String ac : matches.keySet()) {
            for (PeptideUniProtSequenceMatch peptideMatch : matches.get(ac)) {
                String proteinStr = peptideMatch.getProtein().toString();
                if (psm.isDecoy()) proteinStr = "DECOY_"+proteinStr;
                protNames.add(proteinStr);
            }
        }

        psm.addProteinAcc(protNames);
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

    @Override
    public void accept(String s, List<PeptideSpectrumMatch> peptideSpectrumMatchList) {

        for (PeptideSpectrumMatch peptideSpectrumMatch : peptideSpectrumMatchList) {
            add(peptideSpectrumMatch);
        }
    }

    public UniProtDB getUniProtDB() {
        return uniProtDB;
    }
}

