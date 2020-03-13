/*
Copyright (C) SIB - Swiss Institute of Bioinformatics, Lausanne, Switzerland
Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
*/

package newance.psmcombiner;

import newance.psmconverter.PeptideSpectrumMatch;

import java.util.ArrayList;
import java.util.List;
import java.util.function.BiFunction;

/**
 * @author Markus Müller
 */

public class Psm2PeptideStringFunction implements BiFunction<String, List<PeptideSpectrumMatch>, List<String>> {

    public enum StringMode {SEQUENCE, MODIF};

    private final StringMode mode;

    public Psm2PeptideStringFunction(StringMode mode) {
        this.mode = mode;
    }

    @Override
    public List<String> apply(String specID, List<PeptideSpectrumMatch> peptideSpectrumMatchData) {

        List<String> peptides = new ArrayList<>();
        for (PeptideSpectrumMatch psm : peptideSpectrumMatchData) {

            String peptide = (mode==StringMode.SEQUENCE)?psm.getPeptide().toSymbolString():psm.getPeptide().toString();

            peptides.add(peptide);
        }

        return peptides;
    }

}
