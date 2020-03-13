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

import newance.proteinmatch.UniProtDB;
import newance.psmconverter.PeptideSpectrumMatch;
import newance.util.PsmGrouper;

import java.util.HashSet;
import java.util.Set;

/**
 * @author Markus Müller
 */

public class UniProtProteinGrouper extends PsmGrouper {
    private final UniProtDB uniProtDB;
    private final String masterGroup;
    private final String otherGroup;

    public UniProtProteinGrouper(UniProtDB uniProtDB, String otherGroup) {
        this.uniProtDB = uniProtDB;
        this.masterGroup = "UniProt";
        this.otherGroup = otherGroup;
    }

    public UniProtProteinGrouper(UniProtDB uniProtDB, String otherGroup, String masterGroup) {
        this.uniProtDB = uniProtDB;
        this.masterGroup = masterGroup;
        this.otherGroup = otherGroup;
    }

    @Override
    public String apply(String s, PeptideSpectrumMatch psm) {
        String seq = psm.getPeptide().toSymbolString();
        if (psm.isDecoy()) seq = reverse(seq);

        if (uniProtDB.contains(seq.toCharArray()))
            return masterGroup;
        else
            return otherGroup;
    }

    @Override
    public String getMasterGroup() {
        return masterGroup;
    }

    @Override
    public Set<String> getGroups() {
        Set<String> groups = new HashSet<>();
        groups.add(masterGroup);
        groups.add(otherGroup);

        return groups;
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
