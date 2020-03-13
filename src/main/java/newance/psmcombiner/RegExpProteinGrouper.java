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
import newance.util.PsmGrouper;

import java.util.HashSet;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author Markus Müller
 */

public class RegExpProteinGrouper extends PsmGrouper {
    private final Pattern regExp;
    private final Set<String> exclude;
    private final String masterGroup;
    private final String otherGroup;

    public RegExpProteinGrouper(Pattern regExp, String masterGroup, String otherGroup) {
        this.regExp = regExp;
        this.masterGroup = masterGroup;
        this.otherGroup = otherGroup;
        this.exclude = null;
    }

    public RegExpProteinGrouper(Pattern regExp, Set<String> exclude, String masterGroup, String otherGroup) {
        this.regExp = regExp;
        this.masterGroup = masterGroup;
        this.otherGroup = otherGroup;
        this.exclude = exclude;
    }

    @Override
    public String apply(String specID, PeptideSpectrumMatch psm) {

        Set<String> proteins = psm.getProteinAcc();

        for (String protein : proteins) {

            if (exclude!=null) {
                for (String excludedProt : exclude) {
                    if (protein.contains(excludedProt))
                        return otherGroup;
                }
            }

            Matcher m = regExp.matcher(protein);
            if (m.find( ))
                return masterGroup;
        }

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

}
