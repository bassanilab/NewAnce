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
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

/**
 * @author Markus Müller
 */

public class ProteinModifGrouper extends PsmGrouper {
    private final Map<String,Set<String>> exclude;
    private final List<Pattern> regExps;
    private final List<String> groups;
    private final String modifTag;

    public ProteinModifGrouper(List<Pattern> regExps, List<String> groups) {
        this.regExps = regExps;
        this.groups = groups;
        this.exclude = null;
        this.modifTag = "_mod";
    }

    public ProteinModifGrouper(List<Pattern> regExps, List<String> groups, Map<String,Set<String>> exclude) {
        this.regExps = regExps;
        this.groups = groups;
        this.exclude = exclude;
        this.modifTag = "_mod";
    }

    @Override
    public String apply(String specID, PeptideSpectrumMatch psm) {

        if (!psm.getGroup().isEmpty()) return psm.getGroup();

        boolean isModif = psm.getPeptide().hasModifications();

        List<String> proteins = psm.getProteinIDs();

        if (exclude!=null) {
            for (String protein : proteins) {
                for (String excludedProtGrp : exclude.keySet()) {
                    for (String excludedProt : exclude.get(excludedProtGrp)) {
                        if (protein.contains(excludedProt)) {
                            String grp = excludedProt;
                            if (isModif) grp = grp+modifTag;
                            psm.setGroup(grp);
                            return grp;
                        }
                    }
                }
            }
        }

        for (int i=0;i<regExps.size();i++) {

            for (String protein : proteins) {

                 if (regExps.get(i).matcher(protein).find()) {
                     String grp = groups.get(i);
                     if (isModif) grp = grp+modifTag;
                     psm.setGroup(grp);
                     return grp;
                 }
            }
        }

        String grp = groups.get(groups.size()-1);
        if (isModif) grp = grp+modifTag;
        psm.setGroup(grp);
        return grp;
    }

    @Override
    public String getMasterGroup() {
        return groups.get(0);
    }

    @Override
    public Set<String> getGroups() {
        Set<String> groups_mod = new HashSet<>();
        for (String grp : groups) {
            groups_mod.add(grp);
            groups_mod.add(grp+modifTag);
        }
        return groups_mod;
    }

}
