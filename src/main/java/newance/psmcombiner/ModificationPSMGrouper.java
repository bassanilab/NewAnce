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

import newance.util.PsmGrouper;
import org.expasy.mzjava.proteomics.mol.modification.ModAttachment;
import org.expasy.mzjava.proteomics.mol.modification.Modification;
import org.expasy.mzjava.proteomics.mol.modification.ModificationList;
import newance.psmconverter.PeptideSpectrumMatch;

import java.util.HashSet;
import java.util.Set;

/**
 * @author Markus Müller
 */

public class ModificationPSMGrouper extends PsmGrouper {

    public ModificationPSMGrouper() {
    }

    @Override
    public String apply(String s, PeptideSpectrumMatch psm) {

        ModificationList mods = psm.getPeptide().getModifications(ModAttachment.sideChainSet);

        for (Modification mod : mods) {

            if (mod.getLabel().equalsIgnoreCase("phospho"))
                return "Phospho_STY";
            else if (mod.getLabel().equalsIgnoreCase("cysteinyl"))
                return "Cysteinyl_C";
            else if (mod.getLabel().equalsIgnoreCase("oxidation"))
                return "Oxidation_M";
            else if (mod.getLabel().equalsIgnoreCase("carbamidomethyl"))
                return "CAM_C";
            else if (mod.getLabel().equalsIgnoreCase("deamidated"))
                return "Deamid_NQ";
            else if (mod.getLabel().equalsIgnoreCase("acetyl"))
                return "Acetyl_Nt";
        }

        return "None";
    }

    @Override
    public String getMasterGroup() {
        return "None";
    }

    @Override
    public Set<String> getGroups() {
        Set<String> groups = new HashSet<>();
        groups.add("None");
        groups.add("Phospho_STY");
        groups.add("Cysteinyl_C");
        groups.add("Oxidation_M");
        groups.add("CAM_C");
        groups.add("Deamid_NQ");
        groups.add("Acetyl_Nt");

        return groups;
    }

}
