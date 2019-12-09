package newance.psmcombiner;

import newance.util.PsmGrouper;
import org.expasy.mzjava.proteomics.mol.modification.ModAttachment;
import org.expasy.mzjava.proteomics.mol.modification.Modification;
import org.expasy.mzjava.proteomics.mol.modification.ModificationList;
import newance.psmconverter.PeptideMatchData;

import java.util.HashSet;
import java.util.Set;

/**
 * Created by markusmueller on 08.05.18.
 */
public class ModificationPSMGrouper extends PsmGrouper {

    public ModificationPSMGrouper() {
    }

    @Override
    public String apply(String s, PeptideMatchData psm) {

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
