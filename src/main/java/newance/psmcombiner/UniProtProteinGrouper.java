package newance.psmcombiner;

import newance.proteinmatch.UniProtDB;
import newance.psmconverter.PeptideMatchData;
import newance.util.PsmGrouper;

import java.util.HashSet;
import java.util.Set;

/**
 * Created by markusmueller on 08.05.18.
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
    public String apply(String s, PeptideMatchData psm) {
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
