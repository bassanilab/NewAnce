package newance.psmcombiner;

import newance.psmconverter.PeptideMatchData;
import newance.util.PsmGrouper;

import java.util.HashSet;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Copyright (C) 2019
 * @author Markus Müller
 * @Institutions: SIB, Swiss Institute of Bioinformatics; Ludwig Institute for Cancer Research
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
    public String apply(String specID, PeptideMatchData psm) {

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
