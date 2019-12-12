package newance.psmconverter;

import newance.proteinmatch.PeptideUniProtSequenceMatch;
import newance.proteinmatch.UniProtDB;

import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.BiConsumer;

/**
 * Copyright (C) 2019
 * @author Markus Müller
 * @Institutions: SIB, Swiss Institute of Bioinformatics; Ludwig Institute for Cancer Research
 */

public class AddUniProtIds implements BiConsumer<String,List<PeptideMatchData>> {

    private final UniProtDB uniProtDB;

    public AddUniProtIds(UniProtDB uniProtDB) {
        this.uniProtDB = uniProtDB;
    }

    public void add(PeptideMatchData psm) {

        String peptideStr = psm.toSymbolString();

        if (psm.isDecoy()) {
            peptideStr = reverse(peptideStr);
        }

        Map<String, List<PeptideUniProtSequenceMatch>> matches = uniProtDB.findPeptide(peptideStr);


        Set<String> protNames = new HashSet<>();

        for (String ac : matches.keySet()) {
            for (PeptideUniProtSequenceMatch peptideMatch : matches.get(ac)) {
                String proteinStr = peptideMatch.getProtein().getFastaUniProtName();
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
    public void accept(String s, List<PeptideMatchData> peptideMatchDataList) {

        for (PeptideMatchData peptideMatchData : peptideMatchDataList) {
            add(peptideMatchData);
        }
    }

}

