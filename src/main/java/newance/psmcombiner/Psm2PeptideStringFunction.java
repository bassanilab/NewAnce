package newance.psmcombiner;

import newance.psmconverter.PeptideMatchData;

import java.util.ArrayList;
import java.util.List;
import java.util.function.BiFunction;

/**
 * Copyright (C) 2019
 * @author Markus Müller
 * @Institutions: SIB, Swiss Institute of Bioinformatics; Ludwig Institute for Cancer Research
 */

public class Psm2PeptideStringFunction implements BiFunction<String, List<PeptideMatchData>, List<String>> {

    public enum StringMode {SEQUENCE, MODIF};

    private final StringMode mode;

    public Psm2PeptideStringFunction(StringMode mode) {
        this.mode = mode;
    }

    @Override
    public List<String> apply(String specID, List<PeptideMatchData> peptideMatchData) {

        List<String> peptides = new ArrayList<>();
        for (PeptideMatchData psm : peptideMatchData) {

            String peptide = (mode==StringMode.SEQUENCE)?psm.getPeptide().toSymbolString():psm.getPeptide().toString();

            peptides.add(peptide);
        }

        return peptides;
    }

}
