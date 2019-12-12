package newance.psmconverter;

import gnu.trove.impl.unmodifiable.TUnmodifiableObjectDoubleMap;
import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import org.expasy.mzjava.proteomics.mol.Peptide;

import java.util.Collections;
import java.util.Set;

/**
 * Copyright (C) 2019
 * @author Markus Müller
 * @Institutions: SIB, Swiss Institute of Bioinformatics; Ludwig Institute for Cancer Research
 */

public class PeptideMatchData {

    private final Peptide peptide;
    private final Set<String> proteinAcc;
    private final TObjectDoubleMap<String> scoreMap;
    private final int charge;
    private final boolean isDecoy;

    public PeptideMatchData(Peptide peptide, Set<String> proteinAccs, TObjectDoubleMap<String> scoreMap, int charge) {

        this.peptide = peptide;
        this.proteinAcc = proteinAccs;
        this.scoreMap = new TUnmodifiableObjectDoubleMap<>(new TObjectDoubleHashMap<>(scoreMap));
        this.charge = charge;
        this.isDecoy = false;
    }

    public PeptideMatchData(Peptide peptide, Set<String> proteinAccs,
                            TObjectDoubleMap<String> scoreMap, int charge, boolean isDecoy) {

        this.peptide = peptide;
        this.proteinAcc = proteinAccs;
        this.scoreMap = new TUnmodifiableObjectDoubleMap<>(new TObjectDoubleHashMap<>(scoreMap));
        this.charge = charge;
        this.isDecoy = isDecoy;
    }

    public Peptide getPeptide() {

        return peptide;
    }

    public Set<String> getProteinAcc() {

        return Collections.unmodifiableSet(proteinAcc);
    }

    public TObjectDoubleMap<String> getScoreMap() {

        return scoreMap;
    }

    public double getScore(String score) {

        return scoreMap.get(score);
    }

    public int getCharge() {

        return charge;
    }

    public String toSymbolString() {

        return peptide.toSymbolString();
    }

    public boolean isDecoy() {

        return isDecoy;
    }

    public void addProteinAcc(Set<String> newProteinAccs) {
        this.proteinAcc.addAll(newProteinAccs);
    }
}
