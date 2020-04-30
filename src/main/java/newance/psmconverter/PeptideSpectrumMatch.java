/*
Copyright (C) SIB - Swiss Institute of Bioinformatics, Lausanne, Switzerland
Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
*/

package newance.psmconverter;

import gnu.trove.map.TObjectDoubleMap;
import newance.mzjava.mol.Peptide;

import java.util.Collections;
import java.util.List;
import java.util.Set;

/**
 * @author Markus Müller
 */

public class PeptideSpectrumMatch {

    private final String spectrumFile;
    private final Peptide peptide;
    private final Set<String> proteinAcc;
    private final TObjectDoubleMap<String> scoreMap;
    private final int charge;
    private final int rank;
    private final boolean isDecoy;
    private final float retentionTime;
    private final double neutralPrecMass;
    private final int scanNr;
    private final boolean isVariant;
    private final List<Integer> variantPositions;
    private final List<Character> variantWTAAs;
    private String group;

    public PeptideSpectrumMatch(String spectrumFile, Peptide peptide, Set<String> proteinAccs, TObjectDoubleMap<String> scoreMap, int charge, int rank, float retentionTime,
                                int scanNr, double neutralPrecMass, boolean isDecoy, boolean isVariant, List<Integer> variantPositions, List<Character> variantWTAAs) {

        this.spectrumFile = spectrumFile;
        this.peptide = peptide;
        this.proteinAcc = proteinAccs;
        this.scoreMap = scoreMap;
        this.charge = charge;
        this.isDecoy = isDecoy;
        this.rank = rank;
        this.retentionTime = retentionTime;
        this.scanNr = scanNr;
        this.neutralPrecMass = neutralPrecMass;
        this.isVariant = isVariant;
        this.variantPositions = variantPositions;
        this.variantWTAAs = variantWTAAs;
        this.group = "";
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

    public String getSpectrumFile() {
        return spectrumFile;
    }

    public int getRank() {
        return rank;
    }

    public float getRetentionTime() {
        return retentionTime;
    }

    public double getNeutralPrecMass() {
        return neutralPrecMass;
    }

    public int getScanNr() {
        return scanNr;
    }

    public boolean isVariant() {
        return isVariant;
    }

    public List<Integer> getVariantPositions() {
        return variantPositions;
    }

    public List<Character> getVariantWTAAs() {
        return variantWTAAs;
    }

    public String getGroup() {
        return group;
    }

    public void setGroup(String group) {
        this.group = group;
    }
}
