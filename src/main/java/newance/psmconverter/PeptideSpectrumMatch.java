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
import newance.proteinmatch.SequenceVariant;
import newance.proteinmatch.VariantInfo;

import java.util.*;

/**
 * @author Markus Müller
 */

public class PeptideSpectrumMatch {

    private final String spectrumFile;
    private final Peptide peptide;
    private final List<String> proteinIDs;
    private final TObjectDoubleMap<String> scoreMap;
    private final int charge;
    private final int rank;
    private final boolean isDecoy;
    private final float retentionTime;
    private final double neutralPrecMass;
    private final int scanNr;
    private final String firstProteinAC;
    private boolean isVariant;
    private final List<SequenceVariant> variants;
    private int peptideStart;
    private String wtSequence;
    private String group;

    public PeptideSpectrumMatch(String spectrumFile, Peptide peptide, List<String> proteinIDs,
                                TObjectDoubleMap<String> scoreMap, int charge, int rank, float retentionTime,
                                int scanNr, double neutralPrecMass, boolean isDecoy, boolean isVariant,
                                List<SequenceVariant> variants) {

        this.spectrumFile = spectrumFile;
        this.peptide = peptide;
        this.proteinIDs = proteinIDs;
        this.scoreMap = scoreMap;
        this.charge = charge;
        this.isDecoy = isDecoy;
        this.rank = rank;
        this.retentionTime = retentionTime;
        this.scanNr = scanNr;
        this.neutralPrecMass = neutralPrecMass;

        if (proteinIDs != null && !proteinIDs.isEmpty())
            this.firstProteinAC = proteinIDs.get(0);
        else
            this.firstProteinAC = "";

        this.isVariant = isVariant;
        if (isVariant) this.variants = variants;
        else this.variants = null;

        this.peptideStart = -1;
        this.wtSequence = "";
        this.group = "";
    }

    public PeptideSpectrumMatch(String spectrumFile, Peptide peptide, List<String> proteinIDs,
                                TObjectDoubleMap<String> scoreMap, int charge, int rank, float retentionTime,
                                int scanNr, double neutralPrecMass, boolean isDecoy, boolean isVariant) {

        this.spectrumFile = spectrumFile;
        this.peptide = peptide;
        this.proteinIDs = proteinIDs;
        this.scoreMap = scoreMap;
        this.charge = charge;
        this.isDecoy = isDecoy;
        this.rank = rank;
        this.retentionTime = retentionTime;
        this.scanNr = scanNr;
        this.neutralPrecMass = neutralPrecMass;

        if (proteinIDs != null && !proteinIDs.isEmpty())
            this.firstProteinAC = proteinIDs.get(0);
        else
            this.firstProteinAC = "";

        this.isVariant = isVariant;
        if (isVariant) this.variants = new ArrayList<>();
        else this.variants = null;

        this.peptideStart = -1;
        this.wtSequence = "";
        this.group = "";
    }

    public Peptide getPeptide() {

        return peptide;
    }

    public List<String> getProteinIDs() {

        return Collections.unmodifiableList(proteinIDs);
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

    public void addScores(TObjectDoubleMap<String> scoreMap) {
        this.scoreMap.putAll(scoreMap);
    }

    public void addProteinAcc(Set<String> newProteinAccs) {
        this.proteinIDs.addAll(newProteinAccs);
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

    public List<SequenceVariant> getVariants() {
        return variants;
    }

    public String getGroup() {
        return group;
    }

    public void setGroup(String group) {
        this.group = group;
    }

    public String getFirstProteinAC() {
        return firstProteinAC;
    }

    public void setSequenceVariants(List<SequenceVariant> variants) {
        this.variants.clear();
        if (variants != null) this.variants.addAll(variants);
        else isVariant = false;
    }

    public void setVariant(boolean variant) {
        isVariant = variant;
    }

    public int getPeptideStart() {
        return peptideStart;
    }

    public void setPeptideStart(int peptideStart) {
        this.peptideStart = peptideStart;
    }

    public String getWtSequence() {
        return wtSequence;
    }

    public void setWtSequence(String wtSequence) {
        this.wtSequence = wtSequence;
    }
}
