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

import com.google.common.base.Optional;
import com.google.common.base.Preconditions;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import newance.mzjava.mol.AminoAcid;
import newance.mzjava.mol.Peptide;
import newance.mzjava.mol.modification.*;
import newance.proteinmatch.SequenceVariant;

import java.util.*;

import static com.google.common.base.Preconditions.*;

/**
 *
 * @author Markus Muller
 */

public class PeptideMatchDataWrapper  {

    private int rank;
    private int numMissedCleavages;
    private final List<AminoAcid> aminoAcids;
    private final String sequence;
    private Map<Integer, List<ModificationMatch>> sideChainMatchMap;
    private Map<ModAttachment, List<ModificationMatch>> termMatchMap;
    private boolean isDecoy;
    private boolean isVariant;
    private List<String> proteins;
    private String leadingProtein;

    private final TObjectDoubleMap<String> scoreMap = new TObjectDoubleHashMap<>();

    private transient final ModificationMatchResolver modificationMatchResolver = new ModificationMatchResolver() {
        @Override
        public Optional<Modification> resolve(ModificationMatch modMatch) {

            if (modMatch == null) return Optional.absent();

            if (modMatch.getCandidateCount() == 0) return Optional.absent();
            else return Optional.of(modMatch.getModificationCandidate(0));
        }
    };

    public PeptideMatchDataWrapper(String peptideSequence) {

        rank = -1;
        numMissedCleavages = -1;
        isDecoy = false;
        isVariant = false;
        proteins = null;
        sideChainMatchMap = null;
        termMatchMap = null;

        sequence = peptideSequence;
        aminoAcids = new ArrayList<>();
        for (char aa : peptideSequence.toCharArray()) {

            try {
                aminoAcids.add(AminoAcid.valueOf(Character.toString(aa)));
            } catch (IllegalArgumentException e) {

                throw new IllegalStateException("peptideSequence = " + peptideSequence, e);
            }
        }
    }

    public PeptideMatchDataWrapper(List<AminoAcid> aminoAcids) {

        rank = -1;
        numMissedCleavages = -1;
        isDecoy = false;
        isVariant = false;
        proteins = null;
        sideChainMatchMap = null;
        termMatchMap = null;

        this.aminoAcids = aminoAcids;

        StringBuilder buff = new StringBuilder();

        for (AminoAcid aa : aminoAcids) {

            buff.append(aa.getSymbol());
        }

        sequence = buff.toString();
    }

    public void addScore(String name, double value) {

        Preconditions.checkNotNull(name);

        scoreMap.put(name, value);
    }

    public TObjectDoubleMap<String> getScoreMap() {

        return scoreMap;
    }

    public double getScore(String scoreName) {

        Preconditions.checkArgument(scoreMap.containsKey(scoreName));

        return scoreMap.get(scoreName);
    }

    public int getRank() {

        return rank;
    }

    public void setRank(int rank) {

        this.rank = rank;
    }

    public int getNumMissedCleavages() {

        return numMissedCleavages;
    }

    public void setNumMissedCleavages(int numMissedCleavages) {

        this.numMissedCleavages = numMissedCleavages;
    }

    public List<String> getProteins() {

        return proteins;
    }

    public void setProteins(List<String> proteins) {
        this.proteins = proteins;
    }

    public String getLeadingProtein() {
        return leadingProtein;
    }

    public void setLeadingProtein(String leadingProtein) {
        this.leadingProtein = leadingProtein;
    }

    public boolean isDecoy() {
        return isDecoy;
    }

    public void setDecoy(boolean decoy) {
        isDecoy = decoy;
    }

    public boolean isVariant() {
        return isVariant;
    }

    public void setVariant(boolean variant) {
        isVariant = variant;
    }

    public String getSequence() {
        return sequence;
    }

    public AminoAcid getAminoAcid(int index) {
        return aminoAcids.get(index);
    }

    /**
     * Add a the <code>modificationMatch</code> to the side chain of residue at <code>index</code>
     *
     * @param index             the index of the reside where the ModificationMatch is to be added
     * @param modificationMatch the ModificationMatch that is to be added
     */
    public void addModificationMatch(int index, ModificationMatch modificationMatch) {

        checkNotNull(modificationMatch);
        checkElementIndex(index, aminoAcids.size());

        if (sideChainMatchMap==null) sideChainMatchMap = new HashMap<>();
        sideChainMatchMap.putIfAbsent(index, new ArrayList<>());
        sideChainMatchMap.get(index).add(modificationMatch);
    }

    /**
     * Convenience method for adding a ModificationMatch.
     *
     * @param index the index of the modification
     * @param modificationMass the mass difference due to the modification
     * @return the ModificationMatch that was added
     */
    public ModificationMatch addModificationMatch(int index, double modificationMass) {

        ModificationMatch modificationMatch = new ModificationMatch(
                modificationMass,
                aminoAcids.get(index),
                index,
                ModAttachment.SIDE_CHAIN
        );
        addModificationMatch(index, modificationMatch);
        return modificationMatch;
    }

    /**
     * Convenience method for adding a ModificationMatch.
     *
     * @param index the index of the modification
     * @param modification the modification
     * @return the ModificationMatch that was added
     */
    public ModificationMatch addModificationMatch(int index, Modification modification) {

        ModificationMatch modificationMatch = new ModificationMatch(
                modification,
                aminoAcids.get(index),
                index,
                ModAttachment.SIDE_CHAIN
        );
        addModificationMatch(index, modificationMatch);
        return modificationMatch;
    }

    /**
     * Add a the <code>modificationMatch</code> to either end of the this peptide match. The end is
     * specified busing the <code>modAttachment</code>
     *
     * @param modAttachment     specifies which end of the peptide the modificationMatch is to be added
     * @param modificationMatch the ModificationMatch that is to be added
     */
    public void addModificationMatch(ModAttachment modAttachment, ModificationMatch modificationMatch) {

        checkNotNull(modAttachment);
        checkNotNull(modificationMatch);
        checkArgument(modAttachment != ModAttachment.SIDE_CHAIN, "To add a side chain modification a residue index is required. Use addModificationMatch(int index, ModificationMatch modificationMatch) instead");

        if (termMatchMap==null) termMatchMap = new HashMap<>();
        termMatchMap.putIfAbsent(modAttachment, new ArrayList<>());
        termMatchMap.get(modAttachment).add(modificationMatch);
    }

    /**
     * Convenience method for adding a ModificationMatch.
     *
     * @param modAttachment the mod attachment is required to be either ModAttachment.N_TERM or ModAttachment.C_TERM
     * @param modificationMass the mass difference due to the modification
     * @return the ModificationMatch that was added
     */
    public ModificationMatch  addModificationMatch(ModAttachment modAttachment, double modificationMass) {

        int index;
        switch (modAttachment) {
            case N_TERM:

                index = 0;
                break;
            case C_TERM:

                index = aminoAcids.size() - 1;
                break;
            case SIDE_CHAIN:
            default:
                throw new IllegalStateException("Cannot add a terminal modification with attachment " + modAttachment);
        }

        ModificationMatch modificationMatch = new ModificationMatch(
                modificationMass,
                aminoAcids.get(index),
                index,
                modAttachment
        );
        addModificationMatch(modAttachment, modificationMatch);
        return modificationMatch;
    }

    /**
     * Convenience method for adding a ModificationMatch.
     *
     * @param modAttachment the mod attachment is required to be either ModAttachment.N_TERM or ModAttachment.C_TERM
     * @param modification the mass modification
     * @return the ModificationMatch that was added
     */
    public ModificationMatch  addModificationMatch(ModAttachment modAttachment, Modification modification) {

        int index;
        switch (modAttachment) {
            case N_TERM:

                index = 0;
                break;
            case C_TERM:

                index = aminoAcids.size() - 1;
                break;
            case SIDE_CHAIN:
            default:
                throw new IllegalStateException("Cannot add a terminal modification with attachment " + modAttachment);
        }

        ModificationMatch modificationMatch = new ModificationMatch(
                modification,
                aminoAcids.get(index),
                index,
                modAttachment
        );
        addModificationMatch(modAttachment, modificationMatch);
        return modificationMatch;
    }

    public int getModificationCount() {

        int cnt = 0;
        if (sideChainMatchMap!=null) cnt = sideChainMatchMap.size();
        if (termMatchMap!=null) cnt += termMatchMap.size();

        return cnt;
    }

    /**
     * Convert this PeptideMatchDataWrapper to a peptide
     *
     * @return the peptide
     * @throws UnresolvableModificationMatchException if a modification match cannot be resolved to a modification
     */
    public Peptide toPeptide() {

        return toPeptide(modificationMatchResolver);
    }

    /**
     * Convert this PeptideMatchDataWrapper to a peptide using the supplied <code>modificationMatchResolver</code>
     * to convert ModificationMatch instances to Modifications.
     *
     * @return the peptide
     * @throws UnresolvableModificationMatchException if a modification match cannot be resolved to a modification
     */
    public Peptide toPeptide(ModificationMatchResolver modMatchResolver) {

        Multimap<Integer, Modification> sideChainModMap = ArrayListMultimap.create();
        Multimap<ModAttachment, Modification> termModMap = ArrayListMultimap.create();
        if (sideChainMatchMap!=null) {
            for (Integer index : sideChainMatchMap.keySet()) {

                for (ModificationMatch modMatch : sideChainMatchMap.get(index)) {

                    Optional<Modification> modOpt = modMatchResolver.resolve(modMatch);
                    if (modOpt.isPresent()) {

                        sideChainModMap.put(index, modOpt.get());
                    } else {
                        throw new UnresolvableModificationMatchException(modMatch);
                    }
                }
            }
        }

        if (termMatchMap!=null) {
            for (ModAttachment attachment : termMatchMap.keySet()) {

                for (ModificationMatch modMatch : termMatchMap.get(attachment)) {

                    Optional<Modification> modOpt = modMatchResolver.resolve(modMatch);
                    if (modOpt.isPresent()) {

                        termModMap.put(attachment, modOpt.get());
                    } else {

                        throw new UnresolvableModificationMatchException(modMatch);
                    }
                }
            }
        }

        return new Peptide(aminoAcids, sideChainModMap, termModMap);
    }

    public Peptide toBarePeptide() {

        return new Peptide(aminoAcids);
    }

    public void addVariants(List<SequenceVariant> variants) {
        variants.addAll(variants);
    }
}
