package newance.psmconverter;

/**
 * Copyright (C) 2019, SIB/LICR. All rights reserved
 *
 * SIB, Swiss Institute of Bioinformatics
 * Ludwig Institute for Cancer Research (LICR)
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer. Redistributions in binary
 * form must reproduce the above copyright notice, this list of conditions and
 * the following disclaimer in the documentation and/or other materials provided
 * with the distribution. Neither the name of the SIB/LICR nor the names of
 * its contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL SIB/LICR BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

import com.google.common.base.Optional;
import com.google.common.base.Preconditions;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import org.expasy.mzjava.proteomics.mol.AminoAcid;
import org.expasy.mzjava.proteomics.mol.modification.ModAttachment;
import org.expasy.mzjava.proteomics.mol.modification.Modification;
import org.expasy.mzjava.proteomics.mol.Peptide;
import org.expasy.mzjava.proteomics.ms.ident.ModificationMatch;
import org.expasy.mzjava.proteomics.ms.ident.ModificationMatchResolver;
import org.expasy.mzjava.proteomics.ms.ident.UnresolvableModificationMatchException;

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
    private List<Integer> variantPositions;
    private List<Character> variantWTAAs;
    private Set<String> proteins;

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
        variantPositions = null;
        variantWTAAs = null;
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
        variantPositions = null;
        variantWTAAs = null;
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

    public Set<String> getProteins() {

        return proteins;
    }

    public void setProteins(Set<String> proteins) {
        this.proteins = proteins;
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

    public List<Integer> getVariantPositions() {
        return variantPositions;
    }

    public void setVariantPositions(List<Integer> variantPositions) {
        this.variantPositions = variantPositions;
    }

    public List<Character> getVariantWTAAs() {
        return variantWTAAs;
    }

    public void setVariantWTAAs(List<Character> variantWTAAs) {
        this.variantWTAAs = variantWTAAs;
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

        return sideChainMatchMap.size() + termMatchMap.size();
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


}
