/**
 * Copyright (c) 2010, SIB. All rights reserved.
 *
 * SIB (Swiss Institute of Bioinformatics) - http://www.isb-sib.ch Host -
 * http://mzjava.expasy.org
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer. Redistributions in binary
 * form must reproduce the above copyright notice, this list of conditions and
 * the following disclaimer in the documentation and/or other materials provided
 * with the distribution. Neither the name of the SIB/GENEBIO nor the names of
 * its contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL SIB/GENEBIO BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package newance.mzjava.mol;

import com.google.common.base.Optional;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import com.google.common.collect.Multimap;
import com.google.common.collect.UnmodifiableListIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.procedure.TIntObjectProcedure;
import gnu.trove.procedure.TIntProcedure;
import gnu.trove.procedure.TObjectProcedure;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import newance.mzjava.mol.modification.ModAttachment;
import newance.mzjava.mol.modification.Modification;
import newance.mzjava.mol.modification.ModificationList;
import newance.mzjava.mol.modification.ModificationLists;

import java.util.*;

import static com.google.common.base.Preconditions.*;

/**
 * Abstract amino acid sequence base class. This class holds the amino acid sequence
 * along with the modifications to the amino acids.
 * <p/>
 * Modifications can be attached to the side chains of any amino acid residue or the
 * sequence n-term or c-term.
 * <p/>
 * Each modification site can have 0 or more Modifications
 *
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public abstract class AminoAcidSequence implements SymbolSequence<AminoAcid> {

    private final AminoAcid[] aminoAcids;
    private TIntObjectMap<List<Modification>> sideChainModMap;
    private Map<ModAttachment, List<Modification>> termModMap;
    private final boolean hasAmbiguousAA;

    private int hash = 0;

    /**
     * Construct an AminoAcidSequence by copying the residues
     *
     * @param firstResidue the first residue
     * @param residues     the rest of the residues
     */
    public AminoAcidSequence(AminoAcid firstResidue, AminoAcid... residues) {

        aminoAcids = new AminoAcid[residues.length + 1];

        aminoAcids[0] = firstResidue;
        System.arraycopy(residues, 0, aminoAcids, 1, residues.length);

        hasAmbiguousAA = hasAmbiguousAA();
    }

    /**
     * Constructs an AminoAcidSequence by copying the residues.
     *
     * @param residues the residues, the array cannot be empty
     */
    public AminoAcidSequence(AminoAcid[] residues) {

        Preconditions.checkNotNull(residues);
        int size = residues.length;
        if (size <= 0) throw new IllegalStateException("Peptides cannot be empty");

        aminoAcids = new AminoAcid[size];

        System.arraycopy(residues, 0, aminoAcids, 0, residues.length);

        hasAmbiguousAA = hasAmbiguousAA();
    }

    /**
     * Constructor that copies the residues and modifications from src.
     * The new AminoAcidSequence begins at the specified <code>fromInclusive</code>
     * and extends to the residue at index <code>toExclusive - 1</code>.
     * <p/>
     * Thus the length of the new AminoAcidSequence is <code>toExclusive-fromInclusive</code>.
     *
     * @param src           the AminoAcidSequence from which to copy the sequence and modifications
     * @param fromInclusive the beginning index, inclusive.
     * @param toExclusive   the ending index, exclusive.
     */
    public AminoAcidSequence(final AminoAcidSequence src, final int fromInclusive, final int toExclusive) {

        Preconditions.checkNotNull(src);
        final int size = toExclusive - fromInclusive;
        if (size <= 0) throw new IllegalStateException("Peptides cannot be empty");
        aminoAcids = new AminoAcid[size];

        checkElementIndex(fromInclusive, src.size());
        Preconditions.checkPositionIndex(toExclusive, src.size());

        System.arraycopy(src.aminoAcids, fromInclusive, aminoAcids, 0, aminoAcids.length);

        this.sideChainModMap = copySideChainModifications(fromInclusive, toExclusive, src.sideChainModMap);
        this.termModMap = copyTermModifications(src, fromInclusive, toExclusive, src.termModMap);
        hasAmbiguousAA = hasAmbiguousAA();
    }

    /**
     * Constructor that copies the source AminoAcidSequence and replaces any amino acids in oldAAs
     * with newAA.
     *
     * @param src           the AminoAcidSequence from which to copy the sequence and modifications
     * @param oldAAs        set containing the old amino acids
     * @param newAA         the new amino acid
     */
    public AminoAcidSequence(final AminoAcidSequence src, Set<AminoAcid> oldAAs, AminoAcid newAA) {

        Preconditions.checkNotNull(src);
        Preconditions.checkNotNull(newAA);
        Preconditions.checkNotNull(oldAAs);

        final int size = src.size();
        if (size <= 0) throw new IllegalStateException("Peptides cannot be empty");
        aminoAcids = new AminoAcid[size];

        System.arraycopy(src.aminoAcids, 0, aminoAcids, 0, aminoAcids.length);
        for (int i = 0; i < aminoAcids.length; i++) {

            AminoAcid aa = aminoAcids[i];
            if(oldAAs.contains(aa))
                aminoAcids[i] = newAA;
        }

        this.sideChainModMap = copySideChainModifications(0, size, src.sideChainModMap);
        this.termModMap = copyTermModifications(src, 0, size, src.termModMap);
        hasAmbiguousAA = hasAmbiguousAA();
    }

    private Map<ModAttachment, List<Modification>> copyTermModifications(AminoAcidSequence src, int fromInclusive, int toExclusive, Map<ModAttachment, List<Modification>> srcTermModMap) {

        if (srcTermModMap == null)
            return null;

        Map<ModAttachment, List<Modification>> newTermModMap = new EnumMap<>(ModAttachment.class);
        if (fromInclusive == 0 && src.termModMap.containsKey(ModAttachment.N_TERM)) {

            newTermModMap.put(ModAttachment.N_TERM, src.termModMap.get(ModAttachment.N_TERM));
        }
        if (toExclusive == src.size() && src.termModMap.containsKey(ModAttachment.C_TERM)) {

            newTermModMap.put(ModAttachment.C_TERM, src.termModMap.get(ModAttachment.C_TERM));
        }

        return newTermModMap.isEmpty() ? null : newTermModMap;
    }

    private TIntObjectMap<List<Modification>> copySideChainModifications(final int fromInclusive, final int toExclusive, TIntObjectMap<List<Modification>> srcSideChainModMap) {

        if (srcSideChainModMap == null)
            return null;

        final TIntObjectMap<List<Modification>> newSideChainModMap = new TIntObjectHashMap<>(srcSideChainModMap.keySet().size());
        srcSideChainModMap.forEachEntry(new TIntObjectProcedure<List<Modification>>() {
            @Override
            public boolean execute(int index, List<Modification> modList) {

                if (index >= fromInclusive && index < toExclusive) {

                    newSideChainModMap.put(index - fromInclusive, modList);
                }
                return true;
            }
        });
        return newSideChainModMap.isEmpty() ? null : newSideChainModMap;
    }

    /**
     * Constructs an AminoAcidSequence from the list of residues with the modifications contained in
     * side chain and term mod map.
     *
     * @param residues        the residues, this list must have at least one amino acid
     * @param sideChainModMap the side chain modifications, this map can be empty
     * @param termModMap      the terminal modifications, this map can be empty
     */
    public AminoAcidSequence(List<AminoAcid> residues, Multimap<Integer, Modification> sideChainModMap, Multimap<ModAttachment, Modification> termModMap) {

        final int size = residues.size();
        if (size <= 0) throw new IllegalStateException("Peptides cannot be empty");
        aminoAcids = new AminoAcid[size];

        checkNotNull(sideChainModMap);
        checkNotNull(termModMap);
        checkArgument(!termModMap.containsKey(ModAttachment.SIDE_CHAIN), "termModMap cannot contain side chain modifications");

        residues.toArray(aminoAcids);

        for (int key : sideChainModMap.keySet()) {

            checkElementIndex(key, size());
        }

        if (!sideChainModMap.isEmpty()) {

            this.sideChainModMap = new TIntObjectHashMap<>();
            for (int index : sideChainModMap.keySet()) {

                this.sideChainModMap.put(index, ImmutableList.copyOf(sideChainModMap.get(index)));
            }
        }
        if (!termModMap.isEmpty()) {

            this.termModMap = new EnumMap<>(ModAttachment.class);
            for (ModAttachment modAttachment : termModMap.keySet()) {

                this.termModMap.put(modAttachment, ImmutableList.copyOf(termModMap.get(modAttachment)));
            }
        }

        hasAmbiguousAA = hasAmbiguousAA();
    }


    private boolean hasAmbiguousAA() {

        for (AminoAcid aminoAcid : aminoAcids) {

            if (!aminoAcid.isUnambiguous()) {

                return true;
            }
        }

        return false;
    }


    /**
     * Returns true if this AminoAcidSequence has modifications.
     *
     * @return true if this AminoAcidSequence has modifications
     */
    public boolean hasModifications() {

        return (sideChainModMap != null && !sideChainModMap.isEmpty()) || (termModMap != null && !termModMap.isEmpty());
    }

    public String toString() {

        final StringBuilder sb = new StringBuilder();

        if (termModMap != null && (termModMap.containsKey(ModAttachment.N_TERM))) {

            //noinspection unchecked
            appendMods(sb, termModMap.get(ModAttachment.N_TERM));
            sb.append('_');
        }

        for (int i = 0; i < size(); i++) {

            sb.append(aminoAcids[i]);

            if (sideChainModMap != null && sideChainModMap.containsKey(i)) {

                //noinspection unchecked
                appendMods(sb, sideChainModMap.get(i));
            }
        }

        if (termModMap != null && (termModMap.containsKey(ModAttachment.C_TERM))) {

            sb.append('_');
            //noinspection unchecked
            appendMods(sb, termModMap.get(ModAttachment.C_TERM));
        }

        return sb.toString();
    }

    /**
     * Used by to string
     *
     * @param sb       the StringBuilder
     * @param modLists the list of mod lists
     */
    private void appendMods(StringBuilder sb, List<Modification>... modLists) {

        if (modLists.length == 0) return;

        sb.append('(');
        for (List<Modification> mods : modLists) {

            if (mods == null) mods = Collections.emptyList();
            for (Modification mod : mods) {

                sb.append(mod.getLabel());
                sb.append(", ");
            }
        }
        int length = sb.length();
        sb.delete(length - 2, length);
        sb.append(')');
    }

    /**
     * Returns true if this and <code>sequence</code> have the same AA sequence even if
     * the two sequences have different modifications.
     *
     * @param sequence the sequence to check against
     * @return true if this and <code>sequence</code> have the same AA sequence
     */
    public boolean hasSameSequence(AminoAcidSequence sequence) {

        return Arrays.equals(aminoAcids, sequence.aminoAcids);
    }

    /**
     * Returns true if this and <code>sequence</code> have the same AA sequence.
     *
     * @param sequence the sequence to check against
     * @return true if this and <code>sequence</code> have the same AA sequence
     */
    public boolean hasSameSequence(SymbolSequence<AminoAcid> sequence) {

        if (size() != sequence.size()) return false;

        int seqSize = sequence.size();
        for (int i = 0; i < seqSize; i++) {

            AminoAcid thisAA = getSymbol(i);
            AminoAcid thatAA = sequence.getSymbol(i);
            if (!thisAA.equals(thatAA)) return false;
        }

        return true;
    }

    @Override
    public boolean equals(Object o) {

        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        AminoAcidSequence peptide = (AminoAcidSequence) o;

        return Arrays.equals(aminoAcids, peptide.aminoAcids) &&
                !(sideChainModMap != null ? !sideChainModMap.equals(peptide.sideChainModMap) : peptide.sideChainModMap != null) &&
                !(termModMap != null ? !termModMap.equals(peptide.termModMap) : peptide.termModMap != null);
    }

    @Override
    public int hashCode() {

        int result = hash;

        if (result == 0) {

            result = Arrays.hashCode(aminoAcids);
            result = 31 * result + (sideChainModMap != null ? sideChainModMap.hashCode() : 0);
            result = 31 * result + (termModMap != null ? termModMap.hashCode() : 0);
        }

        this.hash = result;

        return result;
    }

    /**
     * Return the number of modifications that this peptide has.
     *
     * @return the number of modifications that this peptide has
     */
    public int getModificationCount() {

        final Counter count = new Counter();

        if (sideChainModMap != null) {

            sideChainModMap.forEachValue(new TObjectProcedure<List<Modification>>() {
                @Override
                public boolean execute(List<Modification> modifications) {

                    count.increment(modifications.size());
                    return true;
                }
            });
        }


        if (termModMap != null) {
            for (List<Modification> modifications : termModMap.values())
                count.increment(modifications.size());
        }

        return count.getCount();
    }

    /**
     * Returns the number of modifications that are attached to the given attachment
     *
     * @param attachment where the modification is attached
     * @return the number of modifications that are attached to the given attachment
     */
    public int getModificationCount(ModAttachment attachment) {

        switch (attachment) {

            case SIDE_CHAIN:
                return countSideChainModifications();
            case N_TERM:
            case C_TERM:
                return (termModMap != null) ? termModMap.get(attachment).size() : 0;
            default:
                throw new IllegalStateException("cannot get modification count for " + attachment);
        }
    }

    private int countSideChainModifications() {

        if (sideChainModMap == null)
            return 0;

        final Counter counter = new Counter();
        sideChainModMap.forEachValue(new TObjectProcedure<List<Modification>>() {
            @Override
            public boolean execute(List<Modification> modifications) {

                counter.increment(modifications.size());
                return true;
            }
        });
        return counter.getCount();
    }

    /**
     * Returns a set containing all the amino acids that are ambiguous.
     *
     * @return a set containing all the amino acids that are ambiguous
     */
    public Set<AminoAcid> getAmbiguousAminoAcids() {

        if (!hasAmbiguousAA)
            return Collections.emptySet();

        Set<AminoAcid> ambiguousAAs = new HashSet<>();
        for (int i = 0; i < size(); i++) {

            if (!aminoAcids[i].isUnambiguous())
                ambiguousAAs.add(aminoAcids[i]);
        }
        return ambiguousAAs;
    }

    /**
     * Returns true if this AminoAcidSequence has amino acids that are ambiguous such as X or B.
     *
     * @return true if there are ambiguous amino acids false otherwise
     */
    public boolean hasAmbiguousAminoAcids() {

        return hasAmbiguousAA;
    }

    /**
     * Returns a string containing the amino acids and ignores the modifications.
     *
     * @return a string containing the amino acids and ignores the modifications
     */
    @Override
    public String toSymbolString() {

        StringBuilder buff = new StringBuilder();

        for (AminoAcid aa : aminoAcids) {

            buff.append(aa.getSymbol());
        }

        return buff.toString();
    }

    /**
     * Returns the number of residues that have a modification
     *
     * @return the number of residues that have a modification
     */
    public int getModifiedResidueCount() {

        int count = (sideChainModMap != null) ? sideChainModMap.size() : 0;

        if (termModMap != null) {

            if (sideChainModMap == null || !sideChainModMap.containsKey(0))
                count += termModMap.containsKey(ModAttachment.N_TERM) ? 1 : 0;
            if (sideChainModMap == null || !sideChainModMap.containsKey(size() - 1))
                count += termModMap.containsKey(ModAttachment.C_TERM) ? 1 : 0;
        }

        return count;
    }

    /**
     * Returns a ModificationList containing all the modifications for the given attachments.
     * <p><p>
     * For example given the peptide <b>(n-term)_PE(m1)PS(m2, m3)IDE_(c-term1, c-term2)</b>:
     * <pre>
     * getModifications(EnumSet.of(ModAttachment.N_TERM) will return {n-term}
     * getModifications(EnumSet.of(ModAttachment.C_TERM) will return {c-term1, c-term2}
     * getModifications(EnumSet.of(ModAttachment.SIDE_CHAIN) will return {m1, m2, m3}
     * getModifications(EnumSet.of(ModAttachment.all) will return {n-term, m1, m2, m3, c-term1, c-term2}
     * </pre>
     *
     * @param attachments the ModAttachments to include
     * @return a ModificationList containing all the modifications for the given attachments
     */
    public ModificationList getModifications(Set<ModAttachment> attachments) {

        final ModificationListImpl foundMods = new ModificationListImpl();

        if (termModMap != null) {

            for (ModAttachment attachment : attachments) {

                if (termModMap.containsKey(attachment)) foundMods.doAdd(termModMap.get(attachment));
            }
        }

        if (sideChainModMap != null && attachments.contains(ModAttachment.SIDE_CHAIN)) {

            sideChainModMap.forEachValue(new TObjectProcedure<List<Modification>>() {
                @Override
                public boolean execute(List<Modification> mods) {

                    foundMods.doAdd(mods);
                    return true;
                }
            });
        }

        return foundMods;
    }

    /**
     * Returns a ModificationList containing all the modifications for the given index and attachments.
     * <p><p>
     * For example given the peptide <b>(n-term)_PE(m1)PS(m2, m3)IDE_(c-term1, c-term2)</b>:
     * <pre>
     * getModifications(0, EnumSet.of(ModAttachment.N_TERM) will return {n-term}
     * getModifications(6, EnumSet.of(ModAttachment.N_TERM) will return {}
     * getModifications(6, EnumSet.of(ModAttachment.C_TERM) will return {c-term1, c-term2}
     * getModifications(3, EnumSet.of(ModAttachment.SIDE_CHAIN) will return {m2, m3}
     * </pre>
     *
     * @param attachments the ModAttachments to include
     * @return a ModificationList containing all the modifications for the given attachments
     */
    public ModificationList getModificationsAt(int index, Set<ModAttachment> attachments) {

        checkElementIndex(index, size());

        ModificationListImpl foundMods = null;

        if (sideChainModMap != null && attachments.contains(ModAttachment.SIDE_CHAIN) && sideChainModMap.containsKey(index)) {

            foundMods = createAndAddModifications(null, sideChainModMap.get(index));
        }

        if (termModMap != null) {

            if (index <= 0 && attachments.contains(ModAttachment.N_TERM) && termModMap.containsKey(ModAttachment.N_TERM)) {

                foundMods = createAndAddModifications(foundMods, termModMap.get(ModAttachment.N_TERM));
            } else if (index >= size() - 1 && attachments.contains(ModAttachment.C_TERM) && termModMap.containsKey(ModAttachment.C_TERM)) {

                foundMods = createAndAddModifications(foundMods, termModMap.get(ModAttachment.C_TERM));
            }
        }

        return foundMods == null ? ModificationLists.EMPTY_MOD_LIST : foundMods;
    }

    private ModificationListImpl createAndAddModifications(ModificationListImpl foundMods, List<Modification> mods) {

        if(foundMods == null)
            foundMods = new ModificationListImpl();

        foundMods.doAdd(mods);
        return foundMods;
    }

    /**
     * Returns true if the residue at <code>index</code> has a modification.
     *
     * @param index the index
     * @return true if the residue at <code>index</code> has a modification
     */
    public boolean hasModificationAt(int index) {

        final boolean hasSideChainMod = sideChainModMap != null && sideChainModMap.containsKey(index);
        if (index == 0) {

            return hasSideChainMod || (termModMap != null && termModMap.containsKey(ModAttachment.N_TERM));
        } else if (index == aminoAcids.length - 1) {

            return hasSideChainMod || (termModMap != null && termModMap.containsKey(ModAttachment.C_TERM));
        } else {

            return hasSideChainMod;
        }
    }

    /**
     * Returns true if the residue at <code>index</code> has a modification attached by <code>modAttachment</code>.
     *
     * @param index          the index
     * @param modAttachments set of mod attachments to check
     * @return true if the residue at <code>index</code> has a modification
     */
    public boolean hasModificationAt(int index, Set<ModAttachment> modAttachments) {

        boolean hasMods = false;

        for (ModAttachment modAttachment : modAttachments) {

            switch (modAttachment) {

                case SIDE_CHAIN:

                    hasMods = hasMods || sideChainModMap != null && sideChainModMap.containsKey(index);
                    break;
                case N_TERM:

                    hasMods = hasMods || index == 0 && termModMap.containsKey(ModAttachment.N_TERM);
                    break;
                case C_TERM:

                    hasMods = hasMods || index == aminoAcids.length - 1 && termModMap.containsKey(ModAttachment.C_TERM);
                    break;
                default:

                    throw new IllegalStateException(modAttachment + " is an unknown ModAttachment");
            }
        }

        return hasMods;
    }

    /**
     * Returns true if there are modifications attached via the given <code>modAttachment</code>
     *
     * @param modAttachment the mod attachments to check
     * @return true if there are modifications attached via the given <code>modAttachment</code>
     */
    public boolean hasModificationAt(ModAttachment modAttachment) {

        switch (modAttachment) {

            case SIDE_CHAIN:
                return sideChainModMap != null && !sideChainModMap.isEmpty();
            case N_TERM:
            case C_TERM:
                return termModMap != null && termModMap.containsKey(modAttachment);
            default:

                throw new IllegalStateException(modAttachment + " is an unknown ModAttachment");
        }
    }

    /**
     * Returns true if there are modifications attached via any of the ModAttachments in <code>modAttachments</code>
     *
     * @param modAttachments the mod attachments to check
     * @return true if there are modifications attached via the given <code>modAttachment</code>
     */
    public boolean hasModificationAt(Set<ModAttachment> modAttachments) {

        for (ModAttachment modAttachment : modAttachments) {

            if (hasModificationAt(modAttachment)) {

                return true;
            }
        }

        return false;
    }

    /**
     * Returns the number of residues in this amino acid sequence.
     *
     * @return the number of residues in this amino acid sequence
     */
    @Override
    public int size() {

        return aminoAcids.length;
    }

    /**
     * Returns a sorted Array containing the indexes of the modified amino acids
     *
     * @return a sorted Array containing the indexes of the modified amino acids
     */
    public int[] getModificationIndexes(Set<ModAttachment> modAttachments) {

        final TIntSet indexSet = new TIntHashSet();

        if (modAttachments.contains(ModAttachment.N_TERM) && hasModificationAt(ModAttachment.N_TERM)) {
            indexSet.add(0);
        }
        if (modAttachments.contains(ModAttachment.C_TERM) && hasModificationAt(ModAttachment.C_TERM)) {
            indexSet.add(size() - 1);
        }
        if (modAttachments.contains(ModAttachment.SIDE_CHAIN) && sideChainModMap != null) {

            sideChainModMap.forEachKey(new TIntProcedure() {
                @Override
                public boolean execute(int index) {

                    indexSet.add(index);
                    return true;
                }
            });
        }

        int[] indexes = indexSet.toArray();
        Arrays.sort(indexes);
        return indexes;
    }

    /**
     * Returns an array containing the indexes of where <code>symbol</code> is found in this
     * sequence.
     * <p/>
     * For example:
     * <p/>
     * the indexes of S in PSPSK is {1, 3}
     *
     * @param symbol the symbol for which the indexes are to be returned
     * @return an array containing the indexes of where <code>symbol</code> is found in this
     * sequence
     */
    public int[] getSymbolIndexes(AminoAcid symbol) {

        TIntList indices = new TIntArrayList(size());

        for (int i = 0; i < aminoAcids.length; i++) {

            if (aminoAcids[i] == symbol)
                indices.add(i);
        }

        return indices.toArray();
    }

    @Override
    public AminoAcid getSymbol(int index) {

        return aminoAcids[index];
    }

    /**
     * Sums the mass of all the monomers
     *
     * @return the sum of all the monomer mass
     */
    protected double calculateMonomerMassSum() {

        double mass = 0;
        for (int i = 0; i < size(); i++) {

            mass += aminoAcids[i].getMassOfMonomer();
        }

        if (sideChainModMap != null) {

            for (List<Modification> mods : sideChainModMap.valueCollection()) {

                for (Modification mod : mods) {

                    mass += mod.getMolecularMass();
                }
            }
        }

        if (termModMap != null) {

            for (List<Modification> mods : termModMap.values()) {

                for (Modification mod : mods) {
                    mass += mod.getMolecularMass();
                }
            }
        }

        return mass;
    }

    protected double calculateMonomerMassDefect() {

        double massDefect = 0;
        for (int i = 0; i < size(); i++) {

            massDefect += aminoAcids[i].getCompositionOfMonomer().getMassDefect();
        }

        if (sideChainModMap != null) {

            for (List<Modification> mods : sideChainModMap.valueCollection()) {

                for (Modification mod : mods) {

                    massDefect += mod.getMass().getMassDefect();
                }
            }
        }

        if (termModMap != null) {

            for (List<Modification> mods : termModMap.values()) {

                for (Modification mod : mods) {
                    massDefect += mod.getMass().getMassDefect();
                }
            }
        }

        return massDefect;
    }

    /**
     * Counts sums the number of times each of the amino acid in <code>aminoAcidSet</code> is contained in this Peptide.
     * <p/>
     * For example PEPSIDE.countAminoAcidsIn((P, I)) returns 3
     *
     * @param aminoAcidSet the amino acids to check
     * @return the sum of the number of times each of the amino acid in <code>aminoAcidSet</code> is contained in this Peptide
     */
    public int countAminoAcidsIn(Set<AminoAcid> aminoAcidSet) {

        int count = 0;
        for (int i = 0; i < size(); i++) {

            if (aminoAcidSet.contains(getSymbol(i))) count += 1;
        }

        return count;
    }

    /**
     * Calculates and returns the chemical composition of an <code>AminoAcidSequence</code>
     * <p/>
     *
     * @return the chemical composition of an <code>AminoAcidSequence</code>. If there are modifications
     * with unknown composition, Optional.absent() is returned.
     */
    public Optional<Composition> getComposition() {

        Composition.Builder builder = new Composition.Builder();

        for (int i = 0; i < size(); i++) {

            builder.addAll(aminoAcids[i].getCompositionOfMonomer());
        }

        if (sideChainModMap != null) {

            for (List<Modification> mods : sideChainModMap.valueCollection()) {

                for (Modification mod : mods) {

                    try {
                        Composition comp = Composition.parseComposition(mod.getMass().getFormula());
                        builder.addAll(comp);
                    } catch(IllegalArgumentException e) {
                        return Optional.absent();
                    }
                }
            }
        }

        if (termModMap != null) {

            for (List<Modification> mods : termModMap.values()) {

                for (Modification mod : mods) {
                    try {
                        Composition comp = Composition.parseComposition(mod.getMass().getFormula());
                        builder.addAll(comp);
                    } catch(IllegalArgumentException e) {
                        return Optional.absent();
                    }
                }
            }
        }

        return Optional.of(builder.build());
    }

    /**
     * Returns <tt>true</tt> if this AminoAcidSequence contains no residues.
     *
     * @return <tt>true</tt> if this AminoAcidSequence contains no residues
     */
    public boolean isEmpty() {

        return aminoAcids.length == 0;
    }

    private static class ModificationListImpl extends ArrayList<Modification> implements ModificationList, RandomAccess {

        private double weight = 0;

        private ModificationListImpl() {

        }

        private void doAdd(List<Modification> mods) {

            for (Modification mod : mods) {

                weight += mod.getMolecularMass();
            }
            super.addAll(mods);
        }

        @Override
        public double getMolecularMass() {

            return weight;
        }

        @Override
        public Iterator<Modification> iterator() {

            return Iterators.unmodifiableIterator(super.iterator());
        }

        @Override
        public boolean add(Modification modification) {

            throw new UnsupportedOperationException("ModificationList is immutable");
        }

        @Override
        public boolean remove(Object o) {

            throw new UnsupportedOperationException("ModificationList is immutable");
        }

        @Override
        public boolean addAll(Collection<? extends Modification> c) {

            throw new UnsupportedOperationException("ModificationList is immutable");
        }

        @Override
        public boolean addAll(int index, Collection<? extends Modification> c) {

            throw new UnsupportedOperationException("ModificationList is immutable");
        }

        @Override
        public boolean removeAll(Collection<?> c) {

            throw new UnsupportedOperationException("ModificationList is immutable");
        }

        @Override
        public boolean retainAll(Collection<?> c) {

            throw new UnsupportedOperationException("ModificationList is immutable");
        }

        @Override
        public void clear() {

            throw new UnsupportedOperationException("ModificationList is immutable");
        }

        @Override
        public Modification set(int index, Modification element) {

            throw new UnsupportedOperationException("ModificationList is immutable");
        }

        @Override
        public void add(int index, Modification element) {

            throw new UnsupportedOperationException("ModificationList is immutable");
        }

        @Override
        public Modification remove(int index) {

            throw new UnsupportedOperationException("ModificationList is immutable");
        }

        @Override
        public ListIterator<Modification> listIterator() {

            return new WrapperUnmodifiableListIterator(super.listIterator());
        }

        @Override
        public ListIterator<Modification> listIterator(int index) {

            return new WrapperUnmodifiableListIterator(super.listIterator(index));
        }

        @Override
        public List<Modification> subList(int fromIndex, int toIndex) {

            throw new UnsupportedOperationException("ModificationList is immutable");
        }
    }

    private static class WrapperUnmodifiableListIterator extends UnmodifiableListIterator<Modification> {

        private final ListIterator<Modification> it;

        private WrapperUnmodifiableListIterator(ListIterator<Modification> it) {

            this.it = it;
        }

        @Override
        public boolean hasNext() {

            return it.hasNext();
        }

        @Override
        public Modification next() {

            return it.next();
        }

        @Override
        public boolean hasPrevious() {

            return it.hasPrevious();
        }

        @Override
        public Modification previous() {

            return it.previous();
        }

        @Override
        public int nextIndex() {

            return it.nextIndex();
        }

        @Override
        public int previousIndex() {

            return it.previousIndex();
        }
    }
}
