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
import com.google.common.collect.Multimap;
import newance.mzjava.mol.modification.ModAttachment;
import newance.mzjava.mol.modification.Modification;
import newance.mzjava.mol.modification.ModificationResolver;

import java.util.List;
import java.util.Set;

/**
 * @author Oliver Horlacher
 * @author fnikitin
 * @version 1.0
 */
public class
Peptide extends AminoAcidSequence implements Weighable {

    private final double mass;

    /**
     * Construct an Peptide by copying the residues.
     *
     * @param firstResidue the first residue
     * @param residues the rest of the residues
     */
    public Peptide(AminoAcid firstResidue, AminoAcid... residues) {

        super(firstResidue, residues);

        mass = calculateMass();
    }

    /**
     * Constructs an Peptide by copying the residues.
     *
     * @param residues the residues, the array cannot be empty
     */
    public Peptide(List<AminoAcid> residues) {

        super(residues.toArray(new AminoAcid[residues.size()]));

        mass = calculateMass();
    }

    /**
     * Copy constructor.
     *
     * @param src the peptide to copy.
     */
    public Peptide(AminoAcidSequence src) {

        super(src, 0, src.size());

        mass = calculateMass();
    }

    /**
     * Constructs a Peptide by copying the residues and modifications from <code>src</code>.
     * The new Peptide begins at the specified <code>beginIndex</code> and
     * extends to the residue at index <code>endIndex - 1</code>.
     *
     * @param src the AminoAcidSequence from which to copy the sequence and modifications
     * @param beginIndex beginIndex   the beginning index, inclusive.
     * @param endIndex     the ending index, exclusive.
     */
    public Peptide(AminoAcidSequence src, int beginIndex, int endIndex) {

        super(src, beginIndex, endIndex);

        mass = calculateMass();
    }

    /**
     * Constructs an Peptide from the list of residues with the modifications contained in
     * side chain and term mod map.
     *
     * @param residues the residues, this list must have at least one amino acid
     * @param sideChainModMap the side chain modifications, this map can be empty
     * @param termModMap the terminal modifications, this map can be empty
     */
    public Peptide(List<AminoAcid> residues, Multimap<Integer, Modification> sideChainModMap, Multimap<ModAttachment, Modification> termModMap) {

        super(residues, sideChainModMap, termModMap);

        mass = calculateMass();
    }

    public Peptide(Peptide peptide, Set<AminoAcid> oldAAs, AminoAcid newAa) {

        super(peptide, oldAAs, newAa);

        mass = calculateMass();
    }

    private double calculateMass() {

        if (hasAmbiguousAminoAcids()) {

            return -1;
        } else {

            final double peptMass = calculateMonomerMassSum() + PeriodicTable.H_MASS * 2 + PeriodicTable.O_MASS;
            Preconditions.checkState(peptMass > 0);
            return peptMass;
        }
    }

    /**
     * Parse s to build a Peptide.
     * <p/>
     * The expected format is a string of amino acids.  Each amino acid can have an optional list of modifications.
     * The modifications can either be formulas as defined by Composition, or PSI-MS Names found in Uni Mod.
     * <pre>
     * PEPTIDE              : the peptide PEPTIDE
     * PEPTM(O)IDE          : the peptide PEPTIDE with an oxidation on M
     * PEPTM(O, Phospho)IDE : the peptide PEPTIDE with an oxidation and phosphorylation on M
     * (Acetyl)_PEPTIDE_(O) : the peptide PEPTIDE with a n-term Acetylation and c_term oxidation modification
     * </pre>
     *
     * @param s the string
     * @return the Peptide
     * @throws PeptideParseException if s is not a valid string
     */
    public static Peptide parse(String s)  {

        return PeptideParser.getInstance().parsePeptide(s);
    }

    /**
     * Parse <code>s</code> to build a Peptide using the modResolver to convert modification strings to
     * Modification
     *
     * @param s           the string
     * @param modResolver Function to convert modification string to Modification
     * @return the Peptide
     * @throws PeptideParseException if s is not a valid string
     */
    public static Peptide parse(String s, ModificationResolver modResolver) {

        return PeptideParser.getInstance().parsePeptide(s, modResolver);
    }

    /**
     * Returns a new peptide that is a sub sequence of this peptide. The
     * sub sequence begins at the specified <code>beginIndex</code> and
     * extends to the amino acid at index <code>endIndex - 1</code>.
     * Thus the length of the sub sequence is <code>endIndex-beginIndex</code>.
     * <p/>
     *
     * @param beginIndex the beginning index, inclusive.
     * @param endIndex   the ending index, exclusive.
     * @return the specified substring.
     * @throws IndexOutOfBoundsException if the
     *                                   <code>beginIndex</code> is negative, or
     *                                   <code>endIndex</code> is larger than the length of
     *                                   this <code>Peptide</code> object, or
     *                                   <code>beginIndex</code> is larger than
     *                                   <code>endIndex</code>.
     */
    public Peptide subSequence(int beginIndex, int endIndex) {

        Preconditions.checkElementIndex(beginIndex, size());
        Preconditions.checkPositionIndex(endIndex, size());

        return new Peptide(this, beginIndex, endIndex);
    }

    /**
     * Returns a new peptide fragment that is a subsequence of this peptide. The
     * subsequence begins at the specified <code>beginIndex</code> and
     * extends to the amino acid at index <code>endIndex - 1</code>.
     * Thus the length of the subsequence is <code>endIndex-beginIndex</code>.
     * <p/>
     *
     * @param beginIndex   the beginning index, inclusive.
     * @param endIndex     the ending index, exclusive.
     * @param fragmentType the type of the fragment
     * @return a new PeptideFragment with the specified sub sequence.
     * @throws IndexOutOfBoundsException if the
     *                                   <code>beginIndex</code> is negative, or
     *                                   <code>endIndex</code> is larger than the length of
     *                                   this <code>String</code> object, or
     *                                   <code>beginIndex</code> is larger than
     *                                   <code>endIndex</code>.
     */
    public PeptideFragment createFragment(int beginIndex, int endIndex, FragmentType fragmentType) {

        Preconditions.checkNotNull(fragmentType);

        if (fragmentType == FragmentType.FORWARD && beginIndex != 0)

            throw new IllegalArgumentException("begin="+beginIndex+", end="+endIndex+": cannot create "+FragmentType.FORWARD +" fragment!");
        else if (fragmentType == FragmentType.REVERSE && endIndex != size())

            throw new IllegalArgumentException("begin="+beginIndex+", end="+endIndex+": cannot create "+FragmentType.REVERSE +" fragment!");

        int size = size();
        Preconditions.checkElementIndex(beginIndex, size);
        Preconditions.checkElementIndex(endIndex, size + 1);

        return new PeptideFragment(fragmentType, this, beginIndex, endIndex);
    }

    /**
     * Create PeptideFragment that is sub sequence for the give ion type and residue number
     *
     * Given peptide CERVILAS createFragment(y, 1) will return S, createFragment(b, 1) will return C
     * createFragment(b, 3) will return CER.
     *
     * If the IonType is i (immonium) the residue at residueNumber - 1 is returned: createFragment(i, 1) is C,
     * createFragment(i, 2) is E
     *
     * @param ionType the type of the ion, only accepts ions that are FORWARD, REVERSE and MONOMER fragments
     * @param residueNumber the index of the residue. The index of the first residue is 1
     * @return the new PeptideFragment
     */
    public PeptideFragment createFragment(IonType ionType, int residueNumber) {

        Preconditions.checkNotNull(ionType);
        Preconditions.checkArgument(residueNumber > 0);
        if(residueNumber > size()) throw new IllegalStateException(ionType.toString() + residueNumber + " from " + toString());

        FragmentType fragmentType = ionType.getFragmentType();
        switch (fragmentType) {

            case FORWARD:
                return createFragment(0, residueNumber, fragmentType);
            case REVERSE:
                return createFragment(size() - residueNumber, size(), fragmentType);
            case MONOMER:
                return createFragment(residueNumber - 1, residueNumber, fragmentType);
            default:
                throw new IllegalArgumentException("Cannot create a fragment for " + fragmentType + " from ionType " + ionType);
        }
    }

    /**
     * Create a sub sequence for the give ion type and residue number
     *
     * Given peptide CERVILAS subSequence(y, 1) will return S, subSequence(b, 1) will return C
     * subSequence(b, 3) will return CER.
     *
     * If the IonType is i (immonium) the residue at residueNumber - 1 is returned: subSequence(i, 1) is C,
     * subSequence(i, 2) is E
     *
     * @param ionType the type of the ion, only accepts ions that are FORWARD, REVERSE and MONOMER fragments
     * @param residueNumber the index of the residue. The index of the first residue is 1
     * @return the sub sequence
     */
    public Peptide subSequence(IonType ionType, int residueNumber) {

        Preconditions.checkNotNull(ionType);
        Preconditions.checkArgument(residueNumber > 0, "Residue number has to be > 0, was " + residueNumber);
        Preconditions.checkArgument(residueNumber <= size(), ionType.toString() + residueNumber + " from " + toString());

        switch (ionType.getFragmentType()) {

            case FORWARD:
                return subSequence(0, residueNumber);
            case REVERSE:
                return subSequence(size() - residueNumber, size());
            case MONOMER:
                return subSequence(residueNumber - 1, residueNumber);
            default:
                throw new IllegalArgumentException("Cannot create a fragment for " + ionType.getFragmentType() + " from ionType " + ionType);
        }
    }

    /**
     * Calculate the m/z of this peptide at the given charge
     *
     * @param charge the charge
     * @return the m/z of this peptide
     */
    public double calculateMz(int charge) {

        Preconditions.checkArgument(charge != 0, "Cannot calculate m/z, charge was 0");
        if(mass == -1)
            throw new IllegalStateException("m/z cannot be calculated for peptides that have ambiguous amino acids. The peptide is " + toString() + " which contains " + getAmbiguousAminoAcids());

        return AAMassCalculator.getInstance().calculateMz(mass, charge, IonType.p);
    }

    /**
     * Returns the mass of this peptide
     *
     * @return the mass of this peptide
     */
    public double getMolecularMass() {

        if(mass == -1)
            throw new IllegalStateException("molecular mass cannot be calculated for peptides that have ambiguous amino acids. The peptide is " + toString() + " which contains " + getAmbiguousAminoAcids());

        return mass;
    }

    public Peptide replace(Set<AminoAcid> oldAAs, AminoAcid newAa){

        return new Peptide(this, oldAAs, newAa);
    }

    public Optional<Composition> getComposition() {

        Optional<Composition> comp = super.getComposition();

        if (comp.isPresent()) {
            Composition.Builder builder = new Composition.Builder().add(AtomicSymbol.H,2).add(AtomicSymbol.O,1);
            return Optional.of(new Composition(comp.get(),builder.build()));
        } else
            return comp;
    }

}
