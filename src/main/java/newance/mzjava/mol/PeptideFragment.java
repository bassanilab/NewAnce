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

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import newance.mzjava.mol.modification.ModAttachment;
import newance.mzjava.mol.modification.Modification;
import newance.mzjava.mol.modification.ModificationResolver;

import java.util.ArrayList;
import java.util.List;

import static com.google.common.base.Preconditions.checkArgument;
import static com.google.common.base.Preconditions.checkNotNull;

/**
 * @author Oliver Horlacher
 * @version 1.0
 */
public class PeptideFragment extends AminoAcidSequence {

    private final FragmentType fragmentType;

    /**
     * Construct an PeptideFragment by copying the residues
     *
     * @param fragmentType the fragment type
     * @param firstResidue the first residue
     * @param residues the rest of the residues
     */
    public PeptideFragment(FragmentType fragmentType, AminoAcid firstResidue, AminoAcid... residues) {

        super(firstResidue, residues);

        checkNotNull(fragmentType);
        this.fragmentType = fragmentType;
    }

    /**
     * Constructor that copies the residues and modifications from src.
     * The new PeptideFragment begins at the specified <code>beginIndex</code>
     * and extends to the residue at index <code>endIndex - 1</code>.
     *
     * Thus the length of the new PeptideFragment is <code>endIndex-beginIndex</code>.
     *
     * @param fragmentType the fragment type
     * @param src the PeptideFragment from which to copy the sequence and modifications
     * @param beginIndex the beginning index, inclusive.
     * @param endIndex   the ending index, exclusive.
     */
    public PeptideFragment(FragmentType fragmentType, AminoAcidSequence src, int beginIndex, int endIndex) {

        super(src, beginIndex, endIndex);

        checkNotNull(fragmentType);
        switch (fragmentType) {

            case FORWARD:

                checkArgument(0 == beginIndex, "Forward fragments have to start at 0, begin index was " + beginIndex);
                break;
            case REVERSE:

                checkArgument(src.size() == endIndex, "Reverse fragments have to end at " + src.size() + " , end index was " + endIndex);
                break;
            case MONOMER:

                checkArgument(endIndex - beginIndex == 1);
                break;
            case INTERNAL:
            case INTACT:
                break;
            default:
                throw new IllegalArgumentException("Unknown fragment type " + fragmentType);
        }

        this.fragmentType = fragmentType;
    }

    /**
     * Constructs an PeptideFragment from the list of residues with the modifications contained in
     * side chain and term mod map.
     *
     * @param fragmentType the fragment type
     * @param residues the residues, this list must have at least one amino acid
     * @param sideChainModMap the side chain modifications, this map can be empty
     * @param termModMap the terminal modifications, this map can be empty
     */
    public PeptideFragment(FragmentType fragmentType, List<AminoAcid> residues, Multimap<Integer, Modification> sideChainModMap, Multimap<ModAttachment, Modification> termModMap) {

        super(residues, sideChainModMap, termModMap);

        checkNotNull(fragmentType);
        this.fragmentType = fragmentType;
    }

    /**
     * Parse s to build a PeptideFragment.
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
     * @param fragmentType the fragment type
     * @return the PeptideFragment
     * @throws PeptideParseException if s is not a valid string
     */
    public static PeptideFragment parse(final String s, final FragmentType fragmentType) {

        List<AminoAcid> sequence = new ArrayList<AminoAcid>();
        ListMultimap<Integer, Modification> sideChainMatchMap = ArrayListMultimap.create();
        Multimap<ModAttachment, Modification> termModMap = ArrayListMultimap.create();

        PeptideParser.getInstance().parsePeptide(s, sequence, sideChainMatchMap, termModMap);

        return new PeptideFragment(fragmentType, sequence, sideChainMatchMap, termModMap);
    }

    /**
     * Parse <code>s</code> to build a PeptideFragment using the modResolver to convert modification strings to
     * Modification
     *
     * @param s the string
     * @param fragmentType the fragment type
     * @param modResolver Function to convert modification string to Modification
     * @return the Peptide
     * @throws PeptideParseException if s is not a valid string
     */
    public static PeptideFragment parse(final String s, final FragmentType fragmentType, ModificationResolver modResolver) {

        List<AminoAcid> sequence = new ArrayList<AminoAcid>();
        ListMultimap<Integer, Modification> sideChainMatchMap = ArrayListMultimap.create();
        Multimap<ModAttachment, Modification> termModMap = ArrayListMultimap.create();

        PeptideParser.getInstance().parsePeptide(s, sequence, sideChainMatchMap, termModMap, modResolver);

        return new PeptideFragment(fragmentType, sequence, sideChainMatchMap, termModMap);
    }

    /**
     * Return the fragmentType
     *
     * @return the fragmentType
     */
    public FragmentType getFragmentType() {

        return fragmentType;
    }

    /**
     * Calculate the m/z of this PeptideFragment given the <code>charge</code> and
     * <code>ionType</code>
     *
     * @param ionType the ionType
     * @param charge the charge
     * @return the m/z
     */
    public double calculateMz(IonType ionType, int charge) {

        checkArgument(ionType.getFragmentType() == fragmentType);
        checkNotNull(ionType);
        checkArgument(charge != 0, "Cannot calculate m/z, charge was 0");

        return AAMassCalculator.getInstance().calculateMz(calculateMass(ionType), charge, ionType);
    }

    /**
     * Calculate the mass of this PeptideFragment given the <code>ionType</code>
     *
     * @param ionType the ionType
     * @return the mass
     */
    public double calculateMass(IonType ionType) {

        checkNotNull(ionType);
        checkArgument(fragmentType.equals(ionType.getFragmentType()));

        double mass = calculateMonomerMassSum();

        mass += AAMassCalculator.getInstance().getDeltaMass(ionType);

        return mass;
    }

    @Override
    public boolean equals(Object o) {

        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        if (!super.equals(o)) return false;

        PeptideFragment fragment = (PeptideFragment) o;

        return fragmentType == fragment.fragmentType;

    }

    @Override
    public int hashCode() {

        int result = super.hashCode();
        result = 31 * result + fragmentType.hashCode();
        return result;
    }
}
