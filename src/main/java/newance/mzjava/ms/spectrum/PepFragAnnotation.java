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
package newance.mzjava.ms.spectrum;


import newance.mzjava.mol.*;
import newance.mzjava.ms.peaklist.PeakAnnotation;

import java.text.DecimalFormat;

import static com.google.common.base.Preconditions.checkArgument;
import static com.google.common.base.Preconditions.checkNotNull;

/**
 * This class is used to annotate peptide fragment peaks in a PeakList.
 * <p/>
 * PepFragAnnotation can store data on the ionType, charge, peptide fragment, isotpe composition
 * <p/>
 * The annotation
 *
 * @author Oliver Horlacher
 * @author nikitin
 * @version 1.0
 */
public class PepFragAnnotation implements PeakAnnotation {

    public static final Composition EMPTY_COMPOSITION = new Composition();

    private final IonType ionType;
    private final Composition isotopeComposition;
    private final int charge;
    private final PeptideFragment fragment;
    private final Mass neutralLoss;
    private double theoreticalMz;

    /**
     * Copy constructor
     *
     * @param src the PepFragAnnotation to copy
     */
    public PepFragAnnotation(PepFragAnnotation src) {

        ionType = src.ionType;
        isotopeComposition = src.isotopeComposition;
        charge = src.charge;
        fragment = src.fragment;
        neutralLoss = src.neutralLoss;
        theoreticalMz = src.theoreticalMz;
    }

    /**
     * Construct a PepFragAnnotation
     *
     * @param ionType  the ion type
     * @param charge   the charge
     * @param sequence the sequence of the
     */
    public PepFragAnnotation(IonType ionType, int charge, AminoAcidSequence sequence) {

        this(ionType, charge, sequence, EMPTY_COMPOSITION, Mass.ZERO);
    }

    private PepFragAnnotation(IonType ionType, int charge, AminoAcidSequence sequence, Composition isotopeComposition, Mass neutralLoss) {

        checkNotNull(ionType);
        checkNotNull(sequence);
        checkNotNull(neutralLoss);

        if (sequence instanceof PeptideFragment && ((PeptideFragment) sequence).getFragmentType() == ionType.getFragmentType()) {

            fragment = (PeptideFragment) sequence;
        } else {

            fragment = new PeptideFragment(ionType.getFragmentType(), sequence, 0, sequence.size());
        }
        this.isotopeComposition = isotopeComposition;
        this.ionType = ionType;
        this.neutralLoss = neutralLoss;
        this.charge = charge;
        this.theoreticalMz = -1;
    }

    /**
     * Return the peaks ion type
     *
     * @return the peaks ion type
     */
    public IonType getIonType() {

        return ionType;
    }

    /**
     * Returns true if this annotation annotates a peak with a neutral loss, false otherwise
     *
     * @return true if this annotation has a neutral loss, false otherwise
     */
    public boolean hasNeutralLoss() {

        return neutralLoss.getMolecularMass() != 0;
    }

    /**
     * Return the mass shift that is associated with this annotation
     *
     * @return the mass shift that is associated with this annotation
     */
    public Mass getNeutralLoss() {

        return neutralLoss;
    }

    /**
     * Return the peptide fragment associated with this annotation
     *
     * @return the peptide fragment associated with this annotation
     */
    public PeptideFragment getFragment() {

        return fragment;
    }

    /**
     * Return the number of isotopes
     *
     * @return the number of isotopes
     */
    public int getIsotopeCount() {

        return isotopeComposition.size();
    }

    /**
     * Return the Composition that contains the isotopes.
     *
     * @return the Composition that contains the isotopes
     */
    public Composition getIsotopeComposition() {

        return isotopeComposition;
    }

    public int getCharge() {

        return charge;
    }

    @Override
    public boolean equals(Object o) {

        if (this == o) return true;
        if (!(o instanceof PepFragAnnotation)) return false;

        PepFragAnnotation that = (PepFragAnnotation) o;

        return charge == that.charge &&
                isotopeComposition.equals(that.isotopeComposition) &&
                fragment.equals(that.fragment) &&
                ionType == that.ionType &&
                neutralLoss.equals(that.neutralLoss);
    }

    @Override
    public int hashCode() {

        int result = ionType.hashCode();
        result = 31 * result + isotopeComposition.hashCode();
        result = 31 * result + charge;
        result = 31 * result + fragment.hashCode();
        result = 31 * result + neutralLoss.hashCode();
        return result;
    }

    @Override
    public PepFragAnnotation copy() {

        return new PepFragAnnotation(this);
    }

    /**
     * Returns the theoretical m/z for this annotation.
     *
     * @return the theoretical m/z for this annotation
     */
    public double getTheoreticalMz() {

        if (theoreticalMz < 0)
            theoreticalMz = calcTheoreticalMz();

        return theoreticalMz;
    }

    private double calcTheoreticalMz() {

        double isotopeDelta;
        if (isotopeComposition.isEmpty()) {

            isotopeDelta = 0;
        } else {

            isotopeDelta = MassCalculator.calcIsotopeDelta(isotopeComposition);
        }

        return AAMassCalculator.getInstance().calculateMz(fragment.calculateMass(ionType) + isotopeDelta + neutralLoss.getMolecularMass(), charge, ionType);
    }

    /**
     * Convert this annotation to a SptxtString
     *
     * @return a string representation of this annotation in the Sptxt format
     */
    public String toSptxtString() {

        StringBuilder buff = new StringBuilder();

        buff.append(ionType);
        buff.append(fragment.size());
        for (int i = 0; i < isotopeComposition.size(); i++) {

            buff.append('i');
        }
        if (!neutralLoss.equals(Mass.ZERO)) {

            double mw = neutralLoss.getMolecularMass();
            if (mw > 0) buff.append('+');
            buff.append(DecimalFormat.getInstance().format(mw));
        }
        if (charge != 1) {

            buff.append('^');
            buff.append(charge);
        }

        return buff.toString();
    }

    @Override
    public String toString() {

        return "FragmentAnnotation{" +
                "ionType=" + ionType +
                ", isotopeComposition=" + isotopeComposition +
                ", charge=" + charge +
                ", fragment=" + fragment +
                ", neutralLoss=" + neutralLoss +
                '}';
    }

    /**
     * PepFragAnnotation Builder
     */
    public static class Builder {

        private static final PeriodicTable periodicTable = PeriodicTable.getInstance();

        private final IonType ionType;
        private int charge;
        private AminoAcidSequence fragment;
        private Mass neutralLoss = Mass.ZERO;

        private final Composition.Builder compositionBuilder = new Composition.Builder();

        /**
         * Construct a new Builder
         *
         * @param ionType  the ion type
         * @param charge   the charge
         * @param fragment the amino acid sequence of the fragment
         */
        public Builder(IonType ionType, int charge, AminoAcidSequence fragment) {

            this.ionType = ionType;
            this.charge = charge;
            this.fragment = fragment;
        }

        /**
         * Copy constructor
         *
         * @param src the PepFragAnnotation from which to copy the values
         */
        public Builder(PepFragAnnotation src) {

            this.ionType = src.getIonType();
            this.charge = src.getCharge();
            this.fragment = src.getFragment();
        }

        /**
         * Set the fragment
         *
         * @param fragment the fragment
         * @return this builder
         */
        public Builder setFragment(AminoAcidSequence fragment) {

            this.fragment = fragment;
            return this;
        }

        /**
         * Set the neutral loss
         *
         * @param neutralLoss the neutral loss
         * @return this builder
         */
        public Builder setNeutralLoss(Mass neutralLoss) {

            this.neutralLoss = neutralLoss;
            return this;
        }

        /**
         * Add <code>count</code> C13 isotopes
         *
         * @param count the number of C13 isotopes
         * @return this builder
         */
        public Builder addC13(int count) {

            compositionBuilder.add(periodicTable.getAtom(AtomicSymbol.C, 13), count);
            return this;
        }

        /**
         * Add <code>count</code> H2 isotpes
         *
         * @param count the number of H2 to add
         * @return this builder
         */
        public Builder addH2(int count) {

            compositionBuilder.add(periodicTable.getAtom(AtomicSymbol.H, 2), count);
            return this;
        }

        /**
         * Add the <code>isotope</code>
         *
         * @param isotope the isotope to add
         * @return this builder
         */
        public Builder addIsotope(Atom isotope) {

            checkArgument(!isotope.isDefaultIsotope());

            compositionBuilder.add(isotope, 1);
            return this;
        }

        /**
         * Add all the atoms from <code>isotopeComposition</code>
         *
         * @param isotopeComposition the atoms to add
         * @return this builder
         */
        public Builder addIsotopeComposition(Composition isotopeComposition) {

            compositionBuilder.addAll(isotopeComposition);
            return this;
        }

        /**
         * Build the PepFragAnnotation annotation.
         *
         * @return the new PepFragAnnotation
         */
        public PepFragAnnotation build() {

            Composition isotopeComposition = compositionBuilder.isEmpty() ? PepFragAnnotation.EMPTY_COMPOSITION : compositionBuilder.build();

            return new PepFragAnnotation(ionType, charge, fragment, isotopeComposition, neutralLoss);
        }

        public void setCharge(int charge) {

            this.charge = charge;
        }

        public int getCharge() {

            return charge;
        }
    }
}
