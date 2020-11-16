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

import com.google.common.base.Predicate;
import newance.mzjava.mol.modification.ModAttachment;
import newance.mzjava.mol.modification.Modification;
import newance.mzjava.mol.modification.ModificationList;
import newance.mzjava.ms.spectrum.AnnotatedPeak;
import newance.mzjava.ms.spectrum.PepFragAnnotation;

import java.util.*;

import static com.google.common.base.Preconditions.checkArgument;
import static com.google.common.base.Preconditions.checkNotNull;

/**
 * Generates a neutral loss for a peptide fragment if the <code>modificationPredicate</code> accepts
 * the fragment.
 * <p>There are convenience constructors for generating neutral losses depending on the presence
 * of a set of amino acids. To create a neutral loss generator that depends on the presence of a modified
 * amino acid, such as a phosphorylated serine, the ModifiedAaPresencePredicate can be used</p>
 * <p>There are also static methods to generate neutral losses that are common. Such as water, ammonia and phosphate loss</p>
 *
 * @author Oliver Horlacher
 * @version 1.0
 */
public class PeptideNeutralLossPeakGenerator implements PeptidePeakGenerator<PepFragAnnotation> {

    private final Mass massShift;
    private final Predicate<PeptideFragment> modificationPredicate;
    private final double peakIntensity;
    private final Map<FragmentType, Set<IonType>> ionTypeMap;

    public PeptideNeutralLossPeakGenerator(Mass massShift, Set<AminoAcid> affectedAA, Set<IonType> ionTypes, double peakIntensity) {

        this(massShift, new AaPresencePredicate(affectedAA), ionTypes, peakIntensity);
    }

    public PeptideNeutralLossPeakGenerator(Mass massShift, Predicate<PeptideFragment> modificationPredicate, Set<IonType> ionTypes, double peakIntensity) {

        this.massShift = massShift;
        this.modificationPredicate = modificationPredicate;
        this.peakIntensity = peakIntensity;

        ionTypeMap = new HashMap<>();
        for (IonType ionType : ionTypes) {

            final FragmentType fragmentType = ionType.getFragmentType();
            Set<IonType> ionTypesForFragment = ionTypeMap.get(fragmentType);
            if (ionTypesForFragment == null) {

                ionTypesForFragment = new HashSet<>();
                ionTypeMap.put(fragmentType, ionTypesForFragment);
            }
            ionTypesForFragment.add(ionType);
        }
    }

    public static PeptideNeutralLossPeakGenerator newWaterLossGenerator(Set<IonType> ionTypes, double peakIntensity) {

        return new PeptideNeutralLossPeakGenerator(Composition.parseComposition("H-2O-1"),
                EnumSet.of(AminoAcid.S, AminoAcid.T, AminoAcid.D, AminoAcid.E),
                ionTypes, peakIntensity);
    }

    public static PeptideNeutralLossPeakGenerator newAmmoniaLossGenerator(Set<IonType> ionTypes, double peakIntensity) {

        return new PeptideNeutralLossPeakGenerator(Composition.parseComposition("H-3N-1"),
                EnumSet.of(AminoAcid.Q, AminoAcid.K, AminoAcid.R, AminoAcid.N),
                ionTypes, peakIntensity);
    }

    public static PeptideNeutralLossPeakGenerator newModificationLossGenerator(Mass loss, Set<AminoAcid> aminoAcids, Set<Modification> modifications, Set<ModAttachment> modAttachments, Set<IonType> ionTypes, double peakIntensity) {

        return new PeptideNeutralLossPeakGenerator(loss,
                new ModifiedAaPresencePredicate(aminoAcids, modifications, modAttachments),
                ionTypes, peakIntensity);
    }

    public static PeptideNeutralLossPeakGenerator newPhosphateLossGenerator(Set<IonType> ionTypes, double peakIntensity) {

        return new PeptideNeutralLossPeakGenerator(
                Composition.parseComposition("H-1P-1O-3"),
                new ModifiedAaPresencePredicate(EnumSet.of(AminoAcid.S, AminoAcid.T), Collections.singleton(Modification.parseModification("HPO3")), EnumSet.of(ModAttachment.SIDE_CHAIN)),
                ionTypes, peakIntensity);
    }

    public static PeptideNeutralLossPeakGenerator newPhosphoricAcidLossGenerator(Set<IonType> ionTypes, double peakIntensity) {

        return new PeptideNeutralLossPeakGenerator(Composition.parseComposition("H-3P-1O-4"),
                new ModifiedAaPresencePredicate(EnumSet.of(AminoAcid.S, AminoAcid.T), Collections.singleton(Modification.parseModification("HPO3")), EnumSet.of(ModAttachment.SIDE_CHAIN)),
                ionTypes, peakIntensity);
    }

    @Override
    public List<AnnotatedPeak<PepFragAnnotation>> generatePeaks(Peptide precursor, PeptideFragment fragment, int[] charges, List<AnnotatedPeak<PepFragAnnotation>> peaks) {

        checkNotNull(peaks);
        checkNotNull(precursor);
        checkNotNull(fragment);
        checkNotNull(charges);
        checkArgument(!fragment.isEmpty(), "fragment has to have a size that is > 0");
        checkArgument(charges.length > 0, "to generate peaks the charges array has to have at leas one element");

        Set<IonType> ionTypes = ionTypeMap.get(fragment.getFragmentType());
        if (ionTypes == null)
            return peaks;

        if (modificationPredicate.apply(fragment)) {

            generatePeaks(fragment, charges, peaks, ionTypes);
        }

        return peaks;
    }

    private void generatePeaks(PeptideFragment fragment, int[] charges, List<AnnotatedPeak<PepFragAnnotation>> peaks, Set<IonType> ionTypes) {

        for (IonType ionType : ionTypes) {

            for (int charge : charges) {

                double mz = AAMassCalculator.getInstance().calculateMz(fragment.calculateMass(ionType) + massShift.getMolecularMass(), charge, ionType);
                PepFragAnnotation annotation = new PepFragAnnotation.Builder(ionType, charge, fragment).setNeutralLoss(massShift).build();

                AnnotatedPeak<PepFragAnnotation> peak = new AnnotatedPeak<>(mz, peakIntensity, charge, annotation);
                peaks.add(peak);
            }
        }
    }

    @Override
    public Set<FragmentType> getFragmentTypes() {

        return ionTypeMap.keySet();
    }

    private static class AaPresencePredicate implements Predicate<PeptideFragment> {

        private final Set<AminoAcid> aminoAcids;

        private AaPresencePredicate(Set<AminoAcid> aminoAcids) {

            this.aminoAcids = aminoAcids;
        }

        @Override
        public boolean apply(PeptideFragment peptideFragment) {

            for (int i = 0; i < peptideFragment.size(); i++) {

                if (aminoAcids.contains(peptideFragment.getSymbol(i)))
                    return true;
            }

            return false;
        }
    }

    /**
     * Accepts a PeptideFragment if the fragment contains an amino acid in <code>aminoAcids</code>
     * that has a modification which is equal to the mass of one modification in <code>modifications</code>.
     *
     * <p>For example if the amino acid set is {S, T}  and the modification set is {HPO3:phosphate} then a
     * neutral loss is generated for PEPTS(HPO3)IDE, PEPTS(HPO3:phosphate)IDE and PEPTS(H1P1O3)IDE but
     * not PEPTS(79.9)IDE or PEPTS(79.9:phosphate)IDE</p>
     */
    public static class ModifiedAaPresencePredicate implements Predicate<PeptideFragment> {

        private final Set<AminoAcid> aminoAcids;
        private final Set<Mass> modificationMasses;
        private final Set<ModAttachment> modAttachments;

        private ModifiedAaPresencePredicate(Set<AminoAcid> aminoAcids, Set<Modification> modifications, Set<ModAttachment> modAttachments) {

            this.aminoAcids = aminoAcids;
            this.modificationMasses = new HashSet<>(modifications.size());
            for(Modification mod : modifications){

                modificationMasses.add(mod.getMass());
            }
            this.modAttachments = modAttachments;
        }

        @Override
        public boolean apply(PeptideFragment fragment) {

            for(int i : fragment.getModificationIndexes(modAttachments)) {

                if(aminoAcids.contains(fragment.getSymbol(i)) && hasModifications(fragment.getModificationsAt(i, modAttachments)))
                    return true;
            }

            return false;
        }

        private boolean hasModifications(ModificationList modificationsAt) {

            for(Modification mod : modificationsAt) {

                if(this.modificationMasses.contains(mod.getMass()))
                    return true;
            }

            return false;
        }
    }
}
