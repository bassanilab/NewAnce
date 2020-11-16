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


import newance.mzjava.ms.peaklist.Peak;
import newance.mzjava.ms.peaklist.PeakList;
import newance.mzjava.ms.spectrum.AnnotatedPeak;
import newance.mzjava.ms.spectrum.PepFragAnnotation;
import newance.mzjava.ms.spectrum.PeptideSpectrum;

import java.io.Serializable;
import java.util.*;

/**
 * @author Oliver Horlacher
 * @version 1.0
 */
public class PeptideFragmenter {

    private final List<PeptidePeakGenerator<PepFragAnnotation>> peakGeneratorList = new ArrayList<PeptidePeakGenerator<PepFragAnnotation>>();
    private final PeakList.Precision precision;

    private final PeakComparator peakComparator = new PeakComparator();

    /**
     * Construct a CID peptide fragmenter that generates peaks for the supplied ion types.
     *
     * @param ionTypes the ion types to generate peaks for
     * @param precision the precision of the generated PeptideSpectrum
     */
    public PeptideFragmenter(Set<IonType> ionTypes, PeakList.Precision precision) {

        this.precision = precision;
        peakGeneratorList.add(new BackbonePeakGenerator(ionTypes, 1));
    }

    /**
     * Constructs a CID peptide fragmenter that generates peaks using the supplied peak generators
     *
     * @param peakGeneratorList the list of peak generators to use for generating peaks
     * @param precision the precision of the generated PeptideSpectrum
     */
    public PeptideFragmenter(List<PeptidePeakGenerator<PepFragAnnotation>> peakGeneratorList, PeakList.Precision precision) {

        this.precision = precision;
        this.peakGeneratorList.addAll(peakGeneratorList);
    }

    public PeptideSpectrum fragment(Peptide peptide, int precursorCharge) {

        int[] charges = new int[precursorCharge];
        for (int i = 0; i < charges.length; i++) {

            charges[i] = i + 1;
        }

        return fragment(peptide, precursorCharge, charges);
    }

    public PeptideSpectrum fragment(Peptide peptide, int precursorCharge, int[] fragmentCharges) {

        Set<FragmentType> fragmentTypes = new HashSet<FragmentType>();
        for(PeptidePeakGenerator<PepFragAnnotation> peakGenerator : peakGeneratorList) {

            fragmentTypes.addAll(peakGenerator.getFragmentTypes());
        }

        List<PeptideFragment> fragmentList = new ArrayList<PeptideFragment>();
        if (fragmentTypes.contains(FragmentType.FORWARD)) {
            generateForward(peptide, fragmentList);
        }
        if (fragmentTypes.contains(FragmentType.REVERSE)) {
            generateReverse(peptide, fragmentList);
        }
        if (fragmentTypes.contains(FragmentType.INTERNAL)) {
            throw new IllegalStateException("Creating internal ion fragments is not supported yet");
        }
        if (fragmentTypes.contains(FragmentType.MONOMER)) {
            generateImmonium(peptide, fragmentList);
        }
        if(fragmentTypes.contains(FragmentType.INTACT)) {
            fragmentList.add(peptide.createFragment(0, peptide.size(), FragmentType.INTACT));
        }

        return fragment(peptide, precursorCharge, fragmentList, fragmentCharges);
    }

    public PeptideSpectrum fragment(Peptide peptide, int precursorCharge, List<PeptideFragment> fragmentList, int[] fragmentCharges) {

        List<AnnotatedPeak<PepFragAnnotation>> peaks = new ArrayList<AnnotatedPeak<PepFragAnnotation>>();

        PeptideSpectrum spectrum = new PeptideSpectrum(peptide, precursorCharge, precision);

        for (PeptideFragment fragment : fragmentList) {

            for (PeptidePeakGenerator<PepFragAnnotation> peakGenerator : peakGeneratorList) {

                peakGenerator.generatePeaks(peptide, fragment, fragmentCharges, peaks);
            }
        }

        //Sort peaks
        Collections.sort(peaks, peakComparator);
        double lastMz = 0;
        for (AnnotatedPeak<PepFragAnnotation> peak : peaks) {

            //Doing this check because sometimes peaks that have exactly the same m/z
            //have the m/z stored in a double that is slightly different
            //to avoid having more than one peaks for these we check for the delta.
            double mz = peak.getMz();
            if(mz - lastMz < 0.000000000001){
                mz = lastMz;
            }
            spectrum.add(mz, peak.getIntensity(), peak.getAnnotations());
            lastMz = mz;
        }

        Peak precursor = spectrum.getPrecursor();
        precursor.setMzAndCharge(peptide.calculateMz(precursorCharge), precursorCharge);
        precursor.setIntensity(1);
        spectrum.trimToSize();
        return spectrum;
    }

    protected void generateForward(Peptide peptide, List<PeptideFragment> fragments) {

        for (int i = 1; i < peptide.size(); i++) {

            fragments.add(peptide.createFragment(0, i, FragmentType.FORWARD));
        }
    }

    protected void generateReverse(Peptide peptide, List<PeptideFragment> fragments) {

        for (int i = peptide.size() - 1; i > 0; i--) {

            fragments.add(peptide.createFragment(i, peptide.size(), FragmentType.REVERSE));
        }
    }

    protected void generateImmonium(Peptide peptide, List<PeptideFragment> fragments) {

        Set<AminoAcid> acids = new HashSet<AminoAcid>();

        for(int i = 0; i < peptide.size(); i++) {

            acids.add(peptide.getSymbol(i));
        }

        for(AminoAcid acid : acids) {

            PeptideFragment fragment = new PeptideFragment(FragmentType.MONOMER, acid);

            fragments.add(fragment);
        }
    }

    private static class PeakComparator implements Comparator<AnnotatedPeak>, Serializable {

        @Override
        public int compare(AnnotatedPeak p1, AnnotatedPeak p2) {

            return Double.compare(p1.getMz(), p2.getMz());
        }
    }
}
