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


import newance.mzjava.ms.spectrum.AnnotatedPeak;
import newance.mzjava.ms.spectrum.PepFragAnnotation;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static com.google.common.base.Preconditions.checkArgument;
import static com.google.common.base.Preconditions.checkNotNull;

/**
 * @author Oliver Horlacher
 * @version 1.0
 */
public class BackbonePeakGenerator implements PeptidePeakGenerator<PepFragAnnotation> {

    private final Set<IonType> ionTypes;
    private final Set<FragmentType> fragmentTypes;
    private final double peakIntensity;

    public BackbonePeakGenerator(Set<IonType> ionTypes, double peakIntensity) {

        this.ionTypes = ionTypes;
        this.peakIntensity = peakIntensity;

        fragmentTypes = new HashSet<FragmentType>(ionTypes.size());
        for(IonType ionType : ionTypes) {

            fragmentTypes.add(ionType.getFragmentType());
        }
    }

    @Override
    public List<AnnotatedPeak<PepFragAnnotation>> generatePeaks(Peptide precursor, PeptideFragment fragment,
                                                                int[] charges, List<AnnotatedPeak<PepFragAnnotation>> peaks) {

        checkNotNull(precursor);
        checkNotNull(fragment);
        checkNotNull(charges);
        checkArgument(!fragment.isEmpty() , "fragment has to have a size that is > 0");
        checkArgument(charges.length > 0, "to generate peaks the charges array has to have at least one element");

        if(peaks == null) peaks = new ArrayList<AnnotatedPeak<PepFragAnnotation>>();

        for(IonType ionType : ionTypes) {

            if (fragment.getFragmentType() == ionType.getFragmentType()) {

                for (int charge : charges) {

                    double mz = fragment.calculateMz(ionType, charge);
                    PepFragAnnotation annotation = new PepFragAnnotation(ionType, charge, fragment);

                    AnnotatedPeak<PepFragAnnotation> peak = new AnnotatedPeak<PepFragAnnotation>(mz, peakIntensity, charge, annotation);
                    peaks.add(peak);
                }
            }
        }

        return peaks;
    }

    @Override
    public Set<FragmentType> getFragmentTypes() {

        return fragmentTypes;
    }
}
