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
package newance.mzjava.ms.consensus;

import com.google.common.base.Optional;
import newance.mzjava.mol.modification.Tolerance;
import newance.mzjava.ms.peaklist.AbstractPeakProcessor;
import newance.mzjava.ms.peaklist.PeakCursor;
import newance.mzjava.ms.spectrum.PepFragAnnotation;
import newance.mzjava.ms.spectrum.PepLibPeakAnnotation;
import newance.mzjava.ms.spectrum.PeptideSpectrum;

import java.util.ArrayList;
import java.util.List;

/**
 * Peak sink that collects the peaks into a LibrarySpectrum and annotates the peaks
 *
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class PeptideFragmentAnnotator extends AbstractPeakProcessor<PepLibPeakAnnotation, PepLibPeakAnnotation> {

    private final PeakCursor<PepFragAnnotation> theoreticalSpectrumCursor;
    private final Tolerance tolerance;


    public PeptideFragmentAnnotator(Tolerance tolerance, PeptideSpectrum theoreticalSpectrum) {

        this.tolerance = tolerance;
        this.theoreticalSpectrumCursor = theoreticalSpectrum.cursor();
    }

    @Override
    public void processPeak(double mz, double intensity, List<PepLibPeakAnnotation> annotations) {

        Tolerance.Location location;
        List<PepLibPeakAnnotation> newAnnotations = new ArrayList<>();
        do {

            if(!theoreticalSpectrumCursor.next()) break;
            location = tolerance.check(mz, theoreticalSpectrumCursor.currMz());
            if (location == Tolerance.Location.WITHIN) {
                // there is only one PepLibPeakAnnotation after merge
                PepLibPeakAnnotation libPeakAnnotation = annotations.get(0);
                List<PepFragAnnotation> fragAnnots = theoreticalSpectrumCursor.currAnnotations();
                for (PepFragAnnotation fragAnnot : fragAnnots) {
                    PepLibPeakAnnotation newAnnot =
                            new PepLibPeakAnnotation(
                                    libPeakAnnotation.getMergedPeakCount(),
                                    libPeakAnnotation.getMzStd(),
                                    libPeakAnnotation.getIntensityStd(),
                                    Optional.of(fragAnnot));
                    newAnnotations.add(newAnnot);
                }

            } else if(location == Tolerance.Location.LARGER){

                theoreticalSpectrumCursor.previous();
            }
        } while (location == Tolerance.Location.SMALLER);

        if (newAnnotations.size()>0) {
            sink.processPeak(mz, intensity, newAnnotations);
        }   else {
            sink.processPeak(mz, intensity, annotations);
        }
    }
}
