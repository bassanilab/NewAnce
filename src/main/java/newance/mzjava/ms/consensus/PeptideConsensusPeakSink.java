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
import newance.mzjava.ms.peaklist.PeakSink;
import newance.mzjava.ms.spectrum.PepFragAnnotation;
import newance.mzjava.ms.spectrum.PepLibPeakAnnotation;

import java.util.List;

/**
 * Peak filter that retains peaks that are within tolerance of a peak in the theoretical spectrum or peaks
 * that have a mergedPeakCount that is >= absoluteMinPeakCount && >= peakFraction*consensusSize.
 *
 * @author Markus Muller
 * @version sqrt -1
 */
public class PeptideConsensusPeakSink implements PeakSink<PepLibPeakAnnotation> {

    private final int minPeakCount;
    private final ConsensusSpectrum<PepLibPeakAnnotation> consensus;

    public PeptideConsensusPeakSink(PeptideConsensusSpectrum consensus, double peakFraction, int minPeakCount) {

        this.consensus = consensus;

        int calculatedPeakCount = Math.max(minPeakCount, (int) Math.floor(consensus.nrOfMembers() * peakFraction));
        if (consensus.nrOfMembers()==1)
            calculatedPeakCount = 1;
        this.minPeakCount = calculatedPeakCount;
    }

    @Override
    public void start(int size) {

        consensus.ensureCapacity(size);
    }

    @Override
    public void processPeak(double mz, double intensity, List<PepLibPeakAnnotation> annotations) {

        int peakCount = annotations.get(0).getMergedPeakCount();
        Optional<PepFragAnnotation> annotation = annotations.get(0).getOptFragmentAnnotation();

        if (annotation.isPresent()) {
            consensus.add(mz, intensity, annotations);
        } else if (peakCount >= minPeakCount) {
            consensus.add(mz, intensity, annotations);
        }
    }

    @Override
    public void end() {

        consensus.trimToSize();
    }
}
