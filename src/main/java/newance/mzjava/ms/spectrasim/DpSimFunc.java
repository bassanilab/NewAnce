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
package newance.mzjava.ms.spectrasim;

import newance.mzjava.mol.modification.Tolerance;
import newance.mzjava.ms.peaklist.PeakAnnotation;
import newance.mzjava.ms.peaklist.PeakList;
import newance.mzjava.ms.spectrasim.peakpairprocessor.PeakListAligner;
import newance.mzjava.ms.spectrasim.peakpairprocessor.PeakPairProcessor;

import java.util.List;

/**
 * Calculates the similarity between two PeakLists using the dot product (dp).
 *
 * A minimum peak count can be set so that NaN is returned for peak list pairs that do not
 * have enough peaks to calculate a reliable dp.
 *
 * @author Oliver Horlacher
 * @version 1.0
 */
public class DpSimFunc<X extends PeakAnnotation, Y extends PeakAnnotation> extends AbstractSimFunc<X, Y> {

    private double aDotB = 0;

    private int peakCount = 0;
    private final int minPeakCount;

    /**
     * Construct a NdpSimFunc that uses the default peak list aligner.
     *
     * @param minPeakCount the minimum number of peaks that are required in order to calculate the similarity
     * @param tolerance the tolerance to use when aligning the two PeakLists
     */
    public DpSimFunc(int minPeakCount, Tolerance tolerance) {

        super(tolerance);
        this.minPeakCount = minPeakCount;
    }

    /**
     * Construct a NdpSimFunc that uses the supplied PeakListAligner and PeakPairProcessor chain to process
     * and align the peak lists.
     *
     * @param minPeakCount the minimum number of peaks that are required in order to calculate the ndp
     * @param peakListAligner the object that is used to allign the PeakList pairs
     * @param chain PeakPairProcessor chain that is used to pre-process the aligned peaks before calculating the similarity
     */
    public DpSimFunc(int minPeakCount, PeakListAligner<X, Y> peakListAligner, PeakPairProcessor<X, Y>... chain) {

        super(peakListAligner, chain);
        this.minPeakCount = minPeakCount;
    }

    @Override
    public double calcSimilarity(PeakList<X> plX, PeakList<Y> plY) {

        aDotB = 0;
        peakCount = 0;

        vectorize(plX, plY);

        return peakCount < minPeakCount ? Double.NaN : aDotB;
    }

    @Override
    public int getTotalPeakCount() {

        return peakCount;
    }

    @Override
    public void processPeakPair(double centroid, double xIntensity, double yIntensity, List<X> xAnnotations, List<Y> yAnnotations) {

        if(xIntensity != 0) peakCount++;
        if(yIntensity != 0) peakCount++;

        aDotB += xIntensity * yIntensity;
    }

    @Override
    public double getBestScore() {

        return 1;
    }

    @Override
    public double getWorstScore() {

        return 0;
    }
}
