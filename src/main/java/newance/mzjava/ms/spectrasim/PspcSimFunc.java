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
 * Poisson shared peak count similarity function.  For each peak the counts from vector x are compared
 * to the counts in vector y using Poisson comparison of two counts.  If the two counts are significantly different
 * the peak with the smaller count is set to zero.
 * <p/>
 * The shared peak count is then calculated from
 *  2*k/(m + n)
 * Where k = the number of shared peaks, m = the number of peaks in x and n = the number of peaks in y
 *
 * <p/>
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class PspcSimFunc<X extends PeakAnnotation, Y extends PeakAnnotation> extends AbstractSimFunc<X, Y> {

    private double m = 0;
    private double n = 0;
    private double k = 0;

    private final int minPeakCount;

    public PspcSimFunc(int minPeakCount, Tolerance tolerance) {

        super(tolerance);
        this.minPeakCount = minPeakCount;
    }

    public PspcSimFunc(int minPeakCount, PeakListAligner<X, Y> peakListAligner, PeakPairProcessor<X, Y>... chain) {

        super(peakListAligner, chain);
        this.minPeakCount = minPeakCount;
    }

    @Override
    public double calcSimilarity(PeakList<X> plX, PeakList<Y> plY) {

        m = 0;
        n = 0;
        k = 0;

        vectorize(plX, plY);

        double spc = 2 * k / (m + n);

        return m + n < minPeakCount ? Double.NaN : spc;
    }

    public int getTotalPeakCount() {

        return (int)(m + n);
    }

    @Override
    public void processPeakPair(double centroid, double x, double y, List<X> xAnnotations, List<Y> yAnnotations) {

        if (PoissonUtils.isSignificantlyDifferent(x, y, 0.80)) {

            //If only y is set to 0, the calculated score may be altered if x was already 0
            if(x < y) x = 0;
            else y = 0;
        }

        if (x > 0 && y > 0) {

            m++;
            n++;
            k++;
        } else if (x > 0) {

            m++;
        } else if (y > 0) {

            n++;
        }
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
