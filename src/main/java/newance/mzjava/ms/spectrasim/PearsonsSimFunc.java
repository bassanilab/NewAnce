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
import org.apache.commons.math3.stat.descriptive.UnivariateStatistic;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Variance;

import java.util.List;

/**
 * Calculates PeakList similarity using the Pearson's correlation coefficient.
 *
 * @author Oliver Horlacher
 * @version 1.0
 */
public class PearsonsSimFunc<X extends PeakAnnotation, Y extends PeakAnnotation> extends AbstractSimFunc<X, Y> {

    private final UnivariateStatistic mean = new Mean();
    private final UnivariateStatistic variance = new Variance();

    private double[] sX = new double[500];
    private int index = 0;
    private double[] sY = new double[500];

    private final int minPeakCount;
    private int peakCount = 0;

    /**
     * Construct a PearsonsSimFunc that uses the default peak list aligner.
     *
     * @param minPeakCount the minimum number of peaks that are required in order to calculate the similarity
     * @param tolerance the tolerance to use when aligning the two PeakLists
     */
    public PearsonsSimFunc(int minPeakCount, Tolerance tolerance) {

        super(tolerance);
        this.minPeakCount = minPeakCount;
    }

    /**
     * Construct a PearsonsSimFunc that uses the supplied PeakListAligner and PeakPairProcessor chain to process
     * and align the peak lists.
     *
     * @param minPeakCount the minimum number of peaks that are required in order to calculate the ndp
     * @param peakListAligner the object that is used to allign the PeakList pairs
     * @param chain PeakPairProcessor chain that is used to pre-process the aligned peaks before calculating the similarity
     */
    public PearsonsSimFunc(int minPeakCount, PeakListAligner<X, Y> peakListAligner, PeakPairProcessor<X, Y>... chain) {

        super(peakListAligner, chain);
        this.minPeakCount = minPeakCount;
    }

    @Override
    public double calcSimilarity(PeakList<X> plX, PeakList<Y> plY) {

        index = 0;
        peakCount = 0;

        vectorize(plX, plY);

        if(index < minPeakCount) return Double.NaN;

        double uX = apply(mean, sX);
        double uY = apply(mean, sY);

        double r = 0;

        double dxdy = getStandardDeviation(sX) * getStandardDeviation(sY);

        for (int i = 0; i < index; i++) {

            double x = sX[i];
            double y = sY[i];

            double v = ((x - uX) * (y - uY)) / dxdy;
            r = r + v;
        }

        return r/(index - 1);
    }

    @Override
    public int getTotalPeakCount() {

        return peakCount;
    }

    public double getStandardDeviation(double[] data) {

        double stdDev = Double.NaN;
        if (index > 0) {
            if (index > 1) {
                stdDev = Math.sqrt(apply(variance, data));
            } else {
                stdDev = 0.0;
            }
        }
        return stdDev;
    }

    public double apply(UnivariateStatistic stat, double[] data) {

        return stat.evaluate(data, 0, index);
    }

    @Override
    public void processPeakPair(double centroid, double xIntensity, double yIntensity, List<X> xAnnotations, List<Y> yAnnotations) {

        if(xIntensity > 0) peakCount++;
        if(yIntensity > 0) peakCount++;

        if(index >= sX.length) {

            sX = grow(sX);
            sY = grow(sY);
        }

        sX[index] = xIntensity;
        sY[index] = yIntensity;

        index++;
    }

    private double[] grow(double[] tmp) {

        double[] arr = new double[tmp.length + 100];
        System.arraycopy(tmp, 0, arr, 0, tmp.length);

        return arr;
    }

    @Override
    public double getBestScore() {

        return 1;
    }

    @Override
    public double getWorstScore() {

        return -1;
    }
}
