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

import gnu.trove.list.array.TDoubleArrayList;
import newance.mzjava.mol.modification.Tolerance;
import newance.mzjava.ms.peaklist.PeakAnnotation;
import newance.mzjava.ms.peaklist.PeakList;
import newance.mzjava.ms.spectrasim.peakpairprocessor.PeakListAligner;
import newance.mzjava.ms.spectrasim.peakpairprocessor.PeakPairProcessor;

import java.util.List;

/**
 * Calculates PeakList similarity using the weighted Pearson's correlation coefficient.
 *
 * @author Oliver Horlacher
 * @version 1.0
 */
public class WeightedPearsonsSimFunc<X extends PeakAnnotation, Y extends PeakAnnotation> extends AbstractSimFunc<X, Y> {

    private final SummingTDoubleArrayList vectorX = new SummingTDoubleArrayList();
    private final SummingTDoubleArrayList vectorY = new SummingTDoubleArrayList();
    private final TDoubleArrayList weights = new TDoubleArrayList();

    private int peakCount;

    public WeightedPearsonsSimFunc(Tolerance tolerance) {

        super(tolerance);
    }

    public WeightedPearsonsSimFunc(PeakListAligner<X, Y> peakListAligner, PeakPairProcessor<X, Y>... chain) {

        super(peakListAligner, chain);
    }

    @Override
    public double calcSimilarity(PeakList<X> plX, PeakList<Y> plY) {

        vectorX.clear();
        vectorY.clear();
        weights.clear();
        peakCount = 0;

        vectorize(plX, plY);

        if(vectorX.isEmpty()) return Double.NaN;

        //Calculate weights
        double summ = vectorX.getSumm() + vectorY.getSumm();

        for(int i = 0; i < vectorX.size(); i++) {

            double total = vectorX.get(i) + vectorY.get(i);

            weights.add(total/summ);
        }

        return WeightedPearsonsCorrelationCoeff.calcR(vectorX, vectorY, weights);
    }

    public int getTotalPeakCount() {

        return peakCount;
    }

    @Override
    public void processPeakPair(double centroid, double xIntensity, double yIntensity, List<X> xAnnotations, List<Y> yAnnotations) {

        if(xIntensity > 0) peakCount++;
        if(yIntensity > 0) peakCount++;

        vectorX.add(xIntensity);
        vectorY.add(yIntensity);
    }

    private static class SummingTDoubleArrayList extends TDoubleArrayList {

        private double summ = 0;

        public double getSumm() {

            return summ;
        }

        @Override
        public boolean add(double v) {

            summ += v;

            return super.add(v);
        }

        @Override
        public void clear() {

            summ = 0;

            super.clear();
        }

        @Override
        public boolean equals(Object other) {

            return super.equals(other);
        }

        @Override
        public int hashCode() {

            return super.hashCode();
        }
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
