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

/**
 * Static methods to calculate the weighted Pearson's correlation coefficient
 *
 * @author Oliver Horlacher
 * @version 1.0
 */
public class WeightedPearsonsCorrelationCoeff {

    private WeightedPearsonsCorrelationCoeff() {}

    public static double calcR(TDoubleArrayList vX, TDoubleArrayList vY, TDoubleArrayList weights) {

        double mX = weightedMean(vX, weights);
        double mY = weightedMean(vY, weights);

        double covXY = weightedCovariance(vX, vY, weights, mX, mY);
        double covXX = weightedCovariance(vX, vX, weights, mX, mX);
        double covYY = weightedCovariance(vY, vY, weights, mY, mY);

        return covXY/Math.sqrt(covXX*covYY);
    }

    private static double weightedMean(TDoubleArrayList values, TDoubleArrayList weights) {

        double sumWX = 0;
        double sumW = 0;
        for(int i = 0; i < values.size(); i++) {

            double w = weights.get(i);
            sumWX += w * values.get(i);
            sumW += w;
        }

        return sumWX/sumW;
    }

    private static double weightedCovariance(TDoubleArrayList vX, TDoubleArrayList vY, TDoubleArrayList vW,
                                             double mX, double mY) {

        double sumW = 0;
        double sumN = 0;
        for(int i = 0; i < vX.size(); i++) {

            double w = vW.get(i);

            sumW += w;

            sumN += w*(vX.get(i) - mX)*(vY.get(i) - mY);
        }

        return sumN/sumW;
    }
}
