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

/**
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class BinomialUtils {

    private BinomialUtils() {}

    /**
     * Computes the Wilson Interval
     * (see http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Wilson_score_interval)
     * Given the total number of events and the number of "correct" events, returns in a double-array
     * in the first component the center of the Wilson interval and in the second component the
     * width of the interval. alpha=95%.
     */
    public static double[] wilson(int total, int correct) {

        double z = 1.96;
        double p = (double) correct / total;
        double center = (p + 1 / 2.0 / total * z * z) / (1 + 1.0 / total * z * z);
        double d = z * Math.sqrt((p * (1 - p) + 1 / 4.0 / total * z * z) / total)
                / (1 + 1.0 / total * z * z);
        return new double[]{center, d};
    }

    /**
     * Calculates the 95% confidence interval for the ratio of the two menas.
     *
     * @param mean1 mean 1
     * @param mean2 mean 2
     * @return array containing the low and high values of the interval
     */
    public static double[] confidenceInterval95(int mean1, int mean2) {

        double[] w = wilson(mean1 + mean2, mean1);

        return new double[]{w[0] - w[1], w[0] + w[1]};
    }
}
