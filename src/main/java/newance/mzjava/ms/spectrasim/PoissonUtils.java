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

import cern.jet.random.StudentT;
import cern.jet.random.engine.MersenneTwister;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;

import java.util.Arrays;

/**
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class PoissonUtils {

    private static final StudentT tDistribution = new StudentT(50, new MersenneTwister());

    private PoissonUtils() {}

    public static boolean isSignificantlyDifferent(double i, double j, double significance) {

        if (significance < 0 || significance > 1)
            throw new IllegalArgumentException("The significance must be a value between 0 and 1. Significance was " + significance);

        double d = testSignificance(i, j);

        return testSignificant(d, significance);
    }

    private static boolean testSignificant(double d, double significance) {

        double cdf = tDistribution.cdf(-d);

        return 2 * cdf < (1 - significance);
    }

    public static boolean isSignificantlyDifferent(double x1, double y1, double x2, double y2, double significance) {

        if (significance < 0 || significance > 1)
            throw new IllegalArgumentException("The significance must be a value between 0 and 1. Significance was " + significance);

        double d = testSignificance(x1, y1, x2, y1);

        return testSignificant(d, significance);
    }

    /**
     * Retruns d for a test of significance between two Poisson counts in samples of equal counts.
     * <p/>
     * see http://books.google.co.nz/books?id=WovYdIgj6HMC&pg=PA34&lpg=PA34&dq=Comparison+of+two+Poisson+counts&source=bl&ots=eazoFoJUBV&sig=XIEtCDGLlXL7NzZdQws_lu1Iky4&hl=en&ei=uJh0Tv_yB8bwsgaB2NnBCw&sa=X&oi=book_result&ct=result&resnum=6&ved=0CFMQ6AEwBQ#v=onepage&q&f=false
     * see http://www.biology.ed.ac.uk/research/groups/jdeacon/statistics/tress10.html#Poisson%20distribution%20for%20count%20data
     *
     * @param x1 mean of sample 1
     * @param x2 mean of sample 2
     * @return the d value
     */
    public static double testSignificance(double x1, double x2) {

        return (Math.abs(x1 - (x1 + x2) / 2) - 0.5) / Math.sqrt((x1 + x2) / 4);
    }

    /**
     * Retruns d for a test of significance between two Poisson counts where the sample sizes are unequal.
     * <p/>
     * see http://books.google.co.nz/books?id=WovYdIgj6HMC&pg=PA34&lpg=PA34&dq=Comparison+of+two+Poisson+counts&source=bl&ots=eazoFoJUBV&sig=XIEtCDGLlXL7NzZdQws_lu1Iky4&hl=en&ei=uJh0Tv_yB8bwsgaB2NnBCw&sa=X&oi=book_result&ct=result&resnum=6&ved=0CFMQ6AEwBQ#v=onepage&q&f=false
     * see http://www.biology.ed.ac.uk/research/groups/jdeacon/statistics/tress10.html#Poisson%20distribution%20for%20count%20data
     *
     * @param x1 mean of sample 1
     * @param y1 size of sample 1
     * @param x2 mean of sample 2
     * @param y2 size of sample 2
     * @return the d value
     */
    public static double testSignificance(double x1, double y1, double x2, double y2) {

        final double numarator = Math.abs(x1 - (x1 + x2) * (y1 / (y1 + y2))) - 0.5;
        final double denominator = Math.sqrt((x1 + x2) * (y1 / (y1 + y2)) * (y2 / (y1 + y2)));
        return numarator / denominator;
    }

    /**
     * Calculates the 95% confidence interval for the ratio of the two means.
     *
     * @param mean1 mean 1
     * @param mean2 mean 2
     * @return array containing the low and high values of the interval
     */
    public static double[] ratioConfidenceInterval95(int mean1, int mean2) {

        double[] w = BinomialUtils.wilson(mean1 + mean2, mean1);

        double l = w[0] - w[1];
        double h = w[0] + w[1];
        return new double[]{l / (1 - l), h / (1 - h)};
    }

    public static double[] hackRatioConfidenceInterval95(int m1, int m2) {

        double[] interval1 = PoissonUtils.confidenceInterval95(m1);
        double[] interval2 = PoissonUtils.confidenceInterval95(m2);

        double[] ratios = new double[]{
                interval1[0] / interval2[0], interval1[0] / interval2[1],
                interval1[1] / interval2[0], interval1[1] / interval2[1]
        };

        Arrays.sort(ratios);

        return new double[]{ratios[0], ratios[3]};
    }

    /**
     * See http://www.math.mcmaster.ca/peter/s743/poissonalpha.html.
     */
    public static double[] confidenceInterval95(int mean) {

        double lower = mean > 0 ? qchisq(0.025, 2 * mean) / 2 : 0;
        double upper = qchisq(0.975, 2 * (mean + 1)) / 2;

        return new double[]{lower, upper};
    }

    /**
     * Calculates a confidence interval for a poisson distribution.
     *
     * @param mean       the mean
     * @param confidence the confidence valid input is from 0 to 1
     * @return the confidence interval
     */
    public static double[] confidenceInterval(int mean, double confidence) {

        double low = (1 - confidence) / 2;
        double high = 1 - low;

        double lower = mean > 0 ? qchisq(low, 2 * mean) / 2 : 0;
        double upper = qchisq(high, 2 * (mean + 1)) / 2;

        return new double[]{lower, upper};
    }

    /**
     * See the R method by the same name
     */
    private static double qchisq(double p, int df) {

        ChiSquaredDistribution chiSquare = new ChiSquaredDistribution(df);
        return chiSquare.inverseCumulativeProbability(p);
    }
}
