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
package newance.mzjava.mol.modification;

/**
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class PpmTolerance extends Tolerance {

    private final double errorPPM;

    public PpmTolerance(double errorPPM) {

        this.errorPPM = errorPPM;
    }

    @Override
    public boolean withinTolerance(double expected, double actual) {

        double error = calcError(expected);

        double min = expected - error;
        double max = expected + error;

        return actual >= min && actual <= max;
    }

    @Override
    public double getError(double expectedMr, double actualMr) {

        double delta = actualMr/expectedMr;
        if (actualMr > expectedMr) {
            return (delta/expectedMr) * 1000000;
        } else {
            return -(delta/expectedMr) * 1000000;
        }
    }

    @Override
    public Location check(double expected, double actual) {

        double error = calcError(expected);

        double min = expected - error;
        double max = expected + error;

        if(actual < min) return Location.SMALLER;
        else if(actual > max) return Location.LARGER;
        else return Location.WITHIN;
    }

    @Override
    public double getMin(double mz) {

        double error = calcError(mz);

        return mz - error;
    }

    @Override
    public double getMax(double mz) {

        double error = calcError(mz);

        return  mz + error;
    }

    private double calcError(double expectedMass) {

        return expectedMass * (errorPPM/1000000);
    }
}
