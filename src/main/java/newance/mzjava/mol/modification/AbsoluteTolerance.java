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
 * Implementation of Tolerance where the + and - is an absolute value
 *
 * Scaling to get around issues with representing some numbers as doubles. For
 * example on windows 125.3 - 125.4 = 0.10000000000000853 instead of the expected 0.1
 *
 * If the allowed error is +-0.1 then withinTolerance(125.3, 125.4) would return false
 * when it should return true
 *
 * @author Oliver Horlacher
 * @version 1.0
 */
public class AbsoluteTolerance extends Tolerance {

    private final double error;
    private final int scale = 100000;
    private final int scaledError;

    public AbsoluteTolerance(double error) {

        this.error = error;
        scaledError = (int)(error*scale);
    }

    public boolean withinTolerance(double expected, double actual) {

        return Math.abs((int)(expected *scale) - (int)(actual *scale)) <= scaledError;
    }

    @Override
    public double getError(double expectedMr, double actualMr) {

        return actualMr - expectedMr;
    }

    public Location check(double expected, double actual) {

        int scaledActual = (int)(actual * scale);

        int min = (int) (expected * scale) - scaledError;
        int max = (int) (expected * scale) + scaledError;

        if(scaledActual < min) return Location.SMALLER;
        else if(scaledActual > max) return Location.LARGER;
        else return Location.WITHIN;
    }

    public double getMin(double mz) {

        return mz - error;
    }

    public double getMax(double mz) {

        return mz + error;
    }
}
