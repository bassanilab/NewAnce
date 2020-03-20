/**
 * Copyright (c) 2011, SIB. All rights reserved.
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
package newance.mzjava.ms.io.mgf;


import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.Locale;

import static com.google.common.base.Preconditions.checkArgument;
import static com.google.common.base.Preconditions.checkNotNull;

/**
 * A factory to create Locale.US {@code NumberFormat}s.
 *
 * @author nikitin
 */
public final class NumberFormatFactory {

    /**
     * exponent value
     */
    private static final int DEFAULT_PRECISION = 2;

    private static final DecimalFormatSymbols US_SYMBOLS =
            new DecimalFormatSymbols(Locale.US);

    private static final NumberFormat DEFAULT_NUMBER_FORMAT =
            valueOf(DEFAULT_PRECISION);

    public static final NumberFormat INT_PRECISION =
            NumberFormatFactory.valueOf(1, 0);

    public static final NumberFormat FLOAT_PRECISION =
            NumberFormatFactory.valueOf(3);

    public static final NumberFormat DOUBLE_PRECISION =
            NumberFormatFactory.valueOf(10);

    private NumberFormatFactory() {

    }

    public static NumberFormat getInstance() {

        return DEFAULT_NUMBER_FORMAT;
    }

    public static NumberFormat valueOf(String pattern) {

        checkNotNull(pattern);
        return new DecimalFormat(pattern, US_SYMBOLS);
    }

    /**
     * Make a new instance of DecimalFormat with one maximum integer digit.
     *
     * @param maxFractionDigits the maximum number of fractional digits.
     * @return instance of DecimalFormat.
     */
    public static NumberFormat valueOf(int maxFractionDigits) {

        return valueOf(createDecimalPattern(1, maxFractionDigits));
    }

    /**
     * Make a new instance of DecimalFormat given the maximum integer and fraction digits.
     *
     * @param maxIntegerDigits  the maximum number of integer digits.
     * @param maxFractionDigits the maximum number of fractional digits.
     * @return instance of DecimalFormat.
     */
    public static NumberFormat valueOf(int maxIntegerDigits,
                                       int maxFractionDigits) {

        checkArgument(maxFractionDigits >= 0 && maxIntegerDigits >= 0);

        return valueOf(createDecimalPattern(maxIntegerDigits, maxFractionDigits));
    }

    /**
     * @return the defaultPrecision
     */
    public static int getDefaultPrecision() {

        return DEFAULT_PRECISION;
    }

    /**
     * Make a string format for decimal printing.
     *
     * @param maxIntegerDigits  the maximum number of integer digits
     * @param maxFractionDigits the maximum number of fractional digits
     * @return a string format.
     */
    private static String createDecimalPattern(int maxIntegerDigits, int maxFractionDigits) {

        StringBuilder sb = new StringBuilder();

        for (int i = 0; i < maxIntegerDigits; i++) {
            sb.append("0");
        }

        if (maxFractionDigits > 0) {
            sb.append(".");
        }

        for (int i = 0; i < maxFractionDigits; i++) {
            sb.append("0");
        }

        return sb.toString();
    }

}
