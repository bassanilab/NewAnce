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
package newance.mzjava.ms.peaklist.peakfilter;

import newance.mzjava.ms.peaklist.AbstractPeakProcessor;
import newance.mzjava.ms.peaklist.PeakAnnotation;

import java.util.List;

/**
 * Peak processor that removes any peaks that have an intensity that is < threshold
 *
 * @author Oliver Horlacher
 * @version 1.0
 */
public class ThresholdFilter<A extends PeakAnnotation> extends AbstractPeakProcessor<A, A> {

    //The threshold
    private final double threshold;
    private final Inequality inequality;

    public enum Inequality{GREATER, GREATER_EQUALS}

    public ThresholdFilter(double threshold) {

        this(threshold, Inequality.GREATER_EQUALS);
    }

    public ThresholdFilter(double threshold, Inequality inequality) {

        this.threshold = threshold;
        this.inequality = inequality;
    }

    @Override
    public void processPeak(double mz, double intensity, List<A> annotations) {

        boolean match;

        if (inequality == Inequality.GREATER)
            match = intensity > threshold;
        else if (inequality == Inequality.GREATER_EQUALS)
            match = intensity >= threshold;
        else
            throw new IllegalStateException("Cannot test for " + inequality);

        if (match)
            sink.processPeak(mz, intensity, annotations);
    }
}

