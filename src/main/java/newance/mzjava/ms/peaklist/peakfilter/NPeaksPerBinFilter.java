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

import java.util.ArrayList;
import java.util.List;

/**
 * Filter that retains the top k peaks per fixed mz range
 *
 * @author Oliver Horlacher
 * @version 1.0
 */
public class NPeaksPerBinFilter<A extends PeakAnnotation> extends AbstractPeakProcessor<A, A> {

    private final double range;
    private double max;

    private final int k;

    private final List<TopNList<A>> topNLists = new ArrayList<TopNList<A>>();
    private TopNList<A> currTopNList;

    /**
     *
     * @param k the number of peaks per mzWindow
     * @param mzWindow the mzWindow
     */
    public NPeaksPerBinFilter(int k, double mzWindow) {

        this.k = k;
        this.range = mzWindow;
    }

    @Override
    public void start(int size) {

        topNLists.clear();
        currTopNList = new TopNList<A>(k);

        max = range;

        sink.start(size);
    }

    @Override
    public void processPeak(double mz, double intensity, List<A> annotations) {

        if (mz >= max) {

            if(!currTopNList.isEmpty()) {

                topNLists.add(currTopNList);
                currTopNList = new TopNList<A>(k);
            }

            while(max < mz) max += range;
        }

        currTopNList.add(mz, intensity, annotations);
    }

    @Override
    public void end() {

        if(!currTopNList.isEmpty()) topNLists.add(currTopNList);

        for(TopNList<A> filter : topNLists) {

            filter.resetCursor();

            while(filter.next()) {

                sink.processPeak(filter.currMz(), filter.currIntensity(), filter.currAnnotations());
            }
        }

        sink.end();
    }
}

