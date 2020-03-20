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
package newance.mzjava.ms.peaklist;

import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TIntObjectHashMap;

import java.util.List;

/**
 * A peak processor that caches all the peaks that are pushed to it by the processPeak method.
 * <p/>
 * This can be used by subclasses that need all the peaks before the processing can begin, for
 * example see the NthPeakNormalizer.
 *
 * @author Oliver Horlacher
 * @version 1.0
 */
public abstract class DelayedPeakProcessor<IN extends PeakAnnotation, OUT extends PeakAnnotation> extends AbstractPeakProcessor<IN, OUT> {

    //The cached m/z's
    private final TDoubleArrayList mzList = new TDoubleArrayList();
    //The cached intensities
    private final TDoubleArrayList intensityList = new TDoubleArrayList();
    //The cached annotations
    private final TIntObjectHashMap<List<IN>> annotationMap = new TIntObjectHashMap<List<IN>>();

    @Override
    public void processPeak(double mz, double intensity, List<IN> annotations) {

        mzList.add(mz);
        intensityList.add(intensity);
        if (!annotations.isEmpty()) {
            annotationMap.put(mzList.size() - 1, annotations);
        }
    }

    @Override
    public void start(int size) {

        mzList.reset();
        intensityList.reset();

        mzList.ensureCapacity(size);
        intensityList.ensureCapacity(size);
        annotationMap.clear();

        sink.start(size);
    }

    @Override
    public void end() {

        if(!mzList.isEmpty())
            processCached(mzList, intensityList, annotationMap);

        sink.end();
    }

    /**
     * This method is called when all the peaks have been cached
     *
     * @param mzList the list of all peak m/z
     * @param intensityList the list of all peak intensities
     * @param annotationMap peak annotations mapped by the peak index
     */
    protected abstract void processCached(TDoubleArrayList mzList, TDoubleArrayList intensityList, TIntObjectHashMap<List<IN>> annotationMap);
}
