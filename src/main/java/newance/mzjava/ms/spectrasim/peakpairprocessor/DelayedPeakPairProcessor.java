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
package newance.mzjava.ms.spectrasim.peakpairprocessor;

import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TIntObjectHashMap;
import newance.mzjava.ms.peaklist.PeakAnnotation;
import newance.mzjava.ms.peaklist.PeakList;

import java.util.Collections;
import java.util.List;

/**
 * A peak pair processor that caches all the peaks that are pushed to it by the processPeakPair method.
 *
 * This can be used by subclasses that need all the peaks before the processing can begin, for
 * example see the .
 *
 * @author Oliver Horlacher
 * @version 1.0
 */
public abstract class DelayedPeakPairProcessor<X extends PeakAnnotation, Y extends PeakAnnotation> extends AbstractPeakPairProcessor<X, Y> {

    private final TDoubleArrayList centroids = new TDoubleArrayList();
    private final TDoubleArrayList vectorX = new TDoubleArrayList();
    private final TDoubleArrayList vectorY = new TDoubleArrayList();
    //The cached annotations
    private final TIntObjectHashMap<List<X>> xAnnotationMap = new TIntObjectHashMap<List<X>>();
    //The cached annotations
    private final TIntObjectHashMap<List<Y>> yAnnotationMap = new TIntObjectHashMap<List<Y>>();

    @Override
    public void processPeakPair(double centroid, double xIntensity, double yIntensity, List<X> xAnnotations, List<Y> yAnnotations) {

        centroids.add(centroid);
        vectorX.add(xIntensity);
        vectorY.add(yIntensity);

        if(!xAnnotations.isEmpty()) {
            xAnnotationMap.put(centroids.size() - 1, xAnnotations);
        }

        if(!yAnnotations.isEmpty()) {
            yAnnotationMap.put(centroids.size() - 1, yAnnotations);
        }
    }

    @Override
    public void begin(PeakList<X> xPeakList, PeakList<Y> yPeakList) {

        centroids.clear();
        vectorX.clear();
        vectorY.clear();

        sink.begin(xPeakList, yPeakList);
    }

    @Override
    public void end() {

        process(centroids, vectorX, vectorY, xAnnotationMap, yAnnotationMap);

        sink.end();
    }

    /**
     * Returns the cached annotations for the peak at index for the given annotationMap
     *
     * @param index the index
     * @param annotationMap the annotationMap
     * @return the cached annotations for the peak at index for the given annotationMap
     */
    protected <T extends PeakAnnotation> List<T> getAnnotations(int index, TIntObjectHashMap<List<T>> annotationMap) {

        List<T> annotations;
        if(annotationMap.contains(index))
            annotations = annotationMap.get(index);
        else
            annotations = Collections.emptyList();

        return annotations;
    }

    /**
     * Preform the processing
     *
     * @param centroids the centroid
     * @param vectorX the intensities of spectrum x
     * @param vectorY the intensities of spectrum y
     * @param xAnnotationMap the annotations of spectrum x peaks
     * @param yAnnotationMap the annotations of spectrum y peaks
     */
    protected abstract void process(TDoubleArrayList centroids, TDoubleArrayList vectorX, TDoubleArrayList vectorY,
                                    TIntObjectHashMap<List<X>> xAnnotationMap, TIntObjectHashMap<List<Y>> yAnnotationMap);
}
