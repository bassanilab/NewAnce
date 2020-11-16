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

import com.google.common.base.Preconditions;
import newance.mzjava.mol.modification.Tolerance;
import newance.mzjava.ms.peaklist.PeakAnnotation;
import newance.mzjava.ms.peaklist.PeakCursor;
import newance.mzjava.ms.peaklist.PeakList;

import java.util.Collections;
import java.util.List;

/**
 * Aligns two peak lists.  This aligner does not resolve many to many conflicts.
 *
 * To determine if there are any of these conflicts a OverlapListener can be set on the DefaultPeakListAligner.
 *
 * @author Oliver Horlacher
 * @version sqrt 1.0
 */
public class MatchedPeakListAligner<X extends PeakAnnotation, Y extends PeakAnnotation> implements PeakListAligner<X, Y> {

    private PeakPairSink<X, Y> sink;

    public interface OverlapListener {

        void overlap();
    }

    private final Tolerance tolerance;
    private OverlapListener overlapListener;
    double lastCentroid;

    public MatchedPeakListAligner(Tolerance error) {

        this.tolerance = error;
    }

    public MatchedPeakListAligner(Tolerance tolerance, PeakPairSink<X, Y> sink) {

        this(tolerance);

        this.sink = sink;
    }

    @Override
    public void align(PeakList<X> xPeakList, PeakList<Y> yPeakList) {

        Preconditions.checkNotNull(sink, "A sink needs to be set on this peak processor list");

        sink.begin(xPeakList, yPeakList);

        PeakCursor<X> xCursor = xPeakList.cursor();
        PeakCursor<Y> yCursor = yPeakList.cursor();

        double currentX, currentY;
        double current = 0;
        lastCentroid = - 1;
        boolean moreX, moreY;

        int maxIterations = xCursor.size() + yCursor.size();
        int iteration = 0;
        while (true) {

            if(iteration > maxIterations) throw new IllegalStateException("Too many iterations");

            iteration++;

            double max = tolerance.getMax(current);

            moreX = xCursor.next(max);
            moreY = yCursor.next(max);

            if (moreX && moreY) {
                currentX = xCursor.currMz();
                currentY = yCursor.currMz();

                current = (currentX <= currentY) ? currentX : currentY;

                double vX = 0.0;
                double vY = 0.0;
                List<X> annotationsX = Collections.emptyList();
                List<Y> annotationsY = Collections.emptyList();

                boolean xInTolerance = (current == currentX || tolerance.withinTolerance(currentX, current)) && bestMatch(currentX, xCursor, yCursor);
                if (xInTolerance) {
                    vX = xCursor.currIntensity();
                    annotationsX = xCursor.currAnnotations();
                }
                boolean yInTolerance = (current == currentY || tolerance.withinTolerance(currentY, current)) && (!xInTolerance || bestMatch(currentY, yCursor, xCursor));
                if (yInTolerance) {
                    vY = yCursor.currIntensity();
                    annotationsY = yCursor.currAnnotations();
                }

                if(xInTolerance) {

                    current = currentX;
                } else {

                    current = currentY;
                }

                if (xInTolerance && yInTolerance) report(current, vX, vY, annotationsX, annotationsY);
            } else {

                break;
            }
        }

        sink.end();
    }

    private boolean bestMatch(double current, PeakCursor<? extends PeakAnnotation> vectorA, PeakCursor<? extends PeakAnnotation> vectorB) {

        double d = Math.abs(current - vectorB.currMz());
        double max = tolerance.getMax(current);
        boolean best = true;

        if(vectorB.canPeek(1)) {

            double peek = vectorB.peekMz(1);
            if(peek <= max && Math.abs(current - peek) < d) {

                vectorA.previous();
                best = false;
            }
        }

        return best;
    }

    private void report(double centroid, double vX, double vY, List<X> annotationsX, List<Y> annotationsY) {

        if(tolerance.withinTolerance(centroid, lastCentroid)) {

            fireOverlap();
        }
        sink.processPeakPair(centroid, vX, vY, annotationsX, annotationsY);
        lastCentroid = centroid;
    }

    private void fireOverlap() {

        if(overlapListener != null) overlapListener.overlap();
    }

    public void setOverlapListener(OverlapListener overlapListener) {

        this.overlapListener = overlapListener;
    }

    @Override
    public void setSink(PeakPairSink<X, Y> sink) {

        this.sink = sink;
    }
}
