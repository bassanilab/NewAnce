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

import com.google.common.base.Preconditions;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;

/**
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class PeakListMerger<A extends PeakAnnotation> implements PeakSource<A> {

    private PeakSink<A> sink;
    private final PeakCursorComparator comparator = new PeakCursorComparator();

    @Override
    public <PS extends PeakSink<A>> PS setSink(PS sink) {

        Preconditions.checkNotNull(sink);
        this.sink = sink;
        return sink;
    }

    public void merge(Collection<? extends PeakList<A>> peakLists) {

        Preconditions.checkNotNull(peakLists);
        Preconditions.checkNotNull(sink, "Sink is null, use setSink() to set a sink to collect the merged peaks");

        PeakCursor[] peakCursors = new PeakCursor[peakLists.size()];

        int totalSize = 0;
        int size = 0;
        for(PeakList peakList : peakLists) {

            if (!peakList.isEmpty()) {

                totalSize += peakList.size();
                final PeakCursor cursor = peakList.cursor();
                cursor.next();
                peakCursors[size++] = cursor;
            }
        }

        doMerge(peakCursors, totalSize, size);
    }

    public void merge(PeakCursor<A> peakCursor1, PeakCursor<A> peakCursor2) {

        Preconditions.checkNotNull(peakCursor1);
        Preconditions.checkNotNull(peakCursor2);
        Preconditions.checkNotNull(sink, "Sink is null, use setSink() to set a sink to collect the merged peaks");

        int totalSize = peakCursor1.size() + peakCursor2.size();
        if(!peakCursor1.isEmpty() && !peakCursor2.isEmpty()) {

            peakCursor1.next();
            peakCursor2.next();
            doMerge(new PeakCursor[]{peakCursor1, peakCursor2}, totalSize, 2);
        } else if(!peakCursor1.isEmpty()) {

            copy(peakCursor1, totalSize);
        } else if(!peakCursor2.isEmpty()) {

            copy(peakCursor2, totalSize);
        }
    }

    private void copy(PeakCursor<A> cursor, int totalSize) {

        sink.start(totalSize);

        while (cursor.next()) {

            sink.processPeak(cursor.currMz(), cursor.currIntensity(), cursor.currAnnotations());
        }

        sink.end();
    }

    private void doMerge(PeakCursor[] peakCursors, int totalSize, int size) {

        Arrays.sort(peakCursors, 0, size, comparator);

        sink.start(totalSize);

        while (size > 0) {

            PeakCursor current = peakCursors[0];

            //noinspection unchecked
            sink.processPeak(current.currMz(), current.currIntensity(), current.currAnnotations());

            if(current.next()) {

                if(size > 1 && current.currMz() > peakCursors[1].currMz()) {

                    int index = Arrays.binarySearch(peakCursors, 1, size, current, comparator);

                    if (index < 0) index = -1 * (index + 1);

                    System.arraycopy(peakCursors, 1, peakCursors, 0, index - 1);
                    peakCursors[index - 1] = current;
                }
            } else {

                if(size > 1) System.arraycopy(peakCursors, 1, peakCursors, 0, size - 1);
                size--;
            }
        }

        sink.end();
    }

    private static class PeakCursorComparator implements Comparator<PeakCursor>, Serializable {

        @Override
        public int compare(PeakCursor peakCursor1, PeakCursor peakCursor2) {

            return Double.compare(peakCursor1.currMz(), peakCursor2.currMz());
        }
    }
}
