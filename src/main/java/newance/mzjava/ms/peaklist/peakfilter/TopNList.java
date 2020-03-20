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

import newance.mzjava.ms.peaklist.PeakAnnotation;

import java.text.MessageFormat;
import java.util.*;

/**
 * List that retains the <code>size</code> peaks that have the largest intensities.
 *
 * The stored peaks can be accessed using the cursor methods.
 *
 * TopNList list = new TopNList(3);
 *
 * list.resetCursor();
 * while(list.hasNext()) {
 *
 *     double mz = list.currMz();
 *     double intensity = list.currIntensity();
 *     List annotations = list.currAnnotations();
 * }
 *
 * @author Oliver Horlacher
 * @version 1.0
 */
public class TopNList<A extends PeakAnnotation> {

    protected LinkedList<Entry<A>> list = new LinkedList<>();
    private final int numberOfPeaksToRetain;

    private final Cursor cursor = new Cursor();

    private final Comparator<Entry> valueComparator = new Comparator<Entry>() {

        public int compare(Entry e1, Entry e2) {

            return Double.compare(e2.getIntensity(), e1.getIntensity());
        }
    };

    private final Comparator<Entry> keyComparator = new Comparator<Entry>() {

        public int compare(Entry e1, Entry e2) {

            return Double.compare(e1.getMz(), e2.getMz());
        }
    };

    /**
     * Constructs a TopNList that retains the top <code>size</code> peaks.
     *
     * @param numberOfPeaksToRetain the number of peaks to retain
     */
    public TopNList(int numberOfPeaksToRetain) {

        this.numberOfPeaksToRetain = numberOfPeaksToRetain;
    }

    /**
     * Add a new peak
     *
     * @param mz the peaks m/z
     * @param intensity the peaks intensity
     * @param annotations the peaks annotations
     */
    public void add(double mz, double intensity, List<A> annotations) {

        if(list.size() < numberOfPeaksToRetain) {

            list.add(new Entry<>(mz, intensity, annotations));

            if(list.size() == numberOfPeaksToRetain) {

                Collections.sort(list, valueComparator);
            }
        } else {

            Entry last = list.getLast();

            if(last.getIntensity() < intensity) {

                Entry<A> current = list.removeLast();
                current.setValues(mz, intensity, annotations);

                ListIterator<Entry<A>> it = list.listIterator();

                boolean added = false;
                while (it.hasNext()) {

                    Entry entry = it.next();

                    if(entry.getIntensity() < intensity) {

                        it.previous();

                        added = true;
                        it.add(current);
                        break;
                    }
                }
                if(!added) list.add(current);
            }
        }
    }

    /**
     * Return the number of peaks that will be retained.
     *
     * @return the number of peaks that will be retained
     */
    public int getNumberOfPeaksToRetain() {

        return numberOfPeaksToRetain;
    }

    /**
     * Remove all peaks
     */
    public void clear() {

        list.clear();
    }

    /**
     * Returns <tt>true</tt> if this contains no peaks.
     *
     * @return <tt>true</tt> if this contains no peaks
     */
    public boolean isEmpty() {

        return list.isEmpty();
    }

    /**
     * Reset the cursor
     */
    public void resetCursor() {

        Collections.sort(list, keyComparator);

        cursor.reset();
    }

    /**
     * Advance the cursor to the next peak
     *
     * @return true if the cursor was advanced
     */
    public boolean next() {

        return cursor.next();
    }

    /**
     * Return the m/z of the current peak
     *
     * @return the m/z of the current peak
     */
    public double currMz() {

        return cursor.getCurrMz();
    }

    /**
     * Return the intensity of the current peak
     *
     * @return the intensity of the current peak
     */
    public double currIntensity() {

        return cursor.getCurrIntensity();
    }

    /**
     * Return the current annotation
     *
     * @return the current annotation
     */
    public List<A> currAnnotations() {

        return cursor.getCurrAnnotations();
    }

    /**
     * Returns the intensity of the last peak
     *
     * @return the intensity of the last peak
     */
    public double getLastIntensity() {

        return list.getLast().getIntensity();
    }

    private class Cursor {

        private ListIterator<Entry<A>> iterator;
        private Entry<A> current;

        void reset() {

            current = null;
            iterator = list.listIterator();
        }

        boolean next() {

            if(iterator.hasNext()) {

                current = iterator.next();
                return true;
            }
            return false;
        }

        double getCurrMz() {

            return current.getMz();
        }

        double getCurrIntensity() {

            return current.getIntensity();
        }

        List<A> getCurrAnnotations() {

            return current.getAnnotations();
        }
    }


    private static class Entry<A> {

        private double mz, intensity;
        private List<A> annotations;

        public Entry(double mz, double value, List<A> annotations) {

            setValues(mz, value, annotations);
        }

        double getMz() {

            return mz;
        }

        double getIntensity() {

            return intensity;
        }

        List<A> getAnnotations() {

            return annotations;
        }

        final void setValues(double mz, double intensity, List<A> annotations) {

            this.mz = mz;
            this.intensity = intensity;
            this.annotations = annotations;
        }

        @Override
        public String toString() {

            return MessageFormat.format("m/z = {0}, intensity = {1}", mz, intensity);
        }
    }
}

