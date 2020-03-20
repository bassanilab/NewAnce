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

import java.util.List;

/**
 * A cursor for traversing a PeakList
 *
 * @author Oliver Horlacher
 * @version 1.0
 */
public interface PeakCursor<A extends PeakAnnotation> {

    /**
     * Returns the number of peaks int this cursor
     *
     * @return the number of peaks int this cursor
     */
    int size();

    /**
     * Returns <tt>true</tt> if this cursor contains no peaks.
     *
     * @return <tt>true</tt> if this cursor contains no peaks
     */
    boolean isEmpty();

    /**
     * Reset this cursor
     */
    void resetCursor();

    /**
     * Moves the cursor to the next peak if there are more peaks.  If this method returns false the
     * cursor is currently pointing to the last peak in the peak list
     *
     * @return true if the cursor was moved, false otherwise
     */
    boolean next();

    /**
     * Moves the cursor to the previous peak if there are any previous peaks.  If this method returns false the
     * cursor is currently pointing to the first peak in the peak list
     *
     * @return true if the cursor was moved, false otherwise
     */
    boolean previous();

    /**
     * Moves the cursor to the next peak only if peak.mz &lt;= <code>mz</code>.
     *
     * @param mz the threshold
     * @return true if there are more peaks, false otherwise
     */
    boolean next(double mz);

    /**
     * Returns the intensity of the peak that this cursor is currently pointing at
     *
     * @return the intensity of the peak that this cursor is currently pointing at
     */
    double currIntensity();

    /**
     * Returns the annotations of the peak that this cursor is currently pointing at
     *
     * @return the annotations of the peak that this cursor is currently pointing at
     */
    List<A> currAnnotations();

    /**
     * Returns the m/z of the peak that this cursor is currently pointing at
     *
     * @return the m/z of the peak that this cursor is currently pointing at
     */
    double currMz();

    /**
     * Returns true if there are n more peaks in the peak list
     *
     * @param n the number of peaks removed from the current peak
     * @return true if there are n more peaks in the peak list
     */
    boolean canPeek(int n);

    /**
     * Returns the intensity of the peak that is n positions removed from the peak that this cursor is currently pointing at
     *
     * @param n the number of peaks removed from the current peak
     * @return the intensity of the peak that is n positions removed from the peak that this cursor is currently pointing at
     */
    double peekIntensity(int n);

    /**
     * Returns the m/z of the peak that is n positions removed from the peak that this cursor is currently pointing at
     *
     * @param n the number of peaks removed from the current peak
     * @return the m/z of the peak that is n positions removed from the peak that this cursor is currently pointing at
     */
    double peekMz(int n);

    /**
     * Returns the m/z of the last peak
     *
     * @return the m/z of the last peak
     */
    double lastMz();

    /**
     * Returns the intensity of the last peak
     *
     * @return the intensity of the last peak
     */
    double lastIntensity();

    /**
     * Returns the m/z of the peak at index
     *
     * @param index the index of the peak
     * @return the m/z of the peak at index
     */
    double getMz(int index);

    /**
     * Returns the intensity of the peak at index
     *
     * @param index the index of the peak
     * @return the intensity of the peak at index
     */
    double getIntensity(int index);

    /**
     * Move the cursor to the m/z that is closest to <code>mz</code>.
     * If <code>mz</code> is exactly in the middle of two peaks the
     * cursor is moved to the peak with the smaller mz
     *
     * @param mz the m/z value to which the cursor is to be moved
     */
    void moveToClosest(double mz);

    /**
     * Move the cursor right before <code>mz</code> (peak.mz &lt; <code>mz</code>)
     *
     * @param mz the m/z value to which the cursor is to be moved before
     */
    void moveBefore(double mz);

    /**
     * Move the cursor right after <code>mz</code> (peak.mz &gt; <code>mz</code>)
     *
     * @param mz the m/z value to which the cursor is to be moved past
     */
    void movePast(double mz);

    /**
     * Returns the index of the peak which is closest to mz. If mz
     * is exactly in the middle of two peaks the index of the peak
     * with the smaller mz is returned.
     *
     * @param mz the m/z for which the closest index is to be found
     * @return the index of the peak which is closest to mz
     */
    int getClosestIndex(double mz);
}
