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

import com.google.common.base.Function;
import com.google.common.base.Optional;

import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.UUID;

/**
 * A list of mz/intensity values.
 *
 * @author Oliver Horlacher
 * @version 1.0
 */
public interface PeakList<A extends PeakAnnotation>  {

    public enum Precision {

        DOUBLE,             //PeakList that has both m/z and intensities as double
        FLOAT,              //PeakList that has both m/z and intensities as float
        DOUBLE_FLOAT,       //PeakList that has m/z as double and intensities as float
        DOUBLE_CONSTANT,    //PeakList that has m/z as double and the intensities as a constant value
        FLOAT_CONSTANT      //PeakList that has m/z as float and the intensities as a constant value
    }

    UUID getId();

    void setId(UUID id);

    /**
     * Add peaks and annotations in <code>peakList</code> to this peak list.
     *
     * @param peakList containing peaks and annotations which are to be added to this peak list
     */
    void addPeaks(PeakList<A> peakList);

    /**
     * Add peaks (mz, intensity pairs) and optional annotations in <code>peakList</code> to this peak list
     *
     * @param peakList            containing peaks which are to be added to this peak list
     * @param annotationConverter convert annotation T to A (annotation will not be added if
     *                            Optional<A> is absent)
     */
    <T extends PeakAnnotation> void addPeaks(PeakList<T> peakList, Function<List<T>, List<A>> annotationConverter);

    /**
     * Add peaks (without annotations) in <code>peakList</code> to this peak list.
     *
     * @param peakList containing peaks which are to be added to this peak list
     */
    <T extends PeakAnnotation> void addPeaksNoAnnotations(PeakList<T> peakList);

    /**
     * Add the pairs of mz and intensities to this peak list. The mzs array has to be sorted.
     * If the m/z array is not sorted an IllegalArgumentException will be thrown
     *
     * @param mzs         the mz values
     * @param intensities the intensities
     * @throws ArrayIndexOutOfBoundsException                               if length is < 0 or larger than the mz or intensity array
     * @throws org.expasy.mzjava.core.ms.peaklist.UnsortedPeakListException if the mzs are not sorted
     */
    void addSorted(double[] mzs, double[] intensities);

    /**
     * Add the pairs of mz and intensities to this peak list. The mzs array has to be sorted.
     * If the m/z array is not sorted an IllegalArgumentException will be thrown
     *
     * @param mzs         the mz values
     * @param intensities the intensities
     * @param length      the number of pairs to add
     * @throws ArrayIndexOutOfBoundsException                               if length is < 0 or larger than the mz or intensity array
     * @throws org.expasy.mzjava.core.ms.peaklist.UnsortedPeakListException if the mzs are not sorted
     */
    void addSorted(double[] mzs, double[] intensities, int length);

    /**
     * Add a mz, intensity pair to the peak list
     *
     * @param mz        the mz to add
     * @param intensity the intensity to add
     * @return the index at which the peak was added
     */
    int add(double mz, double intensity);

    /**
     * Add a mz, intensity pair and associateed annotations to the peak list
     *
     * @param mz         the mz to add
     * @param intensity  the intensity to add
     * @param annotation the peaks annotation
     * @return the index at which the peak was added
     */
    int add(double mz, double intensity, A annotation);

    /**
     * Add a mz, intensity pair and associated annotations to the peak list
     *
     * @param mz          the mz to add
     * @param intensity   the intensity to add
     * @param annotations the peaks annotations
     * @return the index at which the peak was added
     */
    int add(double mz, double intensity, Collection<? extends A> annotations);

    /**
     * Returns the number of peaks in this peak list.
     *
     * @return the number of peaks in this peak list
     */
    int size();

    /**
     * Returns <tt>true</tt> if this peak list contains no peaks.
     *
     * @return <tt>true</tt> if this peak list contains no peaks
     */
    boolean isEmpty();

    /**
     * Returns the mz of the peak at index
     *
     * @param index the index
     * @return the mz of the peak at index
     */
    double getMz(int index);

    /**
     * Returns the intensity of the peak at index
     *
     * @param index the index
     * @return the intensity of the peak at index
     */
    double getIntensity(int index);

    /**
     * Replace the intensity at index with intensity
     *
     * @param intensity the new intensity
     * @param index     the index for which the intensity is to be set
     * @throws IndexOutOfBoundsException if the index is out of range
     *                                   (<tt>index &lt; 0 || index &gt;= size()</tt>)
     */
    public void setIntensityAt(double intensity, int index);

    /**
     * Clears the peak list, the size is set to 0 and objects used to store peaks are unreferenced to allow
     * garbage collection
     */
    void clear();

    /**
     * Returns the sum of all intensities
     *
     * @return the sum of all intensities
     */
    double getTotalIonCurrent();

    /**
     * Returns the index of the peak that is closest to the supplied mz.
     * If mz is exactly equidistant between two peaks the index of the smaller peak is returned
     *
     * @param mz the mz
     * @return the index of the peak that is closest to the supplied mz
     */
    int getClosestIndex(double mz);

    /**
     * Returns the index of the first occurrence of the specified mzKey
     * in this peak list, or (-(<i>insertion point</i>) - 1) if this
     * peak list does not contain the exact mzKey. The
     * <i>insertion point</i> is defined as the point at which the
     * key would be inserted into the peak list: the index of the first
     * mzKey greater than the key, or <tt>size()</tt> if all
     * m/z's in the peak list are less than the specified key.  Note
     * that this guarantees that the return value will be &gt;= 0 if
     * and only if the key is found.
     *
     * @param mzKey the search key
     * @return index of the search key, if it is contained in the peak list;
     * otherwise, <tt>(-(<i>insertion point</i>) - 1)</tt>.  The
     * <i>insertion point</i> is defined as the point at which the
     * key would be inserted into the peak list: the index of the first
     * mzKey greater than the key, or <tt>size()</tt> if all
     * m/z's in the peak list are less than the specified key.  Note
     * that this guarantees that the return value will be &gt;= 0 if
     * and only if the key is found.
     */
    int indexOf(double mzKey);

    /**
     * Returns the index of the first occurrence of the specified mzKey
     * in the specified range within this peak list, or (-(<i>insertion point</i>) - 1)
     * if this peak list does not contain the exact mzKey. The
     * <i>insertion point</i> is defined as the point at which the
     * key would be inserted into the peak list: the index of the first
     * mzKey greater than the key, or <tt>size()</tt> if all
     * m/z's in the peak list are less than the specified key.  Note
     * that this guarantees that the return value will be &gt;= 0 if
     * and only if the key is found.
     *
     * @param fromIndex the index of the first peak (inclusive) to be
     *                  searched
     * @param toIndex   the index of the last peak (exclusive) to be searched
     * @param mzKey     the search key
     * @return index of the search key, if it is contained in the peak list;
     * otherwise, <tt>(-(<i>insertion point</i>) - 1)</tt>.  The
     * <i>insertion point</i> is defined as the point at which the
     * key would be inserted into the peak list: the index of the first
     * mzKey greater than the key, or <tt>size()</tt> if all
     * m/z's in the peak list are less than the specified key.  Note
     * that this guarantees that the return value will be &gt;= 0 if
     * and only if the key is found.
     */
    int indexOf(int fromIndex, int toIndex, double mzKey);

    /**
     * Returns te index of the peak that is most intense in the specified region.
     * Where minMz <= region <= maxMZ.
     *
     * @param minMz the minMz
     * @param maxMz the maxMz
     * @return the index of the most intense peak or -1 if there is no such peak
     */
    int getMostIntenseIndex(double minMz, double maxMz);

    /**
     * Returns the index of the most intense peak in the peak list or -1 if there are no peaks in this peak list.
     *
     * @return the index of the most intense peak in the peak list
     */
    int getMostIntenseIndex();

    /**
     * Get the m/z of the most intense peak
     *
     * @return the m/z of the most intense peak
     */
    double getBasePeakMz();

    /**
     * Get the intensity of the most intense peak
     *
     * @return the intensity of the most intense peak
     */
    double getBasePeakIntensity();

    /**
     * Copies an the mzs from this peak list, beginning at 0,
     * to the destination array.
     * The number of mzs copied is equal to the size of this peak list.
     * The mzs at positions 0 through size() in this peak list are copied into
     * positions 0 through size(), respectively, of the destination array.
     * <p/>
     * If <code>dest</code> is <code>null</code>, then a
     * new array with the correct length is created.
     * <p/>
     * Otherwise, if any of the following is true, an
     * <code>IndexOutOfBoundsException</code> is
     * thrown and the destination is not modified:
     * <ul>
     * <li><code>size()</code> is greater than
     * <code>dest.length</code>, the length of the destination array.
     * </ul>
     * <p/>
     *
     * @param dest the destination array.
     * @return dest or the newly created array if dest was null
     * @throws IndexOutOfBoundsException if copying would cause
     *                                   access of data outside array bounds.
     */
    double[] getMzs(double[] dest);

    /**
     * Copies an the mzs from this peak list, beginning at 0,
     * to the specified position of the destination array.
     * The number of mzs copied is equal to the size of this peak list.
     * The mzs at positions 0 through size() in this peak list are copied into
     * positions <code>destPos</code> through <code>destPos+length-1</code>,
     * respectively, of the destination array.
     * <p/>
     * If <code>dest</code> is <code>null</code>, then a
     * new array with the correct length is created.
     * <p/>
     * Otherwise, if any of the following is true, an
     * <code>IndexOutOfBoundsException</code> is
     * thrown and the destination is not modified:
     * <ul>
     * <li>The <code>destPos</code> argument is negative.
     * <li><code>destPos+size()</code> is greater than
     * <code>dest.length</code>, the length of the destination array.
     * </ul>
     * <p/>
     *
     * @param dest    the destination array.
     * @param destPos starting position in the destination data.
     * @return dest or the newly created array if dest was null
     * @throws IndexOutOfBoundsException if copying would cause
     *                                   access of data outside array bounds.
     */
    double[] getMzs(double[] dest, int destPos);

    /**
     * Copies an the mzs from this peak list, beginning at the specified position,
     * to the specified position of the destination array.
     * A sub-sequence of mzs are copied from this peak list to the destination array
     * referenced by <code>dest</code>. The number of mzs copied is
     * equal to the <code>length</code> argument. The mzs at
     * positions <code>srcPos</code> through
     * <code>srcPos+length-1</code> in this peak list are copied into
     * positions <code>destPos</code> through
     * <code>destPos+length-1</code>, respectively, of the destination
     * array.
     * <p/>
     * If <code>dest</code> is <code>null</code>, then a
     * new array with the correct length is created.
     * <p/>
     * Otherwise, if any of the following is true, an
     * <code>IndexOutOfBoundsException</code> is
     * thrown and the destination is not modified:
     * <ul>
     * <li>The <code>srcPos</code> argument is negative.
     * <li>The <code>destPos</code> argument is negative.
     * <li>The <code>length</code> argument is negative.
     * <li><code>srcPos+length</code> is greater than
     * <code>size()</code>, the size of this.
     * <li><code>destPos+length</code> is greater than
     * <code>dest.length</code>, the length of the destination array.
     * </ul>
     * <p/>
     *
     * @param srcPos  starting position in the mz.
     * @param dest    the destination array.
     * @param destPos starting position in the destination data.
     * @param length  the number of array elements to be copied.
     * @return dest or the newly created array if dest was null
     * @throws IndexOutOfBoundsException if copying would cause
     *                                   access of data outside array bounds.
     */
    double[] getMzs(int srcPos, double[] dest, int destPos, int length);

    /**
     * Copies an the intensities from this peak list, beginning at 0,
     * to the destination array.
     * The number of intensities copied is equal to the size of this peak list.
     * The intensities at positions 0 through size() in this peak list are copied into
     * positions 0 through size(), respectively, of the destination array.
     * <p/>
     * If <code>dest</code> is <code>null</code>, then a
     * new array with the correct length is created.
     * <p/>
     * Otherwise, if any of the following is true, an
     * <code>IndexOutOfBoundsException</code> is
     * thrown and the destination is not modified:
     * <ul>
     * <li><code>size()</code> is greater than
     * <code>dest.length</code>, the length of the destination array.
     * </ul>
     * <p/>
     *
     * @param dest the destination array.
     * @return dest or the newly created array if dest was null
     * @throws IndexOutOfBoundsException if copying would cause
     *                                   access of data outside array bounds.
     */
    double[] getIntensities(double[] dest);

    /**
     * Copies an the intensities from this peak list, beginning at 0,
     * to the specified position of the destination array.
     * The number of intensities copied is equal to the size of this peak list.
     * The intensities at positions 0 through size() in this peak list are copied into
     * positions <code>destPos</code> through <code>destPos+length-1</code>,
     * respectively, of the destination array.
     * <p/>
     * If <code>dest</code> is <code>null</code>, then a
     * new array with the correct length is created.
     * <p/>
     * Otherwise, if any of the following is true, an
     * <code>IndexOutOfBoundsException</code> is
     * thrown and the destination is not modified:
     * <ul>
     * <li>The <code>destPos</code> argument is negative.
     * <li><code>destPos+size()</code> is greater than
     * <code>dest.length</code>, the length of the destination array.
     * </ul>
     * <p/>
     *
     * @param dest    the destination array.
     * @param destPos starting position in the destination data.
     * @return dest or the newly created array if dest was null
     * @throws IndexOutOfBoundsException if copying would cause
     *                                   access of data outside array bounds.
     */
    double[] getIntensities(double[] dest, int destPos);

    /**
     * Copies an the intensities from this peak list, beginning at the specified position,
     * to the specified position of the destination array.
     * A sub-sequence of intensities are copied from this peak list to the destination array
     * referenced by <code>dest</code>. The number of intensities copied is
     * equal to the <code>length</code> argument. The intensities at
     * positions <code>srcPos</code> through
     * <code>srcPos+length-1</code> in this peak list are copied into
     * positions <code>destPos</code> through
     * <code>destPos+length-1</code>, respectively, of the destination
     * array.
     * <p/>
     * If <code>dest</code> is <code>null</code>, then a
     * new array with the correct length is created.
     * <p/>
     * Otherwise, if any of the following is true, an
     * <code>IndexOutOfBoundsException</code> is
     * thrown and the destination is not modified:
     * <ul>
     * <li>The <code>srcPos</code> argument is negative.
     * <li>The <code>destPos</code> argument is negative.
     * <li>The <code>length</code> argument is negative.
     * <li><code>srcPos+length</code> is greater than
     * <code>size()</code>, the size of this.
     * <li><code>destPos+length</code> is greater than
     * <code>dest.length</code>, the length of the destination array.
     * </ul>
     * <p/>
     *
     * @param srcPos  starting position in the mz.
     * @param dest    the destination array.
     * @param destPos starting position in the destination data.
     * @param length  the number of array elements to be copied.
     * @return dest or the newly created array if dest was null
     * @throws IndexOutOfBoundsException if copying would cause
     *                                   access of data outside array bounds.
     */
    double[] getIntensities(int srcPos, double[] dest, int destPos, int length);

    /**
     * Trims the capacity of this <tt>DoublePeakList</tt> instance to be the
     * list's current size.  An application can use this operation to minimize
     * the storage of a <tt>DoublePeakList</tt> instance.
     */
    public void trimToSize();

    /**
     * Increases the capacity of this <tt>DoublePeakList</tt> instance, if
     * necessary, to ensure that it can hold at least the number of peaks
     * specified by the minimum capacity argument.
     *
     * @param minCapacity the desired minimum capacity
     */
    public void ensureCapacity(int minCapacity);

    /**
     * Return the precision that this PeakList uses to store the peak and intensity values
     *
     * @return the precision that this PeakList uses to store the peak and intensity values
     */
    public Precision getPrecision();

    /**
     * Add an annotation to the peak at index
     *
     * @param index      the index of the peak
     * @param annotation the annotation to add
     */
    void addAnnotation(int index, A annotation);

    /**
     * Set the annotations for the peak at index
     *
     * @param index       the index of the peak
     * @param annotations the list of annotations
     */
    void addAnnotations(int index, Collection<A> annotations);

    /**
     * Remove the annotation from the peak at index
     *
     * @param annotation the annotation to remove
     * @param index      the index of the peak
     * @return <tt>true</tt> if an annotation was removed as a result of this call
     */
    boolean removeAnnotation(A annotation, int index);

    /**
     * Remove all annotations for the peak at index
     *
     * @param index the index of the peak
     */
    void clearAnnotationsAt(int index);

    /**
     * Remove all annotations from this PeakList
     */
    void clearAnnotations();

    /**
     * Returns true if the spectrum has annotations
     *
     * @return true if the spectrum has annotations
     */
    boolean hasAnnotations();

    /**
     * Returns an array containing the sorted indexes of the peaks that are annotated.
     * <p/>
     * If this peak list has no annotated peaks an clear array is returned
     *
     * @return an array containing the sorted indexes of the peaks that are annotated
     */
    int[] getAnnotationIndexes();

    /**
     * Returns true if the spectrum has annotations at peak index
     *
     * @param index the index of the peak
     * @return true if the spectrum has annotations at peak index
     */
    boolean hasAnnotationsAt(int index);

    /**
     * Returns an unmodifiable list containing the annotations for the peak at index
     *
     * @param index the index of the peak
     * @return an unmodifiable list containing the annotations for the peak at index
     */
    List<A> getAnnotations(int index);

    /**
     * Return the first annotation in the list of annotations for the peak at index.
     * If the specified index does not have any annotations null is returned.
     *
     * @param index the index of the peak
     * @return the first annotation in the list of annotations, or null if the specified index has no annotations
     */
    Optional<A> getFirstAnnotation(int index);

    /**
     * Sort the annotations for each annotated peak using the <code>comparator</code>
     *
     * @param comparator the comparator that is to be used to sort the annotations
     */
    void sortAnnotations(Comparator<A> comparator);

    /**
     * Return PeakList copy processed with PeakProcessor
     * @param peakProcessor the processor to process the peaks
     * @return the copy of this with the peaks processed by peakProcessor
     */
    PeakList<A> copy(PeakProcessor<A, A> peakProcessor);

    /**
     * Return PeakList copy processed with PeakProcessorChain
     * @param peakProcessorChain the processor chain to process the peaks
     * @return the copy of this with the peaks processed by peakProcessorChain
     */
    PeakList<A> copy(PeakProcessorChain<A> peakProcessorChain);

    /**
     * Returns the precursor peak of this peak list
     *
     * @return the precursor peak of this peak list
     */
    Peak getPrecursor();

    /**
     * Set the precursor peak for this peak list
     *
     * @param precursor the precursor peak
     */
    void setPrecursor(Peak precursor);

    /**
     * Returns a cursor over the peaks in this peak list in proper sequence.
     *
     * @return a cursor over the peaks in this peak list in proper sequence
     */
    PeakCursor<A> cursor();

    /**
     * Calculates the vector length of this peak list
     *
     * @return the vector length of this peak list
     */
    double calcVectorLength();

    /**
     * Apply the peakProcessor to the peaks in this PeakList
     *
     * @param peakProcessor the peak processor
     */
    void apply(PeakProcessor<A, A> peakProcessor);

    /**
     * Apply the peakProcessorChain to the peaks in this PeakList
     *
     * @param peakProcessorChain the peak processor list
     */
    void apply(PeakProcessorChain<A> peakProcessorChain);
}
