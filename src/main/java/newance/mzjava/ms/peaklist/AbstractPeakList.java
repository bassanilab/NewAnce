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
import com.google.common.base.Preconditions;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static com.google.common.base.Preconditions.checkNotNull;
import static com.google.common.base.Preconditions.checkPositionIndex;

/**
 * @author Oliver Horlacher
 * @author nikitin
 * @version 1.0
 */
public abstract class AbstractPeakList<A extends PeakAnnotation> implements PeakList<A> {

    protected double totalIonCurrent = 0.0;

    protected int size = 0;
    protected int loadFactor = 100;

    final TIntObjectMap<List<A>> annotationMap = new TIntObjectHashMap<>();
    private Peak precursor;

    private UUID id;

    protected AbstractPeakList() {

        precursor = new Peak();
    }

    protected AbstractPeakList(AbstractPeakList<A> src) {

        Peak srcPrecursor = src.precursor;
        if (srcPrecursor != null) {

            precursor = srcPrecursor.copy();
        }
        this.id = src.id;
    }

    /**
     * Create a peak list from a string
     *
     * @param str the string
     * @return the new peakList
     */
    public static <A extends PeakAnnotation> PeakList<A> valueOf(String str, Precision precision) {

        PeakList<A> peakList = PeakListFactory.newPeakList(precision, 100);
        valueOf(str, peakList);
        return peakList;
    }

    /**
     * Create a peak list from a string
     *
     * @param str the string
     */
    protected static void valueOf(String str, PeakList peakList) {

        Pattern pat =
                Pattern.compile("([+-]?\\d+\\.?\\d*(?:[eE][-+]?\\d+)?)\\s+(?:Da\\s*)?(\\([+-]?\\d+\\.?\\d*(?:[eE][-+]?\\d+)?\\))?");

        Matcher matcher = pat.matcher(str);

        while (matcher.find()) {
            double mz = Double.parseDouble(matcher.group(1));
            double intensity = 0;
            if (matcher.group(2) != null) {
                intensity =
                        Double.parseDouble(matcher.group(2)
                                .substring(1, matcher.group(2).length() - 1));
            }

            peakList.add(mz, intensity);
        }
    }

    @Override
    public UUID getId() {

        if (id == null) id = UUID.randomUUID();

        return id;
    }

    @Override
    public void setId(UUID id) {

        checkNotNull(id);

        this.id = id;
    }

    /**
     * Checks if the given index is in range.  If not, throws an appropriate
     * runtime exception.  This method does *not* check if the index is
     * negative: It is always used immediately prior to an array access,
     * which throws an ArrayIndexOutOfBoundsException if index is negative.
     *
     * @param index the index
     */
    protected void rangeCheck(int index) {

        if (index >= size) {

            throw new IndexOutOfBoundsException("Index: " + index + ", Size: " + size);
        }
    }

    @Override
    public int add(double mz, double intensity) {

        grow();

        int index;
        if (size > 0 && mz < getMz(size - 1)) {

            index = doInsert(mz, intensity);
        } else {

            index = doAppend(mz, intensity);
        }

        return index;
    }

    @Override
    public int add(double mz, double intensity, A annotation) {

        checkNotNull(annotation);

        int index = add(mz, intensity);
        addAnnotation(index, annotation);

        return index;
    }

    @Override
    public int add(double mz, double intensity, Collection<? extends A> annotations) {

        checkNotNull(annotations);

        int index = add(mz, intensity);
        for (A annotation : annotations) {

            addAnnotation(index, annotation);
        }

        return index;
    }

    /**
     * Shift all annotations that have an index that is >= insertionIndex by 1
     *
     * @param insertionIndex the index where the peak insert occurred
     */
    protected void shiftAnnotations(int insertionIndex) {

        for (int index : annotationMap.keys()) {

            if (index >= insertionIndex) {

                annotationMap.put(index + 1, annotationMap.remove(index));
            }
        }
    }

    /**
     * Insert the peak so that the m/z are still sorted
     *
     * @param mz        the peak m/z
     * @param intensity the peak intensity
     * @return the index at which the peak was added
     */
    protected abstract int doInsert(double mz, double intensity);

    /**
     * Append the peak to the end of the peak list
     *
     * @param mz        the peak m/z
     * @param intensity the peak intensity
     */
    protected abstract int doAppend(double mz, double intensity);

    /**
     * Grow the backing if the array is full.
     */
    protected void grow() {

        ensureCapacity(((size / loadFactor) + 1) * loadFactor);
    }

    @Override
    public void addSorted(double[] mzs, double[] intensities) {

        checkNotNull(mzs);
        checkNotNull(intensities);

        if (mzs.length != intensities.length)
            throw new IllegalStateException("The mz and intensity arrays need to have the same length, the mz array has a length of " + mzs.length + " and the intensity array a length of " + intensities.length);

        addSorted(mzs, intensities, mzs.length);
    }

    @Override
    public void addSorted(double[] mzs, double[] intensities, int length) {

        if (length == 0) return;

        checkNotNull(mzs, "mzs cannot be null");
        checkNotNull(intensities, "intensities cannot be null");
        checkPositionIndex(length, mzs.length);
        checkPositionIndex(length, intensities.length);

        PeakCursor<A> arrayCursor = new ArrayCursor<>(mzs, intensities, length);
        PeakCursor<A> listCursor = cursor();

        PeakListMerger<A> merger = new PeakListMerger<>();
        merger.setSink(newMergeSink());
        merger.merge(listCursor, arrayCursor);
    }

    @Override
    public void addPeaks(PeakList<A> peakList) {

        checkNotNull(peakList);
        if (peakList.isEmpty()) return;

        PeakListMerger<A> merger = new PeakListMerger<>();
        merger.setSink(newMergeSink());
        merger.merge(peakList.cursor(), cursor());
    }


    @Override
    public <T extends PeakAnnotation> void addPeaksNoAnnotations(PeakList<T> peakList) {

        addPeaks(peakList, new RemoveAnnotationFunction<T, A>());
    }

    @Override
    public <T extends PeakAnnotation> void addPeaks(PeakList<T> peakList, Function<List<T>, List<A>> annotationConverter) {

        checkNotNull(peakList);

        PeakListMerger<A> merger = new PeakListMerger<>();
        merger.setSink(newMergeSink());
        merger.merge(new WrapperPeakCursor<>(peakList.cursor(), annotationConverter), cursor());
    }

    abstract PeakSink<A> newMergeSink();

    @Override
    public int getClosestIndex(double mz) {

        if (size == 0) {

            return -1;
        }

        int index = indexOf(0, size, mz);

        if (index < 0) index = -1 * (index + 1);

        if (index == size) {

            return index - 1;
        } else if (index == 0) {

            return index;
        }

        double ds = mz - getMz(index - 1);
        double dl = getMz(index) - mz;

        return ds <= dl ? index - 1 : index;
    }

    @Override
    public int indexOf(double mzKey) {

        return indexOf(0, size, mzKey);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int size() {

        return size;
    }

    @Override
    public boolean isEmpty() {

        return size == 0;
    }

    /**
     * Clears the peak list, the size is set to 0 and the backing arrays are set to an clear array
     */
    @Override
    public void clear() {

        size = 0;
        totalIonCurrent = 0;
        trimToSize();

        annotationMap.clear();
    }

    /**
     * Returns the current load factor, which is the number of elements that are added to the mz and intensity arrays
     * when they are full
     *
     * @return the current load factor
     */
    public int getLoadFactor() {

        return loadFactor;
    }

    /**
     * Set the number of elements that are added to the mz and intensity arrays when the size of this instance exceeds
     * the capacity of these arrays
     *
     * @param loadFactor the new load factor
     */
    public void setLoadFactor(int loadFactor) {

        this.loadFactor = loadFactor;
    }

    /**
     * Returns the sum of all intensities
     *
     * @return the sum of all intensities
     */
    @Override
    public double getTotalIonCurrent() {

        return totalIonCurrent;
    }

    @Override
    public boolean hasAnnotations() {

        return !annotationMap.isEmpty();
    }

    @Override
    public int[] getAnnotationIndexes() {

        int[] indexes = new int[annotationMap.size()];

        int c = 0;
        for (int index : annotationMap.keys()) {

            indexes[c++] = index;
        }
        Arrays.sort(indexes);

        return indexes;
    }

    @Override
    public boolean hasAnnotationsAt(int index) {

        return annotationMap.containsKey(index);
    }

    @Override
    public List<A> getAnnotations(int index) {

        if (annotationMap.containsKey(index)) {

            return Collections.unmodifiableList(annotationMap.get(index));
        } else {

            return Collections.emptyList();
        }
    }

    @Override
    public Optional<A> getFirstAnnotation(int index) {

        List<A> annotationList = annotationMap.get(index);
        if (annotationList != null) {

            return Optional.of(annotationList.get(0));
        } else {

            return Optional.absent();
        }
    }

    @Override
    public void addAnnotation(int index, A annotation) {

        checkNotNull(annotation);
        Preconditions.checkElementIndex(index, size());

        List<A> annotations = annotationMap.get(index);
        if (annotations == null) {

            annotations = new ArrayList<>(1);
            annotationMap.put(index, annotations);
        }

        annotations.add(annotation);
    }

    @Override
    public void addAnnotations(int index, Collection<A> annotations) {

        Preconditions.checkElementIndex(index, size());

        List<A> currAnnotations = annotationMap.get(index);
        if (currAnnotations == null) {

            currAnnotations = new ArrayList<>(annotations.size());
            annotationMap.put(index, currAnnotations);
        }

        currAnnotations.addAll(annotations);
    }

    @Override
    public boolean removeAnnotation(A annotation, int index) {

        List<A> annotations = annotationMap.get(index);
        if (annotations != null) {

            boolean removed = annotations.remove(annotation);

            if (annotations.isEmpty()) annotationMap.remove(index);
            return removed;
        } else {

            return false;
        }
    }

    @Override
    public void clearAnnotationsAt(int index) {

        annotationMap.remove(index);
    }

    @Override
    public void clearAnnotations() {

        annotationMap.clear();
    }

    @Override
    public void sortAnnotations(Comparator<A> comparator) {

        for (List<A> annotations : annotationMap.valueCollection()) {

            Collections.sort(annotations, comparator);
        }
    }

    @Override
    public Peak getPrecursor() {

        return precursor;
    }

    @Override
    public void setPrecursor(Peak precursor) {

        checkNotNull(precursor);

        this.precursor = precursor;
    }

    @Override
    public boolean equals(Object o) {

        if (this == o) return true;
        if (!(o instanceof AbstractPeakList)) return false;

        AbstractPeakList that = (AbstractPeakList) o;

        return size == that.size &&
                Double.compare(that.totalIonCurrent, totalIonCurrent) == 0 &&
                annotationMap.equals(that.annotationMap) &&
                !(precursor != null ? !precursor.equals(that.precursor) : that.precursor != null) &&
                !(id != null ? !id.equals(that.id) : that.id != null);
    }

    @Override
    public int hashCode() {

        int result;
        long temp;
        temp = totalIonCurrent != +0.0d ? Double.doubleToLongBits(totalIonCurrent) : 0L;
        result = (int) (temp ^ (temp >>> 32));
        result = 31 * result + size;
        result = 31 * result + annotationMap.hashCode();
        result = 31 * result + (precursor != null ? precursor.hashCode() : 0);
        result = 31 * result + (id != null ? id.hashCode() : 0);
        return result;
    }

    @Override
    public void apply(PeakProcessor<A, A> peakProcessor) {

        checkNotNull(peakProcessor);

        apply(this, this, peakProcessor);
    }

    protected void apply(PeakList<A> src, PeakList<A> dest, PeakProcessor<A, A> peakProcessor) {

        checkNotNull(peakProcessor);

        PeakCollectorSink<A> peakSink = new PeakCollectorSink<>();
        peakSink.start(src.size());
        peakSink.setPeakList(dest);

        peakProcessor.setSink(peakSink);

        int peakCount = src.size();
        peakProcessor.start(peakCount);
        for (int i = 0; i < peakCount; i++) {

            peakProcessor.processPeak(src.getMz(i), src.getIntensity(i), src.getAnnotations(i));
        }
        peakProcessor.end();
        dest.trimToSize();
    }

    @Override
    public void apply(PeakProcessorChain<A> peakProcessorChain) {

        checkNotNull(peakProcessorChain);

        apply(this, this, peakProcessorChain);
    }

    protected void apply(PeakList<A> org, PeakList<A> dest, PeakProcessorChain<A> peakProcessorChain) {

        checkNotNull(peakProcessorChain);
        if (peakProcessorChain.isEmpty()) return;

        PeakCollectorSink<A> peakSink = new PeakCollectorSink<>();
        peakSink.setPeakList(dest);

        peakProcessorChain.process(org, peakSink);
        trimToSize();
    }

    @Override
    public double getBasePeakMz() {

        return getMz(getMostIntenseIndex());
    }

    @Override
    public double getBasePeakIntensity() {

        return getIntensity(getMostIntenseIndex());
    }

    /**
     * Finds the index that has of the peak that has an mz which is equal to or larger than mz
     *
     * @param mz           the mz to test
     * @param defaultIndex the index to return if there is no peak in this peak lis that has a mz that is >= mz
     * @return the index that has an mz that is equal to or larger than mz or defaultIndex if no such index exists
     */
    protected int indexEqualOrLarger(double mz, int defaultIndex, double[] mzList) {

        int index = getClosestIndex(mz);

        if (index == -1) {

            return defaultIndex;
        } else if (mzList[index] < mz) {

            return index + 1;
        } else {

            return index;
        }
    }

    /**
     * Finds the index that has of the peak that has an mz which is equal to or larger than mz
     *
     * @param mz           the mz to test
     * @param defaultIndex the index to return if there is no peak in this peak lis that has a mz that is >= mz
     * @return the index that has an mz that is equal to or larger than mz or defaultIndex if no such index exists
     */
    protected int indexEqualOrLarger(float mz, int defaultIndex, float[] mzList) {

        int index = getClosestIndex(mz);

        if (index == -1) {

            return defaultIndex;
        } else if (mzList[index] < mz) {

            return index + 1;
        } else {

            return index;
        }
    }

    @Override
    public String toString() {

        return getClass().getSimpleName() + "{" +
                "precursor=" + precursor +
                ", size=" + size +
                '}';
    }

    public PeakCursor<A> cursor() {

        return new PeakListCursor<>(this);
    }

    protected abstract static class AbstractMergePeakSink<A extends PeakAnnotation, L extends AbstractPeakList<A>> implements PeakSink<A> {

        protected int tmpSize = 0;
        protected double tmpTotalIonCurrent = 0;
        protected final Map<Integer, List<A>> tmpAnnotationMap = new HashMap<>();

        protected final L peakList;

        protected AbstractMergePeakSink(L peakList) {

            this.peakList = peakList;
        }

        @Override
        public void processPeak(double mz, double intensity, List<A> annotations) {

            if (tmpSize > 0 && mz == getMz(tmpSize - 1)) {

                addIntensity(intensity, tmpSize - 1);

                if (!annotations.isEmpty()) {

                    List<A> currentAnnotations = tmpAnnotationMap.get(tmpSize - 1);
                    if (currentAnnotations == null) {

                        tmpAnnotationMap.put(tmpSize - 1, new ArrayList<>(annotations));
                    } else {

                        currentAnnotations.addAll(annotations);
                    }
                }
            } else {

                if (!(tmpSize == 0 || mz >= getMz(tmpSize - 1))) {

                    throw new UnsortedPeakListException(tmpSize);
                }

                addToArray(mz, intensity, tmpSize);

                if (!annotations.isEmpty()) {

                    tmpAnnotationMap.put(tmpSize, new ArrayList<>(annotations));
                }

                tmpSize++;
            }
            tmpTotalIonCurrent += intensity;
        }

        protected abstract void addIntensity(double intensity, int index);

        @Override
        public void end() {

            peakList.size = tmpSize;
            peakList.totalIonCurrent = tmpTotalIonCurrent;
            setArrays();

            peakList.annotationMap.clear();
            peakList.annotationMap.putAll(tmpAnnotationMap);
        }

        protected abstract void setArrays();

        protected abstract double getMz(int index);

        protected abstract void addToArray(double mz, double intensity, int tmpSize);
    }

    /**
     * A wrapper for a peak cursor that transforms annotations using a Function.
     *
     * @param <T> the input PeakAnnotation type
     * @param <A> the output PeakAnnotation type
     * @see WrapperPeakCursor#currAnnotations()
     */
    static class WrapperPeakCursor<T extends PeakAnnotation, A extends PeakAnnotation> implements PeakCursor<A> {

        private final PeakCursor<T> cursor;
        private final Function<List<T>, List<A>> annotationFunction;

        WrapperPeakCursor(PeakCursor<T> cursor, Function<List<T>, List<A>> annotationFunction) {

            this.cursor = cursor;
            this.annotationFunction = annotationFunction;
        }

        @Override
        public int size() {

            return cursor.size();
        }

        @Override
        public boolean isEmpty() {

            return cursor.isEmpty();
        }

        @Override
        public void resetCursor() {

            cursor.resetCursor();
        }

        @Override
        public boolean next() {

            return cursor.next();
        }

        @Override
        public boolean previous() {

            return cursor.previous();
        }

        @Override
        public boolean next(double mz) {

            return cursor.next(mz);
        }

        @Override
        public double currIntensity() {

            return cursor.currIntensity();
        }

        /**
         * Returns an empty annotation list
         *
         * @return an empty annotation list
         */
        public List<A> currAnnotations() {

            return annotationFunction.apply(cursor.currAnnotations());
        }

        @Override
        public double currMz() {

            return cursor.currMz();
        }

        @Override
        public boolean canPeek(int n) {

            return cursor.canPeek(n);
        }

        @Override
        public double peekIntensity(int n) {

            return cursor.peekIntensity(n);
        }

        @Override
        public double peekMz(int n) {

            return cursor.peekMz(n);
        }

        @Override
        public double lastMz() {

            return cursor.lastMz();
        }

        @Override
        public double lastIntensity() {

            return cursor.lastIntensity();
        }

        @Override
        public double getMz(int index) {

            return cursor.getMz(index);
        }

        @Override
        public double getIntensity(int index) {

            return cursor.getIntensity(index);
        }

        @Override
        public void moveToClosest(double mz) {

            cursor.moveToClosest(mz);
        }

        @Override
        public void moveBefore(double mz) {

            cursor.moveBefore(mz);
        }

        @Override
        public void movePast(double mz) {

            cursor.movePast(mz);
        }

        @Override
        public int getClosestIndex(double mz) {

            return cursor.getClosestIndex(mz);
        }
    }

    static class RemoveAnnotationFunction<IN extends PeakAnnotation, OUT extends PeakAnnotation> implements Function<List<IN>, List<OUT>> {

        private final List<OUT> emptyList = Collections.emptyList();

        @Override
        public List<OUT> apply(List<IN> input) {

            return emptyList;
        }
    }
}
