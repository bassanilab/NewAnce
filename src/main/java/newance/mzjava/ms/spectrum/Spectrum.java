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
package newance.mzjava.ms.spectrum;

import com.google.common.base.Function;
import com.google.common.base.Optional;
import com.google.common.base.Preconditions;
import newance.mzjava.ms.peaklist.*;

import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.UUID;

/**
 * @author Oliver Horlacher
 * @version 1.0
 */
public abstract class Spectrum<A extends PeakAnnotation> implements PeakList<A> {

    private int msLevel;
    private final PeakList<A> peakList;

    /**
     * Copy constructor
     *
     * @param src the Spectrum to copy
     * @param peakProcessor the processor to process the peaks
     */
    protected Spectrum(Spectrum<A> src, PeakProcessor<A, A> peakProcessor) {

        peakList = src.peakList.copy(peakProcessor);
        msLevel = src.msLevel;
    }

    /**
     * Copy constructor
     *
     * @param src the Spectrum to copy
     * @param peakProcessorChain the processor chain to process the peaks
     */
    protected Spectrum(Spectrum<A> src, PeakProcessorChain<A> peakProcessorChain) {

        peakList = src.peakList.copy(peakProcessorChain);
        msLevel = src.msLevel;
    }

    /**
     * Construct a Spectrum that has an initial capacity of <code>initialCapacity</code> and a precision
     * of <code>precision</code>
     *
     * @param initialCapacity the initial capacity
     * @param precision the precision
     */
    public Spectrum(int initialCapacity, Precision precision) {

        switch (precision) {
            case DOUBLE:

                peakList = new DoublePeakList<>(initialCapacity);
                break;
            case FLOAT:

                peakList = new FloatPeakList<>(initialCapacity);
                break;
            case DOUBLE_FLOAT:

                peakList = new DoubleFloatPeakList<>(initialCapacity);
                break;
            case DOUBLE_CONSTANT:

                peakList = new DoubleConstantPeakList<>(1.0, initialCapacity);
                break;
            case FLOAT_CONSTANT:

                peakList = new FloatConstantPeakList<>(1.0, initialCapacity);
                break;
            default:

                throw new IllegalArgumentException("Cannot create a peak list for precision " + precision);
        }
    }

    /**
     * Construct a Spectrum that has a constant intensity. The precision needs to be either
     * PeakList.Precision.DOUBLE_CONSTANT or PeakList.Precision.FLOAT_CONSTANT. The peak lists initial
     * capacity is set to <code>initialCapacity</code>
     *
     * @param initialCapacity   the initial capacity of the peak list
     * @param constantIntensity the constant intensity if
     * @param precision         the precision that the peak list is stored at
     */
    public Spectrum(int initialCapacity, double constantIntensity, Precision precision) {

        switch (precision) {

            case DOUBLE_CONSTANT:
                peakList = new DoubleConstantPeakList<>(constantIntensity, initialCapacity);
                break;
            case FLOAT_CONSTANT:
                peakList = new FloatConstantPeakList<>(constantIntensity, initialCapacity);
                break;
            default:
                throw new IllegalArgumentException("A constant intensity can only be set for DOUBLE_CONSTANT or FLOAT_CONSTANT precision.  Use new AbstractSpectrum(int initialCapacity, M spectrumMetaData, Precision precision)");
        }
    }

    @Override
    public UUID getId() {

        return peakList.getId();
    }

    @Override
    public void setId(UUID id) {

        peakList.setId(id);
    }

    @Override
    public void addPeaks(PeakList<A> peakList) {

        this.peakList.addPeaks(peakList);
    }

    @Override
    public <T extends PeakAnnotation> void addPeaksNoAnnotations(PeakList<T> peakList) {

        this.peakList.addPeaksNoAnnotations(peakList);
    }

    @Override
    public <T extends PeakAnnotation> void addPeaks(PeakList<T> peakList, Function<List<T>, List<A>> annotationConverter) {

        this.peakList.addPeaks(peakList, annotationConverter);
    }

    @Override
    public void addSorted(double[] mzs, double[] intensities) {

        peakList.addSorted(mzs, intensities);
    }

    @Override
    public void addSorted(double[] mzs, double[] intensities, int length) {

        peakList.addSorted(mzs, intensities, length);
    }

    @Override
    public int add(double mz, double intensity) {

        return peakList.add(mz, intensity);
    }

    @Override
    public int add(double mz, double intensity, A annotation) {

        return peakList.add(mz, intensity, annotation);
    }

    @Override
    public int add(double mz, double intensity, Collection<? extends A> annotations) {

        return peakList.add(mz, intensity, annotations);
    }

    @Override
    public int size() {

        return peakList.size();
    }

    @Override
    public boolean isEmpty() {

        return peakList.isEmpty();
    }

    @Override
    public double getMz(int index) {

        return peakList.getMz(index);
    }

    @Override
    public double getIntensity(int index) {

        return peakList.getIntensity(index);
    }

    @Override
    public void clear() {

        peakList.clear();
    }

    @Override
    public double getTotalIonCurrent() {

        return peakList.getTotalIonCurrent();
    }

    @Override
    public int getClosestIndex(double mz) {

        return peakList.getClosestIndex(mz);
    }

    @Override
    public int indexOf(double mzKey) {

        return peakList.indexOf(mzKey);
    }

    @Override
    public int indexOf(int fromIndex, int toIndex, double mzKey) {

        return peakList.indexOf(fromIndex, toIndex, mzKey);
    }

    @Override
    public int getMostIntenseIndex(double minMz, double maxMz) {

        return peakList.getMostIntenseIndex(minMz, maxMz);
    }

    @Override
    public int getMostIntenseIndex() {

        return peakList.getMostIntenseIndex();
    }

    @Override
    public double getBasePeakMz() {

        return peakList.getBasePeakMz();
    }

    @Override
    public double getBasePeakIntensity() {

        return peakList.getBasePeakIntensity();
    }

    @Override
    public double[] getMzs(double[] dest) {

        return peakList.getMzs(dest);
    }

    @Override
    public double[] getMzs(double[] dest, int destPos) {

        return peakList.getMzs(dest, destPos);
    }

    @Override
    public double[] getMzs(int srcPos, double[] dest, int destPos, int length) {

        return peakList.getMzs(srcPos, dest, destPos, length);
    }

    @Override
    public double[] getIntensities(double[] dest) {

        return peakList.getIntensities(dest);
    }

    @Override
    public double[] getIntensities(double[] dest, int destPos) {

        return peakList.getIntensities(dest, destPos);
    }

    @Override
    public double[] getIntensities(int srcPos, double[] dest, int destPos, int length) {

        return peakList.getIntensities(srcPos, dest, destPos, length);
    }

    @Override
    public Precision getPrecision() {

        return peakList.getPrecision();
    }

    @Override
    public void setIntensityAt(double intensity, int index) {

        peakList.setIntensityAt(intensity, index);
    }

    @Override
    public void trimToSize() {

        peakList.trimToSize();
    }

    @Override
    public void ensureCapacity(int minCapacity) {

        peakList.ensureCapacity(minCapacity);
    }

    @Override
    public void addAnnotation(int index, A annotation) {

        Preconditions.checkElementIndex(index, size());

        peakList.addAnnotation(index, annotation);
    }

    @Override
    public void addAnnotations(int index, Collection<A> annotations) {

        Preconditions.checkElementIndex(index, size());

        peakList.addAnnotations(index, annotations);
    }

    @Override
    public boolean removeAnnotation(A annotation, int index) {

        return peakList.removeAnnotation(annotation, index);
    }

    @Override
    public void clearAnnotationsAt(int index) {

        peakList.clearAnnotationsAt(index);
    }

    @Override
    public void clearAnnotations() {

        peakList.clearAnnotations();
    }

    @Override
    public boolean hasAnnotations() {

        return peakList.hasAnnotations();
    }

    @Override
    public int[] getAnnotationIndexes() {

        return peakList.getAnnotationIndexes();
    }

    @Override
    public boolean hasAnnotationsAt(int index) {

        return peakList.hasAnnotationsAt(index);
    }

    @Override
    public List<A> getAnnotations(int index) {

        return peakList.getAnnotations(index);
    }

    @Override
    public Optional<A> getFirstAnnotation(int index) {

        return peakList.getFirstAnnotation(index);
    }

    @Override
    public void sortAnnotations(Comparator<A> comparator) {

        peakList.sortAnnotations(comparator);
    }

    @Override
    public Peak getPrecursor() {

        return peakList.getPrecursor();
    }

    @Override
    public void setPrecursor(Peak precursor) {

        peakList.setPrecursor(precursor);
    }

    /**
     * Returns the ms level that this peak list was measured at
     *
     * @return the ms level that this peak list was measured at
     */
    public int getMsLevel() {

        return msLevel;
    }

    /**
     * Set the ms level that this peak list was measured at
     *
     * @param msLevel the ms level
     */
    public void setMsLevel(int msLevel) {

        this.msLevel = msLevel;
    }

    @Override
    public PeakCursor<A> cursor() {

        return peakList.cursor();
    }

    @Override
    public double calcVectorLength() {

        return peakList.calcVectorLength();
    }

    @Override
    public void apply(PeakProcessor<A, A> peakProcessor) {

        peakList.apply(peakProcessor);
    }

    @Override
    public void apply(PeakProcessorChain<A> peakProcessorChain) {

        peakList.apply(peakProcessorChain);
    }

    @Override
    public boolean equals(Object o) {

        if (this == o) return true;
        if (!(o instanceof Spectrum)) return false;

        Spectrum spectrum = (Spectrum) o;

        return msLevel == spectrum.msLevel &&
                peakList.equals(spectrum.peakList);

    }

    @Override
    public int hashCode() {
        int result = msLevel;
        result = 31 * result + peakList.hashCode();
        return result;
    }
}