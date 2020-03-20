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


import newance.mzjava.ms.peaklist.Peak;
import newance.mzjava.ms.peaklist.PeakAnnotation;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * A peak that has a list of peak annotations.
 *
 * @author Oliver Horlacher
 * @version 1.0
 */
public class AnnotatedPeak<A extends PeakAnnotation> extends Peak {

    private final List<A> annotations;

    /**
     * Default constructor. Default values are:
     * <pre>
     *     mz = 0.0
     *     intensity = 0.0
     *     charge = 0
     *     annotations = empty list
     * </pre>
     */
    public AnnotatedPeak() {

        super();

        annotations = new ArrayList<A>(1);
    }

    /**
     * Constructs a AnnotatedPeak.
     *
     * @param mz the peak m/z
     * @param intensity the peak
     * @param charge the peak charge
     * @param annotations the annotations
     */
    public AnnotatedPeak(double mz, double intensity, int charge, A... annotations) {

        super(mz, intensity, charge);
        this.annotations = new ArrayList<A>(annotations.length);
        Collections.addAll(this.annotations, annotations);
    }

    /**
     * Copy constructor.
     *
     * @param src the AnnotatedPeak to copy
     */
    public AnnotatedPeak(AnnotatedPeak<A> src) {

        super(src);

        annotations = new ArrayList<A>();
        for(A annotation : src.annotations) {

            //noinspection unchecked
            annotations.add((A)annotation.copy());
        }
    }

    /**
     * Add <code>annotation</code> to this AnnotatedPeak
     *
     * @param annotation the annotation to add
     * @return true if the annotation was added, false otherwise
     */
    public boolean add(A annotation) {

        return annotations.add(annotation);
    }

    /**
     * Remove <code>annotation</code> from this AnnotatedPeak
     *
     * @param annotation the annotation to remove
     * @return true if the annotation was removed, false otherwise
     */
    public boolean remove(A annotation) {

        return annotations.remove(annotation);
    }

    /**
     * Removes the annotation at the specified position in this peak.
     * Returns the annotation that was removed.
     *
     * @param index the index of the annotation to be removed
     * @return the annotation previously at the specified position     *
     * @throws IndexOutOfBoundsException if the index is out of range
     *         (<tt>index &lt; 0 || index &gt;= getAnnotationCount()</tt>)
     */
    public A remove(int index) {

        return annotations.remove(index);
    }

    /**
     * Returns the annotation at the specified position in this peak.
     *
     * @param index index of the annotation to return
     * @return the annotation at the specified position in this peak
     * @throws IndexOutOfBoundsException if the index is out of range
     *         (<tt>index &lt; 0 || index &gt;= getAnnotationCount()</tt>)
     */
    public A getAnnotation(int index) {

        return annotations.get(index);
    }

    /**
     * Returns the number of annotations in this peak.  If this peak contains
     * more than <tt>Integer.MAX_VALUE</tt> annotation, returns
     * <tt>Integer.MAX_VALUE</tt>.
     *
     * @return the number of annotations in this list
     */
    public int getAnnotationCount() {

        return annotations.size();
    }

    /**
     * Returns true this peak has annotations, false otherwise.
     *
     * @return true this peak has annotations, false otherwise
     */
    public boolean hasAnnotations() {

        return !annotations.isEmpty();
    }

    /**
     * Returns an unmodifiable list containing this peaks annotations
     *
     * @return an unmodifiable list containing this peaks annotations
     */
    public List<A> getAnnotations() {

        return Collections.unmodifiableList(annotations);
    }

    /**
     * Use the <code>comparator</code> to sort this peaks annotations
     *
     * @param comparator the comparator that is used to sort the annotations
     */
    public void sortAnnotations(Comparator<A> comparator) {

        Collections.sort(annotations, comparator);
    }

    @Override
    public AnnotatedPeak<A> copy() {

        return new AnnotatedPeak<A>(this);
    }

    @Override
    public String toString() {

        return "AnnotatedPeak{mz=" + getMz() + ", charge=" + getCharge() + ", intensity=" + getIntensity() + ", annotations=" + annotations + "}";
    }

    @Override
    public boolean equals(Object o) {

        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        if (!super.equals(o)) return false;

        AnnotatedPeak that = (AnnotatedPeak) o;

        return annotations.equals(that.annotations);
    }

    @Override
    public int hashCode() {

        int result = super.hashCode();
        result = 31 * result + annotations.hashCode();
        return result;
    }
}
