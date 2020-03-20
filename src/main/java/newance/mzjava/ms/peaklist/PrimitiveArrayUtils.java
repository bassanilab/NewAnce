/**
 * Copyright (c) 2012, SIB. All rights reserved.
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

import static com.google.common.base.Preconditions.checkArgument;
import static com.google.common.base.Preconditions.checkNotNull;

/**
 * Provides static generic methods on primitive arrays.
 *
 * @author nikitin
 */
public final class PrimitiveArrayUtils {

    private PrimitiveArrayUtils() {

    }

    /**
     * Grow the array arr by loadFactor
     *
     * @param arr        the array to grow
     * @param loadFactor the amount to grow the array
     * @return the new array
     */
    public static double[] grow(double[] arr, int loadFactor) {

        checkNotNull(arr);
        checkArgument(loadFactor > 0);

        double[] tmp = new double[arr.length + loadFactor];

        System.arraycopy(arr, 0, tmp, 0, arr.length);

        return tmp;
    }

    /**
     * Grow the array arr by loadFactor
     *
     * @param arr        the array to grow
     * @param loadFactor the amount to grow the array
     * @return the new array
     */
    public static float[] grow(float[] arr, int loadFactor) {

        checkNotNull(arr);
        checkArgument(loadFactor > 0);

        float[] tmp = new float[arr.length + loadFactor];

        System.arraycopy(arr, 0, tmp, 0, arr.length);

        return tmp;
    }

    /**
     * Insert the value into the array at index.
     *
     * @param value the value to insert
     * @param index the index in which the value is to be inserted
     * @param arr   the array to insert the value
     * @param size  the size of the peak list
     */
    public static void insert(double value, int index, double[] arr, int size) {

        System.arraycopy(arr, index, arr, index + 1, size - index);
        arr[index] = value;
    }

    /**
     * Insert the value into the array at index.
     *
     * @param value the value to insert
     * @param index the index in which the value is to be inserted
     * @param arr   the array to insert the value
     * @param size  the size of the peak list
     */
    public static void insert(float value, int index, float[] arr, int size) {

        System.arraycopy(arr, index, arr, index + 1, size - index);
        arr[index] = value;
    }

    /**
     * Wraps System.arrayCopy and creates a new array if required
     *
     * @param src     the src array
     * @param srcPos  the src position
     * @param dest    the dest array
     * @param destPos the dest position
     * @param length  the number of elements to copy
     * @return the dest array or the newly created array
     */
    public static double[] copyArray(double[] src, int srcPos, double[] dest, int destPos, int length) {

        if (dest == null) dest = new double[destPos + length];

        System.arraycopy(src, srcPos, dest, destPos, length);

        return dest;
    }

    /**
     * Copies an array from the specified source array, beginning at the
     * specified position, to the specified position of the destination array.
     * <p/>
     * Not using System.arrayCopy because the src and dest array are not the same type.
     *
     * @param src     the src array
     * @param srcPos  the src position
     * @param dest    the dest array
     * @param destPos the dest position
     * @param length  the number of elements to copy
     * @return the dest array or the newly created array
     */
    public static double[] copyArray(float[] src, int srcPos, double[] dest, int destPos, int length) {

        if (dest == null) dest = new double[destPos + length];

        for (int i = 0; i < length; i++) {

            dest[destPos + i] = src[srcPos + i];
        }

        return dest;
    }

    /**
     * Trim the arr to the size of this DoublePeakList
     *
     * @param arr the array to trim
     * @return the array
     */
    public static double[] trim(double[] arr, int size) {

        double[] tmp = new double[size];

        System.arraycopy(arr, 0, tmp, 0, size);

        return tmp;
    }

    /**
     * Trim the arr to the size of this DoublePeakList
     *
     * @param arr the array to trim
     * @return the array
     */
    public static float[] trim(float[] arr, int size) {

        float[] tmp = new float[size];

        System.arraycopy(arr, 0, tmp, 0, size);

        return tmp;
    }

    /**
     * Array hash code that takes into account the peak list size
     *
     * @param a the array
     * @return the hash code
     */
    public static int arrayHashCode(double a[], int size) {

        if (a == null)
            return 0;

        int result = 1;
        for (int i = 0; i < size; i++) {

            double element = a[i];

            long bits = Double.doubleToLongBits(element);
            result = 31 * result + (int) (bits ^ (bits >>> 32));
        }
        return result;
    }

    /**
     * Array hash code that takes into account the peak list size
     *
     * @param a the array
     * @return the hash code
     */
    public static int arrayHashCode(float a[], int size) {

        if (a == null)
            return 0;

        int result = 1;
        for (int i = 0; i < size; i++) {

            result = 31 * result + Float.floatToIntBits(a[i]);
        }

        return result;
    }

    /**
     * Array equals that takes into account the size
     *
     * @param a  the first array
     * @param a2 the second array
     * @return true if the arrays are equals, false otherwise
     */
    public static boolean arrayEquals(double[] a, double[] a2, int size) {

        if (a == a2)
            return true;
        if (a == null || a2 == null)
            return false;

        for (int i = 0; i < size; i++)
            if (Double.doubleToLongBits(a[i]) != Double.doubleToLongBits(a2[i]))
                return false;

        return true;
    }

    /**
     * Array equals that takes into account the size
     *
     * @param a  the first array
     * @param a2 the second array
     * @return true if the arrays are equals, false otherwise
     */
    public static boolean arrayEquals(float[] a, float[] a2, int size) {

        if (a == a2)
            return true;
        if (a == null || a2 == null)
            return false;

        for (int i = 0; i < size; i++)
            if (Float.floatToIntBits(a[i]) != Float.floatToIntBits(a2[i]))
                return false;

        return true;
    }

    /**
     * Load doubles from src to dest arrays.
     *
     * @param src  the source array to load from.
     * @param dest the destination array to load to.
     * @return the dest array.
     */
    public static double[] loadDoubles(double[] src, double[] dest) {

        return loadDoubles(src, dest, 0);
    }

    /**
     * Load doubles from src to dest arrays.
     *
     * @param src     the source array to load from.
     * @param dest    the destination array to load to.
     * @param destPos starting position in the destination data.
     * @return the dest array.
     */
    public static double[] loadDoubles(double[] src, double[] dest, int destPos) {

        checkNotNull(src, "undefined src array to load into dest!");

        if (dest == null || dest.length < src.length)
            dest = new double[src.length];

        System.arraycopy(src, 0, dest, destPos, src.length);

        return dest;
    }

    /**
     * Convert double array to Double array.
     *
     * @param src  the source array to load from.
     * @param dest the destination array to load to.
     * @return the dest array.
     */
    public static Double[] toDoubles(double[] src, Double[] dest) {

        return toDoubles(src, dest, 0);
    }

    /**
     * Convert double array to Double array.
     *
     * @param src     the source array to load from.
     * @param dest    the destination array to load to.
     * @param destPos starting position in the destination data.
     * @return the dest array.
     */
    public static Double[] toDoubles(double[] src, Double[] dest,
                                     int destPos) {

        checkNotNull(src, "undefined src array to load into dest!");

        if (dest == null || dest.length < src.length)
            dest = new Double[src.length];

        for (int i = 0; i < src.length; i++) {

            dest[i + destPos] = src[i];
        }

        return dest;
    }

    /**
     * Get the index of the maximum value
     *
     * @param values the array to search for the max
     * @return an index or -1 if empty
     */
    public static int indexMax(double[] values) {

        Preconditions.checkNotNull(values);

        double valueMax = Double.MIN_VALUE;
        int indexMax = -1;

        for (int i = 0; i < values.length; i++) {

            if (values[i] > valueMax) {

                valueMax = values[i];
                indexMax = i;
            }
        }

        return indexMax;
    }
}

