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

import java.util.Arrays;

import static com.google.common.base.Preconditions.checkArgument;

/**
 * A peak in a mass spectrometer-typed spectrum must be an ion with a mass
 * value.
 *
 * @author nikitin
 * @author Oliver Horlacher
 * @version 1.0
 */
public class Peak implements Copyable<Peak>, Comparable<Peak> {

    private double mz = 0.0;
    private double intensity = 0.0;
    private int[] chargeList = new int[0];
    private double mass = 0.0;
    private Polarity polarity = Polarity.UNKNOWN;

    /**
     * Default constructor. Default values are:
     * <pre>
     *     mz = 0.0
     *     intensity = 0.0
     *     charge = 0
     * </pre>
     */
    public Peak() {

    }

    /**
     * Constructs a Peak.
     *
     * @param mz        the peak m/z
     * @param intensity the peak
     * @param charge    the peak charge
     */
    public Peak(double mz, double intensity, int... charge) {

        setValues(mz, intensity, charge);
    }

    /**
     * Factory method to construct a Peak that has an intensity of 0
     *
     * @param mz     the peaks m/z
     * @param charge the peaks charge
     * @return the new Peak
     */
    public static Peak noIntensity(double mz, int... charge) {

        return new Peak(mz, 0, charge);
    }

    /**
     * Copy constructor
     *
     * @param peak the peak to copy
     */
    public Peak(Peak peak) {

        mz = peak.mz;
        mass = peak.mass;
        intensity = peak.intensity;
        chargeList = peak.chargeList;
        polarity = peak.polarity;
        mass = peak.mass;
    }

    /**
     * Set this peaks values
     *
     * @param mz        the peak m/z
     * @param intensity the peak
     * @param charge    the peak charge
     */
    public void setValues(double mz, double intensity, int... charge) {

        checkArgument(!Double.isInfinite(mz), "m/z cannot be infinite");
        checkArgument(!Double.isNaN(mz), "m/z cannot be NaN");
        checkArgument(!Double.isInfinite(intensity), "intensity cannot be infinite");
        checkArgument(!Double.isNaN(intensity), "intensity cannot be NaN");

        setMzAndCharge(mz, charge);

        this.intensity = intensity;
    }

    /**
     * Return the m/z of the peak
     *
     * @return the peak mz
     */
    public double getMz() {

        return mz;
    }

    /**
     * Returns mass of this peak
     *
     * @return m in the m/z
     */
    public double getMass() {

        if (chargeList.length == 0)
            throw new IllegalStateException("The mass is undefined because the peak does not have any charge");

        return mass;
    }

    /**
     * Return the default chargeList of this peak
     *
     * @return the chargeList
     */
    public int getCharge() {

        return chargeList.length == 0 ? 0 : chargeList[0];
    }

    /**
     * Return the charge polarity.
     *
     * @return the charge polarity
     */
    public Polarity getPolarity() {

        return polarity;
    }

    /**
     * Return an array containing this peaks charges.
     *
     * @return an array containing this peaks charges
     */
    public int[] getChargeList() {

        int[] list = new int[chargeList.length];
        System.arraycopy(chargeList, 0, list, 0, chargeList.length);

        return list;
    }

    /**
     * Return the peak intensity
     *
     * @return the peak intensity
     */
    public double getIntensity() {

        return intensity;
    }

    /**
     * Set the m/z and chargeList of this peak
     *
     * @param mz     the new m/z
     * @param charge the new chargeList
     */
    public void setMzAndCharge(double mz, int... charge) {

        if (charge.length == 1 && charge[0] == 0) charge = new int[0];

        polarity = charge.length > 0 ? Polarity.getPolarity(charge[0]) : Polarity.UNKNOWN;

        for (int i = charge.length - 1; i >= 0; i--) {

            int z = charge[i];
            Preconditions.checkArgument(z != 0, "Charge cannot be 0. An unknown charge is specified by a empty array (new int[0])");
            charge[i] = Math.abs(z);
        }

        if (chargeList == null || chargeList.length != charge.length) {
            chargeList = new int[charge.length];
        }
        System.arraycopy(charge, 0, chargeList, 0, charge.length);

        if (chargeList.length >= 1) {

            mass = mz * chargeList[0];
        } else {

            mass = 0;
        }

        this.mz = mz;
    }

    /**
     * Set the intensity of this peak
     *
     * @param intensity the new intensity
     */
    public void setIntensity(double intensity) {

        this.intensity = intensity;
    }

    @Override
    public boolean equals(Object o) {

        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Peak peak = (Peak) o;

        return polarity == peak.polarity &&
                Double.compare(peak.intensity, intensity) == 0 &&
                Double.compare(peak.mz, mz) == 0 &&
                Arrays.equals(chargeList, peak.chargeList);
    }

    @Override
    public int hashCode() {

        int result;
        long temp;
        temp = mz != +0.0d ? Double.doubleToLongBits(mz) : 0L;
        result = (int) (temp ^ (temp >>> 32));
        temp = intensity != +0.0d ? Double.doubleToLongBits(intensity) : 0L;
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        result = 31 * result + Arrays.hashCode(chargeList);
        result = 31 * result + polarity.hashCode();
        return result;
    }

    @Override
    public Peak copy() {

        return new Peak(this);
    }

    @Override
    public int compareTo(Peak peak) {

        int cmp = Double.compare(getCharge(), peak.getCharge());
        if (cmp != 0) {
            return cmp;
        }
        cmp = Double.compare(mz, peak.mz);
        if (cmp != 0) {
            return cmp;
        }
        return Double.compare(intensity, peak.intensity);
    }

    @Override
    public String toString() {

        //noinspection StringBufferReplaceableByString
        return new StringBuilder("Peak{").append("mz=").append(mz).append(", charge=")
                .append(Arrays.toString(chargeList)).append(", polarity=").append(polarity).append(", intensity=").append(intensity).append('}').toString();
    }
}
