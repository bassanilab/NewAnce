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


import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TIntObjectHashMap;
import newance.mzjava.ms.peaklist.DelayedPeakProcessor;
import newance.mzjava.ms.peaklist.PeakAnnotation;

import java.util.Collections;
import java.util.List;

import static com.google.common.base.Preconditions.checkArgument;
import static com.google.common.base.Preconditions.checkNotNull;

/**
 * Class to select all peaks that are local maxima and replace their
 * mass/intensity values by the centroid values calculated within a
 * neighborhood. It is assumed that the mass/intensity signal is smooth, i.e.
 * that the signal is concave around the peak. The intensity is set to
 * the value at the apex of the centroid.
 *
 * @author mueller
 * @version 1.0
 */
public class CentroidFilter<A extends PeakAnnotation>  extends DelayedPeakProcessor<A, A> {

    public enum IntensityMode {HIGHEST, SUM, AVG}

    /**
     * Internal parameters used in the algorithm.
     * iMin, iMax : lower and upper bound of sliding window indices
     * iLowBound, iUpBound : start and end indices of centroid mz values
     */
    private int iMin, iMax, iLowBound, iUpBound;
    /**
     * Maximum distance between consecutive mz values in same peak group
     */
    private final double mzMaxDiff;

    private final IntensityMode intensityMode;

    /**
     *
     * @param mzMaxDiff Maximum distance between consecutive mz values in same peak group
     * @param intensityMode Mode for the calculation of centroid intensity
     *                      - IntensityMode.HIGHEST : set to highest intensity in peak
     *                      - IntensityMode.SUM : set to summed intensities in peak
     *                      - IntensityMode.AVG : set to average intensity in peak
     */
    public CentroidFilter(double mzMaxDiff, IntensityMode intensityMode) {

        checkArgument(mzMaxDiff > 0);
        checkNotNull(intensityMode);

        this.intensityMode = intensityMode;

        this.mzMaxDiff = mzMaxDiff;
    }

    public static double calcMzSamplingDist(TDoubleArrayList mzList) {

        double minDiff = 1000000.0;
        double diff;

        mzList.sort();
        // Calc minimal difference between 2 consecutive mz values
        for (int i = 1; i < mzList.size(); i++) {
            diff = mzList.get(i) - mzList.get(i-1);
            if (diff < minDiff) {
                minDiff = diff;
            }
        }

        return minDiff;
    }

    /**
     * Finds lower and upper index bound for peaks belonging to centroid at index.
     */
    private void calcBounds(TDoubleArrayList mzList, TDoubleArrayList intensityList, int index)
    {
        iLowBound = Math.max(index,0);
        // find lower bound for centroid, i.e. the point where the intensity becomes 0 or starts raising again.
        while (iLowBound > 0  && mzList.get(iLowBound) - mzList.get(iLowBound-1) < mzMaxDiff) {
            if (intensityList.get(iLowBound-1)<=0 || intensityList.get(iLowBound) < intensityList.get(iLowBound - 1)) {
                break;
            }
            iLowBound--;
        }
        // find upper bound for centroid, i.e. the point where the intensity becomes 0 or starts raising again.
        int size = mzList.size();
        iUpBound = index;
        while (iUpBound < size-1 && mzList.get(iUpBound+1) - mzList.get(iUpBound) < mzMaxDiff) {
            if (intensityList.get(iUpBound+1)<=0 || intensityList.get(iUpBound) < intensityList.get(iUpBound+1)) {
                break;
            }
            iUpBound++;
        }
    }

    /**
     *  Checks whether peak at index is local maximum within +/- mzMaxDiff
     */
    private boolean isMaximum(TDoubleArrayList mzList, TDoubleArrayList intensityList, int index) {
        final double mz = mzList.get(index);
        final double h = intensityList.get(index);

        // adjust left border of sliding window
        while (iMin<mzList.size() && mzList.get(iMin) <  mz-mzMaxDiff) iMin++;
        // adjust right border of sliding window
        while (iMax<mzList.size() && mzList.get(iMax) <= mz+mzMaxDiff) iMax++;
        iMax--;

        // pick N highest peak within window
        boolean isMax = true;
        for (int i = iMin; i <= iMax; i++) {
            if (intensityList.get(i) > h) {
                isMax = false;
                break;
            }
        }

        // store indices of peaks
        return isMax;
    }


    /**
     * Given the index-bounds of a centroid, the centroid mass and intensity is calculated. The intensity is set to
     * the value at the apex of the centroid.
     */
    protected double calcCentroidMz(TDoubleArrayList mzList, TDoubleArrayList intensityList, PeakIntensityInfo neighborhood) {

        double htot = 0.0;
        double m = 0.0;

        for (int j = iLowBound; j <= iUpBound; j++) {
            double h = intensityList.get(j);
            htot += intensityList.get(j);
            m += h * mzList.get(j);
        }

        if (htot > 0) {
            m = m / htot;
        } else {
            m = m /(iUpBound-iLowBound+1);
        }

        neighborhood.setPeakCount(iUpBound-iLowBound+1);
        neighborhood.setTotalIntensity(htot);

        return m;
    }

    @Override
    protected void processCached(TDoubleArrayList mzList, TDoubleArrayList intensityList, TIntObjectHashMap<List<A>> annotationMap) {

        iMin = iMax = 0;

        // collect informations on peaks reduced into centroid
        PeakIntensityInfo neighborhood = new PeakIntensityInfo();

        for (int i = 0; i < mzList.size(); i++) {
            if (isMaximum(mzList,intensityList,i) && intensityList.get(i)>0.0) {

                calcBounds(mzList,intensityList,i);
                double centroidMz = calcCentroidMz(mzList,intensityList, neighborhood);

                neighborhood.setHighestIntensity(intensityList.get(i));

                List<A> annotations;

                if (annotationMap.contains(i))
                    annotations = annotationMap.get(i);
                else
                    annotations = Collections.emptyList();

                sink.processPeak(centroidMz, getCentroidIntensity(neighborhood), annotations);
            }
        }
    }


    protected double getCentroidIntensity(PeakIntensityInfo peakInfo) {

        switch (intensityMode) {
            case HIGHEST:
                return peakInfo.getHighestIntensity();
            case SUM:
                return peakInfo.getTotalIntensity();
            case AVG:
                return peakInfo.getTotalIntensity() / peakInfo.getPeakCount();
            default:
                throw new IllegalStateException("Cannot calc centroid intensity for " + intensityMode);
        }
    }

    public static class PeakIntensityInfo {

        private double highestIntensity;
        private double totalIntensity;
        private int peakCount;

        public double getHighestIntensity() {
            return highestIntensity;
        }

        public void setHighestIntensity(double highestIntensity) {
            this.highestIntensity = highestIntensity;
        }

        public void setTotalIntensity(double totalIntensity) {
            this.totalIntensity = totalIntensity;
        }

        public void setPeakCount(int peakCount) {
            this.peakCount = peakCount;
        }

        public double getTotalIntensity() {
            return totalIntensity;
        }

        public int getPeakCount() {
            return peakCount;
        }
    }
}

