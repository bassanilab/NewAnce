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
import gnu.trove.procedure.TDoubleProcedure;
import newance.mzjava.ms.peaklist.DelayedPeakProcessor;
import newance.mzjava.ms.peaklist.PeakAnnotation;

import java.util.ArrayList;
import java.util.List;

import static com.google.common.base.Preconditions.checkArgument;

/**
 * Class to bin a spectrum. Peaks will still be stored as a PeakList, but the mz values will be
 * regularly spaced on the centers of the bins.
 * @author mueller
 * @version 1.0
 */
public class BinnedSpectrumFilter<A extends PeakAnnotation>  extends DelayedPeakProcessor<A, A> {

    public enum IntensityMode {HIGHEST, SUM, AVG}

    private final double startMz, endMz, binWidth;
    private final IntensityMode intensityMode;

    private final TDoubleArrayList intensityBuffer;
    private final HighestIntensityProcedure highestIntensityProcedure = new HighestIntensityProcedure();
    private final SumIntensityProcedure sumIntensityProcedure = new SumIntensityProcedure();
    private final AverageIntensityProcedure averageIntensityProcedure = new AverageIntensityProcedure();

    /**
     * Constructor class for binned spectrum. Binned mz values will be at
     * startMz + binWidth/2 + i*binWidth, i>0 && i < (endMz-startMz)/binWidth
     *
     * @param startMz  start mz of the first bin
     * @param endMz    end mz of the last bin
     * @param binWidth  binWidth
     * @param intensityMode Mode for the calculation of centroid intensity
     *                      - IntensityMode.HIGHEST : set to highest intensity in bin
     *                      - IntensityMode.SUM : set to summed intensities in bin
     *                      - IntensityMode.AVG : set to average intensity in bin
     */
    public BinnedSpectrumFilter(double startMz, double endMz, double binWidth, IntensityMode intensityMode) {

        checkArgument(binWidth > 0);
        checkArgument(endMz > startMz);

        this.startMz = startMz;
        this.endMz = endMz;
        this.binWidth = binWidth;
        this.intensityMode = intensityMode;

        intensityBuffer = new TDoubleArrayList();
    }


    @Override
    protected void processCached(TDoubleArrayList mzList, TDoubleArrayList intensityList, TIntObjectHashMap<List<A>> annotationMap) {
        int maxIndex = (int)Math.floor((endMz - startMz) / binWidth);

        List<List<A>> annotations = new ArrayList<List<A>>(maxIndex);
        int i = 0;
        double halfBinWidth = binWidth/2.0;
        for (int j = 0; j < maxIndex; j++) {
            double binnedMz = startMz + halfBinWidth + j*binWidth;
            annotations.add(j,new ArrayList<A>());
            intensityBuffer.resetQuick();
            while (i<mzList.size() && mzList.get(i) < binnedMz+halfBinWidth) {
                if (mzList.get(i) >= binnedMz-halfBinWidth) {
                    intensityBuffer.add(intensityList.get(i));
                    List<A> peakAnnots = annotations.get(j);
                    if (annotationMap.contains(i))  {
                        peakAnnots.addAll(annotationMap.get(i));
                    }
                }
                i++;
            }
            sink.processPeak(binnedMz, getBinIntensity(), annotations.get(j));
        }
    }

    protected double getBinIntensity() {

        ResetableProcedure procedure;
        switch (intensityMode) {
            case HIGHEST:

                procedure = highestIntensityProcedure;
                break;
            case SUM:

                procedure = sumIntensityProcedure;
                break;
            case AVG:

                procedure = averageIntensityProcedure;
                break;
            default:
                throw new IllegalStateException("Cannot calc centroid intensity for " + intensityMode);
        }

        procedure.reset();
        intensityBuffer.forEach(procedure);
        return procedure.getResult();
    }

    private interface ResetableProcedure extends TDoubleProcedure {

        void reset();

        double getResult();
    }

    private static class HighestIntensityProcedure implements ResetableProcedure {

        private double max = 0;

        @Override
        public boolean execute(double value) {

            max = Math.max(max, value);

            return true;
        }

        public void reset() {

            max = 0;
        }

        @Override
        public double getResult() {

            return max;
        }
    }

    private static class SumIntensityProcedure implements ResetableProcedure {

        private double sum = 0;

        @Override
        public boolean execute(double value) {

            sum += value;
            return true;
        }

        @Override
        public double getResult() {

            return sum;
        }

        public void reset(){

            sum = 0;
        }
    }

    private static class AverageIntensityProcedure implements ResetableProcedure {

        private int count = 0;
        private double sum = 0;

        @Override
        public void reset() {

            count = 0;
            sum = 0;
        }

        @Override
        public double getResult() {

            double average = sum;
            if (count > 0) {

                average = sum/count;
            }
            return average;
        }

        @Override
        public boolean execute(double value) {

            sum =+ value;
            return true;
        }
    }
}

