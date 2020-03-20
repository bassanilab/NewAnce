package newance.mzjava.ms.peaklist.peakfilter;

import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TIntObjectHashMap;
import newance.mzjava.ms.peaklist.DelayedPeakProcessor;
import newance.mzjava.ms.peaklist.PeakAnnotation;
import newance.mzjava.ms.spectrum.LibPeakAnnotation;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Class to merge peaks that are close together. This class is used to build consensus spectra.
 *
 * @author Markus
 */
public abstract class AbstractMergePeakFilter<IN extends PeakAnnotation, OUT extends LibPeakAnnotation> extends DelayedPeakProcessor<IN, OUT> {

    public enum IntensityMode {MEAN_ALL_INTENSITY, MEAN_POS_INTENSITY, SUM_INTENSITY}

    /**
     * Indices to bound peak groups to be merged
     * iLowBound, iUpBound : lower and upper bound
     * count : peak count
     */
    private int iLowBound, iUpBound,count;

    /**
     * Doubles to store peak group statistics
     */
    private double mzStd,mzMean,iStd,iMean;
    /**
     * Maximum distance between consecutive mz values in same peak group
     */
    private final double mzMaxDiff;

    /**
     * Maximal allowed width for an mz cluster. If the range of a cluster is larger than this value,
     * the cluster is split into sub-clusters.
     */
    private final double maxMzClusterWidth;

    /**
     * Method how to combine intensities
     */
    private final IntensityMode intCombMeth;

    /**
     * Total number of spectra merged
     */
    private final int totSpectrumCount;

    /**
     * Variable to stare temporary maximal density
     */
    private double maxDensity;

    /**
     * Constructor
     * @param mzMaxDiff maximal difference in mz between two consecutive peaks in the same peak group
     * @param maxMzClusterWidth maximal difference in mz between lightest and heaviest peak in the same peak group
     * @param intCombMeth indicates whether to sum or average intensities
     * @param totSpectrumCount total number of spectra merged or maximal peak count
     */
    public AbstractMergePeakFilter(double mzMaxDiff, double maxMzClusterWidth, IntensityMode intCombMeth, int totSpectrumCount) {
        this.mzMaxDiff = mzMaxDiff;
        this.maxMzClusterWidth = maxMzClusterWidth;
        this.intCombMeth = intCombMeth;
        this.totSpectrumCount = totSpectrumCount;
    }

    /**
     * Return the current lower bound peak index
     *
     * @return the current lower bound peak index
     */
    public int getLowerBound() {

        return iLowBound;
    }

    /**
     * Return the current upper bound peak index
     *
     * @return the current upper bound peak index
     */
    public int getUpperBound() {

        return iUpBound;
    }

    private boolean hasNextCluster(TDoubleArrayList mzList)
    {
        return iUpBound<mzList.size()-1;
    }

    private void nextCluster(TDoubleArrayList mzList)
    {
        iUpBound++;
        iLowBound = iUpBound;

        while (iUpBound<mzList.size()-1 && mzList.get(iUpBound+1)-mzList.get(iUpBound)<mzMaxDiff) {
            iUpBound++;
        }

        count = iUpBound-iLowBound+1;

        if (mzList.get(iUpBound)-mzList.get(iLowBound)>maxMzClusterWidth) {
            checkCluster(mzList);
        }
    }

    // Split cluster if the mz density has a pronounced minimum
    private void checkCluster(TDoubleArrayList mzList)
    {
        double[] density = calcMassDensity(mzList);
        int cutIdx = -1;

        int iMax = iLowBound;
        int iMin = iLowBound;
        double mzMin = mzList.get(iLowBound);
        double mzMax = mzList.get(iUpBound);
        double d = 0.0;
        for (int j = iLowBound; j <= iUpBound; j++) {
            double mz =  mzList.get(j);

            if (mz-mzMin <= mzMaxDiff) continue;
            if (mzMax-mz <= mzMaxDiff) break;

            // adjust left border of sliding window
            while (iMin<=iUpBound && mzList.get(iMin)<mz-mzMaxDiff) iMin++;
            // adjust right border of sliding window
            while (iMax<=iUpBound && mzList.get(iMax)<=mz+mzMaxDiff) iMax++;
            iMax--;

            d = density[j-iLowBound];
            boolean isMin = true;
            for (int i=iMin;isMin && i<=iMax;i++) {
                isMin = d<=density[i-iLowBound];
            }

            if (isMin) {
                cutIdx = j;
                break;
            }
        }

        if (cutIdx>=0 && (maxDensity-d)/maxDensity>0.1) {
            iUpBound = cutIdx;
            count = iUpBound-iLowBound+1;
        }
    }

    // Calculates mass density between iUpBound and iLowBound. Density values are calculated as weighted sums, where the
    // weights are linearly decreasing with distance.
    private double[] calcMassDensity(TDoubleArrayList mzList)
    {
        double[] density = new double[count];
        maxDensity = 0.0;
        double a = 0.5/mzMaxDiff;
        Arrays.fill(density,0.0);
        int iMax = iLowBound;
        int iMin = iLowBound;
        for (int j = iLowBound; j <= iUpBound; j++) {
            double mz =  mzList.get(j);
            // adjust left border of sliding window
            while (iMin<=iUpBound && mzList.get(iMin)<mz-mzMaxDiff) iMin++;
            // adjust right border of sliding window
            while (iMax<=iUpBound && mzList.get(iMax)<=mz+mzMaxDiff) iMax++;
            iMax--;

            double d = 0.0;
            for (int i=iMin;i<=iMax;i++) {
                d += Math.max(1.0-Math.abs(mz-mzList.get(i))*a,0.0);
            }
            density[j-iLowBound] = d;
            if (d>maxDensity) maxDensity = d;
        }

        return density;
    }

    private void mergePeaks(TDoubleArrayList mzList, TDoubleArrayList intensityList) {
        iMean = 0.0;
        mzMean = 0.0;
        for (int j = iLowBound; j <= iUpBound; j++) {
            double h = intensityList.get(j);
            iMean += h;
            mzMean += h * mzList.get(j);
        }

        // calc mean mz and intensity
        count = iUpBound-iLowBound+1;
        if (iMean > 0) {
            mzMean = mzMean / iMean;
        } else {
            mzMean = mzMean /count;
        }

        double meanInt = iMean / count;

        switch (intCombMeth) {

            case MEAN_POS_INTENSITY:
                iMean = meanInt;
                break;
            case MEAN_ALL_INTENSITY:
                iMean = iMean / Math.max(totSpectrumCount,count);
                break;
            case SUM_INTENSITY:
                // nothing to do here
                break;
        }

        // calc standard deviation for mz and intensity
        mzStd = 0.0;
        iStd = 0.0;
        for (int j = iLowBound; j <= iUpBound; j++) {
            double tmp =  mzList.get(j)-mzMean;
            mzStd += tmp*tmp;
            tmp = intensityList.get(j)-meanInt;
            iStd += tmp*tmp;
        }
        mzStd = Math.sqrt(mzStd/count);
        iStd = Math.sqrt(iStd/count);
    }

    protected void processCached(TDoubleArrayList mzList, TDoubleArrayList intensityList, TIntObjectHashMap<List<IN>> annotationMap) {

        iLowBound = -1;
        iUpBound = -1;
        while (hasNextCluster(mzList)) {
            nextCluster(mzList);
            mergePeaks(mzList, intensityList);

            OUT annotation = createPeakAnnotation(count, mzStd, iStd);
            sink.processPeak(mzMean, iMean, Collections.singletonList(annotation));
        }
    }

    protected abstract OUT createPeakAnnotation(int peakCount, double mzStd, double intensityStd);
}
