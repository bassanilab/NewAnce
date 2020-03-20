package newance.mzjava.ms.peaklist.peaktransformer;

import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TIntObjectHashMap;
import newance.mzjava.ms.peaklist.DelayedPeakProcessor;
import newance.mzjava.ms.peaklist.PeakAnnotation;

import java.util.Collections;
import java.util.List;

/**
 * As described in Eng et al., Journal Proteome Research, 2008, this transformation subtracts the average intensity
 * within a window of width 2*mvAvgWindow from every peak intensity. This method can produce negative intensities.
 *
 * Created with IntelliJ IDEA.
 * User: Markus
 * Date: 12.06.12
 * Time: 12:09
 * To change this template use File | Settings | File Templates.
 */
public class ContrastEnhancingTransformer<A extends PeakAnnotation>  extends DelayedPeakProcessor<A, A> {

    /** Window width for moving average process */
    private final double mvAvgWindow;

    /**
     *
     * @param mvAvgWindow  half width of sliding window
     */
    public ContrastEnhancingTransformer(double mvAvgWindow) {
        this.mvAvgWindow = mvAvgWindow;
    }

    @Override
    protected void processCached(TDoubleArrayList mzList, TDoubleArrayList intensityList, TIntObjectHashMap<List<A>> annotationMap) {

        double[] smoothedI = subtract(mzList,intensityList);

        int size =  mzList.size();
        for(int i = 0; i <size; i++) {

            List<A> annotations;
            if(annotationMap.contains(i))
                annotations = annotationMap.get(i);
            else
                annotations = Collections.emptyList();

            sink.processPeak(mzList.get(i), smoothedI[i], annotations);
        }
    }

    /**
     * Performs smoothing. For each peak, peaks within a window of +/-mvAvgWindow are collected, their average intensity
     * is calculated and subtracted from the current peak intensity.
     * @param mzList List of mz values
     * @param intensityList List of intensity values
     * @return processed intensity list
     */
    private double[] subtract(TDoubleArrayList mzList, TDoubleArrayList intensityList)
    {
        double[] intensitiesCopy = new double[intensityList.size()];

        int iMin, iMax,cnt;
        double avgI;
        iMax = iMin = cnt = 0;
        avgI = 0.0;
        for (int i = 0; i < mzList.size(); i++) {
            double mz = mzList.get(i);

            // adjust left border of sliding window
            while (iMin<mzList.size() && mzList.get(iMin) <  mz-mvAvgWindow) {
                avgI -= intensityList.get(iMin);
                cnt--;
                iMin++;
            }
            // adjust right border of sliding window
            while (iMax<mzList.size() && mzList.get(iMax) <= mz+mvAvgWindow) {
                avgI += intensityList.get(iMax);
                cnt++;
                iMax++;
            }

            if (cnt>0) intensitiesCopy[i] = intensityList.get(i) - avgI / cnt;
            else intensitiesCopy[i] = 0;
        }

        return intensitiesCopy;
    }
}
