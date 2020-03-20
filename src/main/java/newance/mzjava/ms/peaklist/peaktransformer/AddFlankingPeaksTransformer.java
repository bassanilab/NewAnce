package newance.mzjava.ms.peaklist.peaktransformer;

import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TIntObjectHashMap;
import newance.mzjava.ms.peaklist.DelayedPeakProcessor;
import newance.mzjava.ms.peaklist.PeakAnnotation;

import java.util.Collections;
import java.util.List;

/**
 * Add flanking peaks to all peaks with an intensity higher than a threshold to a binned spectrum. Works only
 * for binned spectra. This is used to pre-process spectra for calculating the xcorr score.
 *
 * @author Markus Muller
 * @version 0.0
 */
public class AddFlankingPeaksTransformer<A extends PeakAnnotation>  extends DelayedPeakProcessor<A, A> {

    /** Relative intensity of flanking peaks */
    private final double relativeWeight;

    /** Minimal intensity of peak */
    private final double minIntensity;

    /**
     *
     * @param relativeWeight  relative intensity of flanking peaks
     * @param minIntensity  Minimal intensity of peak
     */
    public AddFlankingPeaksTransformer(double relativeWeight, double minIntensity) {
        this.relativeWeight = relativeWeight;
        this.minIntensity = minIntensity;
    }

    @Override
    protected void processCached(TDoubleArrayList mzList, TDoubleArrayList intensityList, TIntObjectHashMap<List<A>> annotationMap) {

        double[] smoothedI = addFlankingPeaks(mzList,intensityList);

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

    private double[] addFlankingPeaks(TDoubleArrayList mzList, TDoubleArrayList intensityList)
    {
        double[] intensitiesCopy = new double[intensityList.size()];

        int size = mzList.size();
        for (int i = 0; i <size; i++) {
            double intensity = intensityList.get(i);
            if (intensity>=minIntensity) {
                double addedIntensity = intensity*relativeWeight;
                if (i-1>=0 && addedIntensity>intensitiesCopy[i-1]) intensitiesCopy[i-1] = addedIntensity;
                intensitiesCopy[i] = intensity;
                if (i+1<size && addedIntensity>intensitiesCopy[i+1]) intensitiesCopy[i+1] = addedIntensity;
            }
        }

        return intensitiesCopy;
    }
}
