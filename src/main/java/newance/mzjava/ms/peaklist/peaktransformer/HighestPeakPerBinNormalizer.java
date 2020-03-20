package newance.mzjava.ms.peaklist.peaktransformer;

import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TIntObjectHashMap;
import newance.mzjava.ms.peaklist.DelayedPeakProcessor;
import newance.mzjava.ms.peaklist.PeakAnnotation;

import java.util.Collections;
import java.util.List;

/**
 * Normalises a PeakList by dividing the PeakList into bins and scaling the peaks in each bin so that
 * the most intense peak has intensity = <code>maxBinIntensity</code>
 *
 * @author  Markus Muller
 * @version 1.0
 */
public class HighestPeakPerBinNormalizer<A extends PeakAnnotation>  extends DelayedPeakProcessor<A, A> {

    private final double binWidth;
    private final double maxIntensPerBin;

    /**
     * Construct a HighestPeakPerBinNormalizer.
     *
     * @param binWidth the with of the bin
     * @param maxBinIntensity the intensity to which the highest peak in each bin is to be scaled
     */
    public HighestPeakPerBinNormalizer(double binWidth, double maxBinIntensity) {

        this.binWidth = binWidth;
        this.maxIntensPerBin = maxBinIntensity;
    }

    @Override
    protected void processCached(TDoubleArrayList mzList, TDoubleArrayList intensityList, TIntObjectHashMap<List<A>> annotationMap) {

        TDoubleArrayList newIntens = normalize(mzList, intensityList);
        int size =  mzList.size();
        for(int i = 0; i <size; i++) {

            List<A> annotations;
            if(annotationMap.contains(i))
                annotations = annotationMap.get(i);
            else
                annotations = Collections.emptyList();

            sink.processPeak(mzList.get(i), newIntens.get(i), annotations);
        }
    }

    private TDoubleArrayList normalize(TDoubleArrayList mzList, TDoubleArrayList intensityList) {

        TDoubleArrayList newIntens = new TDoubleArrayList();

        int binIndex = 0;
        int iMin = 0;
        int iMax = 0;
        double maxIntens = 0.0;

        for (int i = 0; i < mzList.size(); i++) {
            int currBinIndex = (int) Math.floor(mzList.get(i)/binWidth);

            if (currBinIndex == binIndex) {
                if (intensityList.get(i) > maxIntens) {
                    maxIntens = intensityList.get(i);
                }
                iMax = i;
            } else {
                for (int j = iMin; j <= iMax; j++) {
                    if (Math.abs(maxIntens)>0.00001)
                        newIntens.add(maxIntensPerBin*intensityList.get(j)/maxIntens);
                    else
                        newIntens.add(0.0);
                }
                iMin = i;
                iMax = i;
                maxIntens = intensityList.get(i);
                binIndex = currBinIndex;
            }
        }
        for (int j = iMin; j <= iMax; j++) {
            if (Math.abs(maxIntens)>0.00001)
                newIntens.add(maxIntensPerBin*intensityList.get(j)/maxIntens);
            else
                newIntens.add(0.0);
        }

        return  newIntens;
    }

}
