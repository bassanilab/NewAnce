package newance.mzjava.ms.peaklist.peakfilter;


import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TIntObjectHashMap;
import newance.mzjava.ms.peaklist.DelayedPeakProcessor;
import newance.mzjava.ms.peaklist.PeakAnnotation;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * This filters applies a sliding window and keeps a peak only if it is among the N-highest peaks in this window.
 * The sliding window is centered around each peak with width +/- mzWindow. All peaks that are close (within +/- mzDist
 * of each peak) are also retained, which may be useful to keep isotopic distributions intact. To switch off this latter
 * option, set mzDist = 0.
 * @author mmuller
 * @author Oliver Horlacher
 */
public class NPeaksPerSlidingWindowFilter<A extends PeakAnnotation> extends DelayedPeakProcessor<A, A> {

    /**
     * Number of highest peaks within window *
     */
    private final int nrOfPeaks;
    /**
     * Half the width of the window *
     */
    private final double mzWindow;
    /**
     * Distance in m/z to highest peak, within which other peaks are retained
     */
    private final double mzDist;

    private int iMin, iMax;

    public NPeaksPerSlidingWindowFilter(final int nrOfPeaks,
                                        final double mzWindow,
                                        final double mzDist) {

        this.nrOfPeaks = nrOfPeaks;
        this.mzWindow = mzWindow;
        this.mzDist = mzDist;
    }

    public NPeaksPerSlidingWindowFilter(final int nrOfPeaks,
                                        final double mzWindow) {
        this(nrOfPeaks, mzWindow, 0.0);
    }


    protected void processCached(TDoubleArrayList mzList, TDoubleArrayList intensityList, TIntObjectHashMap<List<A>> annotationMap) {

        // first get all local top rank peaks
        int size =  mzList.size();
        int[] idxPass = new int[size];
        int[] idxStick = new int[size];
        int cnt1 = 0;
        int cnt2 = 0;
        iMin = iMax = 0;
        for(int i = 0; i <size; i++) {
            if (isPassingFilter(mzList, intensityList, i)) {
                idxPass[cnt1] = i;
                cnt1++;
            } else {
                idxStick[cnt2] = i;
                cnt2++;
            }
        }

        // then get all peak close to them
        int cnt3 = cnt1;
        if (mzDist>0 && cnt1>0 && cnt2>0) {
            // find peaks close to the N-highest peaks
            int j = 0;
            int curr = idxPass[j];
            // only check peaks that have not yet passed
            for (int k = 0; k < cnt2; k++) {
                int i = idxStick[k];
                if (mzList.get(i) < mzList.get(curr) - mzDist)
                    continue;

                if (mzList.get(i) <= mzList.get(curr) + mzDist) {
                    // within mzDist of high ranking peak
                    idxPass[cnt3] = i;
                    cnt3++;
                } else {
                    // outside mzDist of high ranking peak: update j
                    if (j < cnt1 - 1) {
                        j++;
                        curr = idxPass[j];
                        if (Math.abs(mzList.get(i) - mzList.get(curr)) <= mzDist) {
                            // check whether within mzDist of new peak j
                            idxPass[cnt3] = i;
                            cnt3++;
                        }
                    } else {
                        break;
                    }
                }
            }
        }

        Arrays.sort(idxPass,0,cnt3);

        for(int i = 0; i <cnt3; i++) {
            int curr = idxPass[i];
            List<A> annotations;
            if(annotationMap.contains(curr))
                annotations = annotationMap.get(curr);
            else
                annotations = Collections.emptyList();

            sink.processPeak(mzList.get(curr), intensityList.get(curr), annotations);
        }
    }

    /**
     * Selects nrOfPeaks highest peaks within sliding window. Returns the same
     * but processed spectrum object.
     *
     * @param mzList the mz's to process.
     * @param intensityList the intensities to process.
     * @param index    the peak index.
     * @return true if current peak is selected.
     */
    public boolean isPassingFilter(TDoubleArrayList mzList, TDoubleArrayList intensityList, int index) {
        final double mz = mzList.get(index);
        final double h = intensityList.get(index);

        // adjust left border of sliding window
        while (iMin<mzList.size() && mzList.get(iMin) <  mz-mzWindow) iMin++;
        // adjust right border of sliding window
        while (iMax<mzList.size() && mzList.get(iMax) <= mz+mzWindow) iMax++;
        iMax--;

        // pick N highest peak within window
        int nr = 0;
        for (int i = iMin; i <= iMax; i++) {
            if (intensityList.get(i) > h) {
                nr++;
            }
        }

        // store indices of peaks
        return nr < nrOfPeaks;
    }
}
