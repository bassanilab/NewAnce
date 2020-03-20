package newance.mzjava.ms.peaklist.peakfilter;

import newance.mzjava.ms.peaklist.AbstractPeakProcessor;
import newance.mzjava.ms.peaklist.IntervalList;
import newance.mzjava.ms.peaklist.PeakAnnotation;

import java.util.List;

/**
 * Class to select or exclude peaks in certain mz ranges.
 * <p/>
 * Created with IntelliJ IDEA.
 * User: Markus
 * Date: 13.06.12
 * Time: 09:32
 * To change this template use File | Settings | File Templates.
 */

public class MzRangeFilter<A extends PeakAnnotation> extends AbstractPeakProcessor<A, A> {

    /**
     * Object to store intervals
     */
    private final IntervalList intervalList;
    /**
     * flag indicating whether to include or exclude the peaks lying in the intervals
     */
    private final boolean include;

    /**
     * @param intervals list of mz intervals
     * @param include   true, if peak in intervals are included, false if they are excluded
     */
    public MzRangeFilter(IntervalList intervals, boolean include) {

        this.include = include;
        this.intervalList = intervals;
    }

    @Override
    public void processPeak(double mz, double intensity, List<A> annotations) {

        if (intervalList.contains(mz) == include) {

            sink.processPeak(mz, intensity, annotations);
        }
    }
}
