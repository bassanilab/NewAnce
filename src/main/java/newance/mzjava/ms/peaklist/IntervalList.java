package newance.mzjava.ms.peaklist;

import com.google.common.collect.Range;

import java.util.*;

/**
 * This class manages a list of intervals and implements a method to quickly check whether a given value lies within
 * the list of intervals or not. For more than 50 intervals, this class implements a quick access to the intervals
 * using a binned array. The bins are set to specific values indicating whether they are contained in a interval,
 * lie outside all intervals or contain a border of a interval. These methods accelerate the lookup of the intervals
 * only if the binWidth is much smaller than the average interval size. For less than 50 intervals, it turned out
 * that direct checks on the intervals is the faster way to check whether a value is contained in the intervals.
 *
 * @author Markus
 * @author nikitin
 */
public class IntervalList {

    /**
     * A array of bins from minVal to maxVal and width binWidth is used.
     */
    double minVal, maxVal, binWidth;

    /**
     * Binned array value is either 'e' (excluded), 'i' (included) or 'b' (boundary)
     */
    private char[] binArray;

    /**
     * Maps index in binned array to list of intervals that have a border in that bin.
     */
    private Map<Integer, List<Range<Double>>> intervalMap;

    /**
     * Number of intervals
     */
    private int nrIntervals;

    /**
     * List to store up to 50 intervals
     */
    private Range<Double>[] intervals;

    public IntervalList() {

        this.minVal = 0.0;
        this.maxVal = 4000.0;
        this.binWidth = 1.0;

        init();
    }

    /**
     * The binned array indexing the intervals goes from minVal to maxVal with bin size binWidth.
     *
     * @param minVal   Start value of the binned array
     * @param maxVal   End value of the binned array
     * @param binWidth Bind width of the binned array
     */
    public IntervalList(double minVal, double maxVal, double binWidth) {

        this.minVal = minVal;
        this.maxVal = maxVal;
        this.binWidth = binWidth;

        init();
    }

    /**
     * The binned array indexing the intervals goes from minVal to maxVal with bin size binWidth.
     *
     * @param minVal   Start value of the binned array
     * @param maxVal   End value of the binned array
     * @param binWidth Bind width of the binned array
     */
    public void setBinArrayBounds(double minVal, double maxVal, double binWidth) {

        this.minVal = minVal;
        this.maxVal = maxVal;
        this.binWidth = binWidth;

        init();
    }

    public double getMinVal() {

        return minVal;
    }

    public double getMaxVal() {

        return maxVal;
    }

    public double getBinWidth() {

        return binWidth;
    }

    /**
     * Add new interval to list. lBound and rBound are exchanged if lBound > rBound.
     *
     * @param lBound left bound of interval
     * @param rBound right bound of interval
     */
    public void addInterval(double lBound, double rBound) {

        Range<Double> interval;

        if (lBound <= rBound) {
            interval = Range.closed(lBound, rBound);
        } else {
            interval = Range.closed(rBound, lBound);
        }

        if (nrIntervals < 50) {
            intervals[nrIntervals] = interval;
            nrIntervals++;
        }


        int lBoundIdx = (int) Math.round((interval.lowerEndpoint() - minVal) / binWidth);
        int uBoundIdx = (int) Math.round((interval.upperEndpoint() - minVal) / binWidth);

        if (uBoundIdx < 0) uBoundIdx = 0;
        if (lBoundIdx < 0) lBoundIdx = 0;
        if (uBoundIdx > binArray.length - 1) uBoundIdx = binArray.length - 1;
        if (lBoundIdx > binArray.length - 1) lBoundIdx = binArray.length - 1;

        List<Range<Double>> il;
        if (binArray[lBoundIdx] != 'i') {
            binArray[lBoundIdx] = 'b'; // b for boundary

            // store intervals to check boundaries when testing
            // lower bound
            if (intervalMap.containsKey(lBoundIdx)) {
                il = intervalMap.get(lBoundIdx);
            } else {
                il = new ArrayList<Range<Double>>();
                intervalMap.put(lBoundIdx, il);
            }
            il.add(interval);
        }

        if (uBoundIdx != lBoundIdx && binArray[uBoundIdx] != 'i') {
            binArray[uBoundIdx] = 'b'; // b for boundary
            if (intervalMap.containsKey(uBoundIdx)) {
                il = intervalMap.get(uBoundIdx);
            } else {
                il = new ArrayList<Range<Double>>();
                intervalMap.put(uBoundIdx, il);
            }
            il.add(interval);
        }

        // all internal cells can be set to 'i'. they don't have to be checked later
        for (int j = lBoundIdx + 1; j < uBoundIdx; j++) binArray[j] = 'i';  // i for internal
    }

    /**
     * Checks whether value is included in any of the intervals
     *
     * @param value value to be tested
     * @return true if included, false otherwise
     */
    public final boolean contains(double value) {

        boolean match = false;
        if (nrIntervals < 50) {
            for (int i = 0; !match && i < nrIntervals; i++) {
                match = intervals[i].contains(value);
            }
        } else {
            int idx = (int) Math.round((value - minVal) / binWidth);
            if (idx < 0) idx = 0;
            if (idx > binArray.length - 1) idx = binArray.length - 1;

            if (binArray[idx] == 'i') return true;
            else if (binArray[idx] == 'e') return false;

            List<Range<Double>> il;
            // check boundaries
            if (intervalMap.containsKey(idx)) {
                il = intervalMap.get(idx);
                for (Range<Double> i : il) {
                    match = i.contains(value);
                    if (match) break;
                }
            }
        }

        return match;
    }

    private void init() {

        int nrBins = (int) Math.ceil((maxVal - minVal) / binWidth);

        binArray = new char[nrBins];
        Arrays.fill(binArray, 'e');  // e for external

        intervalMap = new HashMap<Integer, List<Range<Double>>>();

        nrIntervals = 0;
        //noinspection unchecked
        intervals = new Range[50];
    }

}
