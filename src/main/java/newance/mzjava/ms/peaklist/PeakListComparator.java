package newance.mzjava.ms.peaklist;

import java.io.Serializable;
import java.util.Comparator;

/**
 * @author Markus Muller
 * @version 0.0
 */
public class PeakListComparator implements Comparator<PeakList>, Serializable {

    @Override
    public int compare(PeakList peakList, PeakList peakList2) {

        return peakList.getPrecursor().compareTo(peakList2.getPrecursor());
    }
}
