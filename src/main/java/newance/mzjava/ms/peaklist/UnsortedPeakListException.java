package newance.mzjava.ms.peaklist;

/**
 * @author fnikitin
 * Date: 3/7/13
 */
public class UnsortedPeakListException extends RuntimeException {

    private static final String DEFAULT_MESSAGE = "mz values are not sorted at index ";

    private final int peakIndex;

    public UnsortedPeakListException(int peakIndex) {

        super(DEFAULT_MESSAGE + peakIndex);
        this.peakIndex = peakIndex;
    }

    public UnsortedPeakListException(String message, int peakIndex) {

        super(message+" - "+ DEFAULT_MESSAGE + peakIndex);
        this.peakIndex = peakIndex;
    }

    public UnsortedPeakListException(String message, UnsortedPeakListException cause) {

        super(message, cause);
        this.peakIndex = cause.peakIndex;
    }

    public int getPeakIndex() {

        return peakIndex;
    }
}
