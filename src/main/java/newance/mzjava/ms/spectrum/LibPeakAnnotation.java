package newance.mzjava.ms.spectrum;

import newance.mzjava.ms.peaklist.PeakAnnotation;

/**
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class LibPeakAnnotation implements PeakAnnotation {
    /**
     * number of peaks merged into consensus peak
     */
    private final int mergedPeakCount;
    /**
     * standard deviation of mz values
     */
    private final double mzStd;
    /**
     * standard deviation of intensities
     */
    private final double intensityStd;

    public LibPeakAnnotation(int mergedPeakCount, double mzStd, double intensityStd) {

        this.mergedPeakCount = mergedPeakCount;
        this.mzStd = mzStd;
        this.intensityStd = intensityStd;
    }

    public LibPeakAnnotation(LibPeakAnnotation src) {

        this.mergedPeakCount = src.mergedPeakCount;
        this.mzStd = src.mzStd;
        this.intensityStd = src.intensityStd;
    }

    public int getMergedPeakCount() {

        return mergedPeakCount;
    }

    public double getMzStd() {

        return mzStd;
    }

    public double getIntensityStd() {

        return intensityStd;
    }

    @Override
    public LibPeakAnnotation copy() {

        return new LibPeakAnnotation(this);
    }

    @Override
    public boolean equals(Object o) {

        if (this == o) return true;
        if (!(o instanceof LibPeakAnnotation)) return false;

        LibPeakAnnotation that = (LibPeakAnnotation) o;

        return Double.compare(that.intensityStd, intensityStd) == 0 &&
                mergedPeakCount == that.mergedPeakCount &&
                Double.compare(that.mzStd, mzStd) == 0;
    }

    @Override
    public int hashCode() {

        int result;
        long temp;
        result = mergedPeakCount;
        temp = mzStd != +0.0d ? Double.doubleToLongBits(mzStd) : 0L;
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        temp = intensityStd != +0.0d ? Double.doubleToLongBits(intensityStd) : 0L;
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        return result;
    }

    @Override
    public int getCharge() {

        return 0;
    }

    @Override
    public String toString() {

        return "LibraryPeakAnnotation{" +
                "mergedPeakCount=" + mergedPeakCount +
                ", mzStd=" + mzStd +
                ", intensityStd=" + intensityStd +
                '}';
    }
}
