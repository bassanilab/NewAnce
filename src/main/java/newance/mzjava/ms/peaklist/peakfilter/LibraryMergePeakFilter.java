package newance.mzjava.ms.peaklist.peakfilter;


import newance.mzjava.ms.peaklist.PeakAnnotation;
import newance.mzjava.ms.spectrum.LibPeakAnnotation;

/**
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class LibraryMergePeakFilter<A extends PeakAnnotation> extends AbstractMergePeakFilter<A, LibPeakAnnotation> {

    public LibraryMergePeakFilter(double mzMaxDiff, double maxMzClusterWidth, IntensityMode intCombMeth, int totSpectrumCount) {

        super(mzMaxDiff, maxMzClusterWidth, intCombMeth, totSpectrumCount);
    }

    @Override
    protected LibPeakAnnotation createPeakAnnotation(int peakCount, double mzStd, double intensityStd) {

        return new LibPeakAnnotation(peakCount, mzStd, intensityStd);
    }
}
