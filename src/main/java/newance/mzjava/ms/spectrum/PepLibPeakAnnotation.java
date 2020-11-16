package newance.mzjava.ms.spectrum;

import com.google.common.base.Optional;

import static com.google.common.base.Preconditions.checkNotNull;

/**
 * Class to store annotations of consensus peaks. It stores statistical information about the consensus peaks.
 *
 * @author Oliver Horlacher
 * @version 1
 */
public class PepLibPeakAnnotation extends LibPeakAnnotation {

    private final PepFragAnnotation fragmentAnnotation;

    /**
     * Copy constructor
     *
     * @param src the annotation to copy
     */
    public PepLibPeakAnnotation(PepLibPeakAnnotation src) {

        super(src);
        PepFragAnnotation srcFrag = src.fragmentAnnotation;
        fragmentAnnotation = srcFrag != null ? srcFrag.copy() : null;
    }

    /**
     * Constructor
     *
     * @param mergedPeakCount number of peaks merged into consensus peak
     * @param mzStd           standard deviation of mz values
     * @param intensityStd    standard deviation of intensities
     */
    public PepLibPeakAnnotation(int mergedPeakCount, double mzStd, double intensityStd) {

        this(mergedPeakCount, mzStd, intensityStd, Optional.<PepFragAnnotation>absent());
    }

    public PepLibPeakAnnotation(int mergedPeakCount, double mzStd, double intensityStd, Optional<PepFragAnnotation> optFragmentAnnotation) {

        super(mergedPeakCount, mzStd, intensityStd);

        checkNotNull(optFragmentAnnotation);
        fragmentAnnotation = optFragmentAnnotation.orNull();
    }

    public Optional<PepFragAnnotation> getOptFragmentAnnotation() {

        return Optional.fromNullable(fragmentAnnotation);
    }

    @Override
    public PepLibPeakAnnotation copy() {

        return new PepLibPeakAnnotation(this);
    }

    @Override
    public boolean equals(Object o) {

        if (this == o) return true;
        if (!(o instanceof PepLibPeakAnnotation)) return false;
        if (!super.equals(o)) return false;

        PepLibPeakAnnotation that = (PepLibPeakAnnotation) o;

        return !(fragmentAnnotation != null ? !fragmentAnnotation.equals(that.fragmentAnnotation) : that.fragmentAnnotation != null);

    }

    @Override
    public int hashCode() {

        int result = super.hashCode();
        result = 31 * result + (fragmentAnnotation != null ? fragmentAnnotation.hashCode() : 0);
        return result;
    }

    @Override
    public String toString() {

        return "PepLibPeakAnnotation{" +
                "mergedPeakCount=" + getMergedPeakCount() +
                ", mzStd=" + getMzStd() +
                ", intensityStd=" + getIntensityStd() +
                ", fragmentAnnotation=" + fragmentAnnotation +
                '}';
    }

    @Override
    public int getCharge() {

        return fragmentAnnotation == null ? super.getCharge() : fragmentAnnotation.getCharge();
    }
}
