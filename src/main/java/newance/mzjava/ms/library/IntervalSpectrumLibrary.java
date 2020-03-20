package newance.mzjava.ms.library;

import com.google.common.base.Preconditions;
import newance.mzjava.ms.peaklist.Peak;
import newance.mzjava.ms.peaklist.PeakList;

import java.util.Collection;

/**
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class IntervalSpectrumLibrary<S extends PeakList> extends AbstractIntervalSpectrumLibrary<S> {

    private final double lowerBound, upperBound;

    public IntervalSpectrumLibrary(double lowerBound, double upperBound, Collection<S> spectra) {

        super(spectra);

        Preconditions.checkArgument(lowerBound < upperBound);
        this.lowerBound = lowerBound;
        this.upperBound = upperBound;
    }

    protected Peak newFromPeak(double mz, int charge) {

        return new Peak(mz + lowerBound, 0, charge);
    }

    protected Peak newToPeak(double mz, int charge) {

        return new Peak(mz + upperBound, Double.MAX_VALUE, charge);
    }
}
