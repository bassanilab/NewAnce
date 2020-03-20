package newance.mzjava.ms.library;

import com.google.common.base.Preconditions;
import newance.mzjava.mol.modification.Tolerance;
import newance.mzjava.ms.peaklist.Peak;
import newance.mzjava.ms.peaklist.PeakList;

import java.util.Collection;

/**
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class DefaultSpectrumLibrary<S extends PeakList> extends AbstractIntervalSpectrumLibrary<S> {

    private final Tolerance tolerance;

    public DefaultSpectrumLibrary(Tolerance tolerance, Collection<S> spectra) {

        super(spectra);
        Preconditions.checkNotNull(tolerance);
        this.tolerance = tolerance;
    }

    protected Peak newFromPeak(double mz, int charge) {

        return new Peak(tolerance.getMin(mz), 0, charge);
    }

    protected Peak newToPeak(double mz, int charge) {

        return new Peak(tolerance.getMax(mz), Double.MAX_VALUE, charge);
    }
}
