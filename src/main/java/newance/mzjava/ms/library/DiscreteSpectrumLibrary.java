package newance.mzjava.ms.library;

import newance.mzjava.mol.modification.Tolerance;
import newance.mzjava.ms.peaklist.Peak;
import newance.mzjava.ms.peaklist.PeakList;

import java.util.*;

import static com.google.common.base.Preconditions.checkArgument;
import static com.google.common.base.Preconditions.checkNotNull;

/**
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class DiscreteSpectrumLibrary<S extends PeakList> extends AbstractMapSpectrumLibrary<S> {

    private final double[] mzOffsets;
    private final Tolerance tolerance;

    public DiscreteSpectrumLibrary(double[] mzOffsets, Tolerance tolerance, Collection<S> spectra) {

        super(spectra);

        checkNotNull(mzOffsets);
        checkNotNull(tolerance);
        checkArgument(mzOffsets.length >= 1);

        Arrays.sort(mzOffsets);
        this.mzOffsets = mzOffsets;
        this.tolerance = tolerance;
    }

    @Override
    public void forEach(Peak precursor, Procedure<S> procedure) {

        final double mz = precursor.getMz();
        final int charge = precursor.getCharge();

        final Set<S> spectra = new HashSet<>();

        for (double offset : mzOffsets) {

            final Peak fromKey = new Peak(tolerance.getMin(mz + offset), 0, charge);
            final Peak toKey = new Peak(tolerance.getMax(mz + offset), Double.MAX_VALUE, charge);

            for (NavigableMap<Peak, S> spectraMap : maps) {

                final NavigableMap<Peak, S> withinTolerance = spectraMap.subMap(fromKey, true, toKey, true);
                for(S spectrum : withinTolerance.values()) {

                    if (spectra.add(spectrum)) {

                        procedure.execute(spectrum);
                    }
                }
            }
        }
    }
}
