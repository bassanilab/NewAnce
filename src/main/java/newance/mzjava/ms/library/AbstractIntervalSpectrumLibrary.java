package newance.mzjava.ms.library;


import newance.mzjava.ms.peaklist.Peak;
import newance.mzjava.ms.peaklist.PeakList;

import java.util.Collection;
import java.util.NavigableMap;

/**
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public abstract class AbstractIntervalSpectrumLibrary<S extends PeakList> extends AbstractMapSpectrumLibrary<S> {

    public AbstractIntervalSpectrumLibrary(Collection<S> spectra) {

        super(spectra);
    }

    @Override
    public void forEach(Peak precursor, Procedure<S> procedure) {

        final double mz = precursor.getMz();
        final int charge = precursor.getCharge();

        final Peak fromKey = newFromPeak(mz, charge);
        final Peak toKey = newToPeak(mz, charge);

        for (NavigableMap<Peak, S> spectraMap : maps) {

            final NavigableMap<Peak, S> withinTolerance = spectraMap.subMap(fromKey, true, toKey, true);
            for(S spectrum : withinTolerance.values()) {

                procedure.execute(spectrum);
            }
        }
    }

    protected abstract Peak newToPeak(double mz, int charge);

    protected abstract Peak newFromPeak(double mz, int charge);
}
