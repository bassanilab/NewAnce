package newance.mzjava.ms.library;

import newance.mzjava.ms.peaklist.Peak;
import newance.mzjava.ms.peaklist.PeakList;

import java.util.*;

/**
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public abstract class AbstractMapSpectrumLibrary<S extends PeakList> implements SpectrumLibrary<S> {

    protected final List<NavigableMap<Peak, S>> maps = new ArrayList<>();

    public AbstractMapSpectrumLibrary(Collection<S> spectra) {

        for (S spectrum : spectra) {

            add(spectrum.getPrecursor(), spectrum);
        }
    }

    private void add(Peak precursor, S spectrum) {

        for(NavigableMap<Peak, S> spectraMap : maps) {

            if (!spectraMap.containsKey(precursor)) {
                spectraMap.put(precursor, spectrum);
                return;
            }
        }

        NavigableMap<Peak, S> spectraMap = new TreeMap<>();
        spectraMap.put(precursor, spectrum);
        maps.add(spectraMap);
    }
}
