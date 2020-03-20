package newance.mzjava.ms.consensus;


import newance.mzjava.ms.peaklist.PeakAnnotation;
import newance.mzjava.ms.peaklist.PeakList;
import newance.mzjava.ms.peaklist.peakfilter.AbstractMergePeakFilter;
import newance.mzjava.ms.spectrum.LibPeakAnnotation;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * @author Markus Muller
 * @version 0.0
 */
public abstract class AbstractConsensusSpectrumBuilder<A extends LibPeakAnnotation, B extends AbstractConsensusSpectrumBuilder<A, B,C >, C extends ConsensusSpectrum<A>> {
    protected final List<PeakList<PeakAnnotation>> spectra;
    protected double peakFraction;
    protected int minPeakCount;
    protected double fragMzTolerance;
    protected double maxMzClusterWidth;
    protected AbstractMergePeakFilter.IntensityMode intCombMeth;
    protected int initialCapacity;
    protected int msLevel;

    protected AbstractConsensusSpectrumBuilder() {
        spectra = new ArrayList<>();
        peakFraction= 0.2;
        minPeakCount = 2;
        fragMzTolerance = 0.3;
        maxMzClusterWidth = 0.4;
        intCombMeth = AbstractMergePeakFilter.IntensityMode.MEAN_ALL_INTENSITY;
        msLevel = 2;
    }

    public B spectra(Collection<? extends PeakList<PeakAnnotation>> spectra) {
        this.spectra.addAll(spectra);
        return (B)this;
    }

    public B initialCapacity(int initialCapacity) {
        this.initialCapacity = initialCapacity;
        return (B)this;
    }

    public B fragMzTolerance(double fragMzTolerance) {
        this.fragMzTolerance = fragMzTolerance;
        return (B) this;
    }

    public B intensityCombMethod(AbstractMergePeakFilter.IntensityMode intCombMethod) {
        this.intCombMeth = intCombMethod;
        return (B) this;
    }

    public B maxMzClusterWidth(double maxMzClusterWidth) {
        this.maxMzClusterWidth = maxMzClusterWidth;
        return (B) this;
    }

    public B setPeakFilterParams(double peakFraction,int minPeakCount) {
        this.peakFraction = peakFraction;
        this.minPeakCount = minPeakCount;
        return (B) this;
    }

    public B msLevel(int msLevel) {
        this.msLevel = msLevel;
        return (B) this;
    }

    protected abstract void calculatePrecursor(ConsensusSpectrum<? extends LibPeakAnnotation> consensus);
    public abstract C build();

}
