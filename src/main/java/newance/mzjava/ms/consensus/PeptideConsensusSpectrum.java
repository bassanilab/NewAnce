package newance.mzjava.ms.consensus;

import com.google.common.base.Optional;
import newance.mzjava.mol.Peptide;
import newance.mzjava.mol.PeptideFragmenter;
import newance.mzjava.mol.modification.AbsoluteTolerance;
import newance.mzjava.mol.modification.Tolerance;
import newance.mzjava.ms.peaklist.*;
import newance.mzjava.ms.peaklist.peakfilter.AbstractMergePeakFilter;
import newance.mzjava.ms.spectrum.*;

import java.net.URI;
import java.util.*;

import static com.google.common.base.Preconditions.checkNotNull;

/**
 * @author Oliver Horlacher
 * @version 1.0
 */
public class PeptideConsensusSpectrum extends ConsensusSpectrum<PepLibPeakAnnotation> {

    public enum Status {UNKNOWN, SINGLETON, NORMAL, INQUORATE_UNCONFIRMED, CONFLICTING_ID, IMPURE, DECOY}

    private Peptide peptide;
    private final Set<String> proteinAccessionNumbers = new LinkedHashSet<>();

    private String comment = "";
    private Status status = Status.UNKNOWN;

    private URI spectrumSource;

    private Optional<RetentionTimeDiscrete> retentionTime = Optional.absent();

    /**
     * Construct a LibrarySpectrum for the <code>peptide</code>
     *
     * @param peptide the peptide
     */
    public PeptideConsensusSpectrum(Peptide peptide) {

        this(peptide, 0, Precision.DOUBLE, new HashSet<UUID>());
    }

    /**
     * Construct a LibrarySpectrum for the <code>peptide</code>
     *
     * @param peptide   the peptide
     * @param memberIds UUIDs of member spectra
     */
    public PeptideConsensusSpectrum(Peptide peptide, Set<UUID> memberIds) {

        this(peptide, 0, Precision.DOUBLE, memberIds);
    }

    /**
     * Construct a Spectrum that has a precision of <code>precision</code>
     *
     * @param peptide   the peptide
     * @param precision the precision
     */
    public PeptideConsensusSpectrum(Peptide peptide, Precision precision, Set<UUID> memberIds) {

        this(peptide, 0, precision, memberIds);
    }

    /**
     * Construct a Spectrum that has a precision of <code>precision</code>
     *
     * @param peptide   the peptide
     * @param precision the precision
     */
    public PeptideConsensusSpectrum(Peptide peptide, Precision precision) {

        this(peptide, 0, precision, new HashSet<UUID>());
    }

    /**
     * Construct a Spectrum that has an initial capacity of <code>initialCapacity</code> and a precision
     * of <code>precision</code>
     *
     * @param peptide         the peptide
     * @param initialCapacity the initial capacity
     * @param precision       the precision
     */
    public PeptideConsensusSpectrum(Peptide peptide, int initialCapacity, Precision precision, Set<UUID> memberIds) {

        super(initialCapacity, precision, memberIds);

        checkNotNull(peptide);

        this.peptide = peptide;

        spectrumSource = URIBuilder.UNDEFINED_URI;

        setMsLevel(2);
    }

    public PeptideConsensusSpectrum(Peptide peptide, Precision precision, Set<UUID> memberIDs, Set<String> proteinAccessionNumbers) {

        this(peptide, precision, memberIDs);

        this.proteinAccessionNumbers.addAll(proteinAccessionNumbers);
    }

    public PeptideConsensusSpectrum(Peptide peptide, Precision precision, Set<UUID> memberIDs, Set<String> proteinAccessionNumbers, RetentionTimeDiscrete retentionTimeDiscrete) {

        this(peptide, precision, memberIDs);

        this.proteinAccessionNumbers.addAll(proteinAccessionNumbers);
        this.retentionTime = Optional.of(retentionTimeDiscrete);
    }

    /**
     * Construct a Spectrum that has an initial capacity of <code>initialCapacity</code> and a precision
     * of <code>precision</code>
     *
     * @param peptide         the peptide
     * @param initialCapacity the initial capacity
     * @param precision       the precision
     */
    public PeptideConsensusSpectrum(Peptide peptide, int initialCapacity, Precision precision) {

        super(initialCapacity, precision, new HashSet<UUID>());

        checkNotNull(peptide);

        this.peptide = peptide;

        spectrumSource = URIBuilder.UNDEFINED_URI;

        setMsLevel(2);
    }

    /**
     * Copy constructor
     *
     * @param src the LibrarySpectrum to copy
     */
    protected PeptideConsensusSpectrum(PeptideConsensusSpectrum src, PeakProcessor<PepLibPeakAnnotation, PepLibPeakAnnotation> peakProcessor) {

        super(src, peakProcessor);

        peptide = src.peptide;
        comment = src.comment;
        status = src.status;
        proteinAccessionNumbers.addAll(src.proteinAccessionNumbers);
        spectrumSource = src.spectrumSource;
        retentionTime = src.retentionTime;
    }

    protected PeptideConsensusSpectrum(PeptideConsensusSpectrum src, PeakProcessorChain<PepLibPeakAnnotation> peakProcessorChain) {

        super(src, peakProcessorChain);

        peptide = src.peptide;
        comment = src.comment;
        status = src.status;
        proteinAccessionNumbers.addAll(src.proteinAccessionNumbers);
        spectrumSource = src.spectrumSource;
        retentionTime = src.retentionTime;
    }

    public static BuilderOld builder() {

        return new BuilderOld();
    }

    @Override
    public boolean equals(Object o) {

        if (this == o) return true;
        if (!(o instanceof PeptideConsensusSpectrum)) return false;
        if (!super.equals(o)) return false;

        PeptideConsensusSpectrum that = (PeptideConsensusSpectrum) o;

        return Objects.equals(peptide, that.peptide) && Arrays.equals(getPrecursor().getChargeList(), that.getPrecursor().getChargeList());
    }

    @Override
    public int hashCode() {

        int result = super.hashCode();
        result = 31 * result + Objects.hashCode(peptide);
        result = 31 * result + Arrays.hashCode(getPrecursor().getChargeList());
        return result;
    }

    @Override
    public PeptideConsensusSpectrum copy(PeakProcessor<PepLibPeakAnnotation, PepLibPeakAnnotation> peakProcessor) {

        return new PeptideConsensusSpectrum(this, peakProcessor);
    }

    @Override
    public PeptideConsensusSpectrum copy(PeakProcessorChain<PepLibPeakAnnotation> peakProcessorChain) {

        return new PeptideConsensusSpectrum(this, peakProcessorChain);
    }

    /**
     * Return the peptide
     *
     * @return the peptide
     */
    public Peptide getPeptide() {

        return peptide;
    }

    //review this breaks the equals method
    public void setPeptide(Peptide peptide) {

        this.peptide = peptide;
    }

    /**
     * Add a protein accession number to the accession number list of this LibrarySpectrum
     *
     * @param acc the protein accession number
     */
    public void addProteinAccessionNumbers(String... acc) {

        Collections.addAll(proteinAccessionNumbers, acc);
    }

    /**
     * Add all the accession numbers from <code>accessionNumbers</code> to the accession number list of this
     * LibrarySpectrum
     *
     * @param accessionNumbers the collection containing the accession numbers to add
     */      //review do we really need this
    public void addProteinAccessionNumbers(Collection<String> accessionNumbers) {

        proteinAccessionNumbers.addAll(accessionNumbers);
    }

    /**
     * Clear the accession number list //review should this class be mutable?
     */
    public void clearProteinAccessionNumber() {

        proteinAccessionNumbers.clear();
    }

    /**
     * Return an unmodifiable list containing the protein accession numbers associated with the peptide
     * tha this LibrarySpectrum represents.
     *
     * @return an unmodifiable list containing the protein accession numbers associated with the peptide
     * tha this LibrarySpectrum represents
     */
    public Set<String> getProteinAccessionNumbers() {

        return Collections.unmodifiableSet(proteinAccessionNumbers);
    }

    @Override
    public String toString() {

        return "LibrarySpectrum{" +
                "peptide=" + peptide +
                " precursor=" + getPrecursor() +
                '}';
    }

    /**
     * Set the status
     *
     * @param status the status
     */
    public void setStatus(Status status) {

        this.status = status;
    }

    /**
     * Return the comment
     *
     * @return the comment
     */
    public String getComment() {

        return comment;
    }

    /**
     * Set the comment
     *
     * @param comment the comment   //review make immutable?
     */
    public void setComment(String comment) {

        this.comment = comment;
    }

    /**
     * Returns this LibrarySpectrum's status.
     *
     * @return this LibrarySpectrum's status
     */
    public Status getStatus() {

        return status;
    }

    public URI getSpectrumSource() {

        return spectrumSource;
    }

    public void setSpectrumSource(URI spectrumSource) {

        this.spectrumSource = spectrumSource;
    }

    public void setRetentionTime(Optional<RetentionTimeDiscrete> retentionTime) {
        this.retentionTime = retentionTime;
    }

    public Optional<RetentionTimeDiscrete> getRetentionTime() {
        return retentionTime;
    }

    /**
     * @author Markus Muller
     * @version 0.0
     */
    public static class BuilderOld extends AbstractConsensusSpectrumBuilder<PepLibPeakAnnotation, BuilderOld, PeptideConsensusSpectrum> {

        private Peptide peptide;
        private Precision precision;
        private PeptideFragmenter fragmenter;
        private URI spectrumSource;

        protected BuilderOld() {

            peptide = null;
            precision = Precision.DOUBLE;
            fragmenter = null;
        }


        public BuilderOld peptide(Peptide peptide) {

            this.peptide = peptide;
            return this;
        }

        public BuilderOld precision(Precision precision) {

            this.precision = precision;
            return this;
        }

        public BuilderOld fragmenter(PeptideFragmenter fragmenter) {

            this.fragmenter = fragmenter;
            return this;
        }

        public BuilderOld uri(URI spectrumSource) {

            this.spectrumSource = spectrumSource;
            return this;
        }

        @Override
        public PeptideConsensusSpectrum build() {

            checkNotNull(peptide);
            checkNotNull(spectra);
            checkNotNull(fragMzTolerance);
            checkNotNull(fragmenter);

            Set<UUID> memberIDs = new HashSet<>();
            if (spectra.isEmpty())
                return new PeptideConsensusSpectrum(peptide, precision, memberIDs);

            for (PeakList<PeakAnnotation> spectrum : spectra) {
                memberIDs.add(spectrum.getId());
            }

            PeptideConsensusSpectrum consensus = new PeptideConsensusSpectrum(peptide, precision, memberIDs);
            calculatePrecursor(consensus);
            consensus.setSpectrumSource(spectrumSource);

            int charge = consensus.getPrecursor().getCharge();
            PeptideSpectrum theoreticalSpectrum = fragmenter.fragment(peptide, charge);

            PeptideFragmentAnnotator fragmentAnnotator = new PeptideFragmentAnnotator(new AbsoluteTolerance(fragMzTolerance), theoreticalSpectrum);

            MergePeakFilter mergePeakFilter = new MergePeakFilter(fragMzTolerance, maxMzClusterWidth, intCombMeth, spectra.size());

            PeptideConsensusPeakSink peakFilter = new PeptideConsensusPeakSink(consensus, peakFraction, minPeakCount);
            fragmentAnnotator.setSink(peakFilter);

            if (spectra.size() > 1) {
                PeakListMerger<PeakAnnotation> peakListMerger = new PeakListMerger<>();
                mergePeakFilter.setSink(fragmentAnnotator);
                peakListMerger.setSink(mergePeakFilter);
                peakListMerger.merge(spectra);
            } else {
                PeakList spectrum = spectra.get(0);
                consensus.ensureCapacity(spectrum.size());
                for (int i = 0; i < spectrum.size(); i++) {
                    fragmentAnnotator.processPeak(spectrum.getMz(i), spectrum.getIntensity(i), Collections.singletonList(mergePeakFilter.createPeakAnnotation(1, 0.0, 0.0)));
                }
                consensus.trimToSize();
            }

            return consensus;
        }

        protected void calculatePrecursor(ConsensusSpectrum<? extends LibPeakAnnotation> consensus) {

            double meanIntensity = 0.0;
            double meanMz = 0.0;

            if (!spectra.isEmpty()) {

                int n = 0;

                for (PeakList pl : spectra) {
                    n++;
                    double delta = pl.getPrecursor().getMz() - meanMz;
                    meanMz += delta / n;
                    delta = pl.getPrecursor().getIntensity() - meanIntensity;
                    meanIntensity += delta / n;
                }

                int charge = spectra.get(0).getPrecursor().getCharge();
                Peak precursor = new Peak(peptide.calculateMz(charge), meanIntensity, charge);

                consensus.setPrecursor(precursor);
                consensus.setPrecursorStats(meanMz, 0.0);
            }

        }

        private class MergePeakFilter extends AbstractMergePeakFilter<PeakAnnotation, PepLibPeakAnnotation> {

            public MergePeakFilter(double mzMaxDiff, double maxMzClusterWidth, IntensityMode intCombMeth, int totSpectrumCount) {

                super(mzMaxDiff, maxMzClusterWidth, intCombMeth, totSpectrumCount);
            }

            @Override
            protected PepLibPeakAnnotation createPeakAnnotation(int peakCount, double mzStd, double intensityStd) {

                return new PepLibPeakAnnotation(peakCount, mzStd, intensityStd);
            }
        }
    }

    public static interface ConsensusParametersSetter {

        AnnotationParametersSetter setConsensusParameters(double maxMzDiff, double maxMzClusterWidth, AbstractMergePeakFilter.IntensityMode intensityMode);
    }

    public static interface AnnotationParametersSetter {

        FilterParametersSetter setAnnotationParameters(Tolerance fragmentTolerance, PeptideFragmenter fragmenter);
    }

    public static interface FilterParametersSetter {

        Builder setFilterParameters(double peakFraction, int minPeakCount);
    }

    public static ConsensusParametersSetter builder(Precision precision, URI spectrumSource) {

        return new Builder(precision, spectrumSource);
    }

    public static class Builder implements ConsensusParametersSetter, AnnotationParametersSetter, FilterParametersSetter {

        protected final Precision precision;
        protected final URI spectrumSource;

        protected double maxMzDiff;
        protected double maxMzClusterWidth;
        protected AbstractMergePeakFilter.IntensityMode intensityMode;

        protected Tolerance fragmentTolerance;
        protected PeptideFragmenter fragmenter;

        protected int minPeakCount;
        protected double peakFraction;

        protected Builder(Precision precision, URI spectrumSource) {

            this.precision = precision;
            this.spectrumSource = spectrumSource;
        }

        @Override
        public FilterParametersSetter setAnnotationParameters(Tolerance fragmentTolerance, PeptideFragmenter fragmenter) {

            checkNotNull(fragmentTolerance);
            checkNotNull(fragmenter);

            this.fragmentTolerance = fragmentTolerance;
            this.fragmenter = fragmenter;
            return this;
        }

        @Override
        public AnnotationParametersSetter setConsensusParameters(double maxMzDiff, double maxMzClusterWidth, AbstractMergePeakFilter.IntensityMode intensityMode) {

            checkNotNull(intensityMode);

            this.maxMzDiff = maxMzDiff;
            this.maxMzClusterWidth = maxMzClusterWidth;
            this.intensityMode = intensityMode;
            return this;
        }

        @Override
        public Builder setFilterParameters(double peakFraction, int minPeakCount) {

            this.minPeakCount = minPeakCount;
            this.peakFraction = peakFraction;
            return this;
        }

        public <A extends PeakAnnotation, S extends PeakList<A>> PeptideConsensusSpectrum buildConsensus(int charge, Peptide peptide, Collection<S> spectra, Set<String> proteinAccessionNumbers, RetentionTimeList retentionTimes) {

            PeptideConsensusSpectrum consensusSpectrum = buildConsensus(charge,peptide,spectra,proteinAccessionNumbers);

            consensusSpectrum.setRetentionTime(calculateRetentionTime(retentionTimes));

            return consensusSpectrum;
        }

        public <A extends PeakAnnotation, S extends PeakList<A>> PeptideConsensusSpectrum buildConsensus(int charge, Peptide peptide, Collection<S> spectra, Set<String> proteinAccessionNumbers) {

            checkNotNull(peptide);
            checkNotNull(spectra);

            Set<UUID> memberIDs = new HashSet<>();
            if (spectra.isEmpty())
                throw new IllegalStateException("Attempting to build a consensus from spectra collection that is empty");

            for (S spectrum : spectra) {

                memberIDs.add(spectrum.getId());
            }

            PeptideConsensusSpectrum consensus = new PeptideConsensusSpectrum(peptide, precision, memberIDs, proteinAccessionNumbers);
            consensus.setPrecursor(calculatePrecursor(charge, peptide, spectra));
            consensus.setSpectrumSource(spectrumSource);

            PeptideSpectrum theoreticalSpectrum = fragmenter.fragment(peptide, charge);

            PeptideFragmentAnnotator fragmentAnnotator = new PeptideFragmentAnnotator(fragmentTolerance, theoreticalSpectrum);

            MergePeakFilter<A> mergePeakFilter = new MergePeakFilter<>(maxMzDiff, maxMzClusterWidth, intensityMode, spectra.size());

            PeptideConsensusPeakSink peakFilter = new PeptideConsensusPeakSink(consensus, peakFraction, minPeakCount);
            fragmentAnnotator.setSink(peakFilter);

            if (spectra.size() > 1) {

                PeakListMerger<A> peakListMerger = new PeakListMerger<>();
                mergePeakFilter.setSink(fragmentAnnotator);
                peakListMerger.setSink(mergePeakFilter);
                peakListMerger.merge(spectra);
            } else {

                PeakList spectrum = spectra.iterator().next();
                consensus.ensureCapacity(spectrum.size());
                for (int i = 0; i < spectrum.size(); i++) {
                    fragmentAnnotator.processPeak(spectrum.getMz(i), spectrum.getIntensity(i), Collections.singletonList(mergePeakFilter.createPeakAnnotation(1, 0.0, 0.0)));
                }
                consensus.trimToSize();
            }

            consensus.setStatus(spectra.size() == 1 ? PeptideConsensusSpectrum.Status.SINGLETON : PeptideConsensusSpectrum.Status.NORMAL);

            return consensus;
        }

        private Peak calculatePrecursor(int charge, Peptide peptide, Collection<? extends PeakList> spectra) {

            double sumIntensity = 0.0;

            for (PeakList pl : spectra) {

                sumIntensity += pl.getPrecursor().getIntensity();
            }

            return new Peak(peptide.calculateMz(charge), sumIntensity, charge);
        }

        private Optional<RetentionTimeDiscrete> calculateRetentionTime(RetentionTimeList retentionTimes) {

            if (retentionTimes==null || retentionTimes.isEmpty()) return Optional.absent();

            double[] rts = new double[retentionTimes.size()];

            for (int i=0;i<retentionTimes.size();i++) {
                rts[i] = retentionTimes.get(i).getTime();
            }

            Arrays.sort(rts);

            return Optional.of(new RetentionTimeDiscrete(rts[retentionTimes.size()/2], TimeUnit.SECOND));
        }



        private static class MergePeakFilter<A extends PeakAnnotation> extends AbstractMergePeakFilter<A, PepLibPeakAnnotation> {

            public MergePeakFilter(double mzMaxDiff, double maxMzClusterWidth, IntensityMode intCombMeth, int totSpectrumCount) {

                super(mzMaxDiff, maxMzClusterWidth, intCombMeth, totSpectrumCount);
            }

            @Override
            protected PepLibPeakAnnotation createPeakAnnotation(int peakCount, double mzStd, double intensityStd) {

                return new PepLibPeakAnnotation(peakCount, mzStd, intensityStd);
            }
        }
    }
}
