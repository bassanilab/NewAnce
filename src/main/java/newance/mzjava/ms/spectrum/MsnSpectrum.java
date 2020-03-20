package newance.mzjava.ms.spectrum;

import com.google.common.base.Function;
import newance.mzjava.ms.peaklist.PeakAnnotation;
import newance.mzjava.ms.peaklist.PeakProcessor;
import newance.mzjava.ms.peaklist.PeakProcessorChain;
import newance.mzjava.ms.peaklist.peaktransformer.IdentityPeakProcessor;

import java.net.URI;
import java.util.Collection;

/**
 * @author Oliver Horlacher
 * @version 1.0
 */
public class MsnSpectrum extends Spectrum<PeakAnnotation> {

    private String fragMethod = "";
    private ScanNumber parentScanNumber = ScanNumber.SCAN_NUMBER_UNKNOWN;
    private final ScanNumberList scanNumbers;
    private final RetentionTimeList retentionTimes;
    private int spectrumIndex = -1;
    private String comment = "";
    private URI spectrumSource = URIBuilder.UNDEFINED_URI;

    /**
     * Default constructor
     */
    public MsnSpectrum() {

        this(0, Precision.DOUBLE);
    }

    /**
     * Construct a MsnSpectrum that stores the peak list at the given precision
     *
     * @param precision the precision that the peak list is stored at
     */
    public MsnSpectrum(Precision precision) {

        this(0, precision);
    }

    /**
     * Construct a MsnSpectrum that stores the peak list at the given precision and
     * has an initial capacity equal to <code>initialCapacity</code>
     *
     * @param initialCapacity the initial capacity of the peak list
     * @param precision       the precision that the peak list is stored at
     */
    public MsnSpectrum(int initialCapacity, Precision precision) {

        super(initialCapacity, precision);

        scanNumbers = new ScanNumberList();
        retentionTimes = new RetentionTimeList();
    }

    /**
     * Construct a MsnSpectrum that has a constant intensity. The precision needs to be either
     * PeakList.Precision.DOUBLE_CONSTANT or PeakList.Precision.FLOAT_CONSTANT. The peak lists initial
     * capacity is set to <code>initialCapacity</code>
     *
     * @param initialCapacity   the initial capacity of the peak list
     * @param constantIntensity the constant intensity if
     * @param precision         the precision that the peak list is stored at
     */
    public MsnSpectrum(int initialCapacity, double constantIntensity, Precision precision) {

        super(initialCapacity, constantIntensity, precision);

        scanNumbers = new ScanNumberList();
        retentionTimes = new RetentionTimeList();
    }

    /**
     * Copy constructor
     *
     * @param src the MsnSpectrum to copy
     */
    protected MsnSpectrum(MsnSpectrum src, PeakProcessor<PeakAnnotation,PeakAnnotation> peakProcessor) {

        super(src, peakProcessor);

        fragMethod = src.fragMethod;
        parentScanNumber = src.parentScanNumber;
        scanNumbers = src.scanNumbers.copy();
        retentionTimes = src.retentionTimes.copy();
        spectrumIndex = src.spectrumIndex;
        comment = src.comment;
        spectrumSource = src.spectrumSource;
    }

    /**
     * Copy constructor
     *
     * @param src the MsnSpectrum to copy
     */
    protected MsnSpectrum(MsnSpectrum src, PeakProcessor<PeakAnnotation,PeakAnnotation> peakProcessor,Function<Double,Double> rtCorrFct) {

        super(src, peakProcessor);

        fragMethod = src.fragMethod;
        parentScanNumber = src.parentScanNumber;
        scanNumbers = src.scanNumbers.copy();
        retentionTimes = new RetentionTimeList();
        for (RetentionTime rt :src.retentionTimes) {
            retentionTimes.add(rtCorrFct.apply(rt.getTime()),TimeUnit.SECOND);
        }
        spectrumIndex = src.spectrumIndex;
        comment = src.comment;
        spectrumSource = src.spectrumSource;
    }

    /**
     * Copy constructor
     *
     * @param src the MsnSpectrum to copy
     */
    protected MsnSpectrum(MsnSpectrum src, PeakProcessorChain<PeakAnnotation> peakProcessorChain) {

        super(src, peakProcessorChain);

        fragMethod = src.fragMethod;
        parentScanNumber = src.parentScanNumber;
        scanNumbers = src.scanNumbers.copy();
        retentionTimes = src.retentionTimes.copy();
        spectrumIndex = src.spectrumIndex;
        comment = src.comment;
        spectrumSource = src.spectrumSource;
    }

    /**
     * Copy constructor
     *
     * @param src the MsnSpectrum to copy
     */
    protected MsnSpectrum(MsnSpectrum src, PeakProcessorChain<PeakAnnotation> peakProcessorChain,Function<Double,Double> rtCorrFct) {

        super(src, peakProcessorChain);

        fragMethod = src.fragMethod;
        parentScanNumber = src.parentScanNumber;
        scanNumbers = src.scanNumbers.copy();
        retentionTimes = new RetentionTimeList();
        for (RetentionTime rt : src.retentionTimes) {
            retentionTimes.add(rtCorrFct.apply(rt.getTime()),TimeUnit.SECOND);
        }
        spectrumIndex = src.spectrumIndex;
        comment = src.comment;
        spectrumSource = src.spectrumSource;
    }

    @Override
    public MsnSpectrum copy(PeakProcessor<PeakAnnotation, PeakAnnotation> peakProcessor) {

        return new MsnSpectrum(this, peakProcessor);
    }

    @Override
    public MsnSpectrum copy(PeakProcessorChain<PeakAnnotation> peakProcessorChain) {

        return new MsnSpectrum(this, peakProcessorChain);
    }

    public MsnSpectrum copy(PeakProcessor<PeakAnnotation, PeakAnnotation> peakProcessor, Function<Double,Double> rtCorrFct) {

        return new MsnSpectrum(this, peakProcessor, rtCorrFct);
    }

    public MsnSpectrum copy(PeakProcessorChain<PeakAnnotation> peakProcessorChain,Function<Double,Double> rtCorrFct) {

        return new MsnSpectrum(this, peakProcessorChain, rtCorrFct);
    }

    public MsnSpectrum copy(Function<Double,Double> rtCorrFct) {

        return new MsnSpectrum(this, new IdentityPeakProcessor(), rtCorrFct);
    }

    /**
     * Return the fragmentation method
     *
     * @return the fragmentation method
     */
    public String getFragMethod() {

        return fragMethod;
    }

    /**
     * Set the fragmentation method
     *
     * @param fragMethod the fragmentation method
     */
    public void setFragMethod(String fragMethod) {

        this.fragMethod = fragMethod;
    }

    /**
     * Returns the scan number of the parent scan.
     *
     * @return the scan number of the parent scan
     */
    public ScanNumber getParentScanNumber() {

        return parentScanNumber;
    }

    /**
     * Set the parent scan number
     * @param parentScanNumber the new value for the parent scan number
     */
    public void setParentScanNumber(ScanNumber parentScanNumber) {

        this.parentScanNumber = parentScanNumber;
    }

    /**
     * Returns the scan number of this MsnSpectrum
     *
     * @return the scan number of this MsnSpectrum
     */
    public ScanNumberList getScanNumbers() {

        return scanNumbers;
    }

    /**
     * Add a all scan numbers from <code>scanNumbers</code> to this MsnSpectrum
     *
     * @param scanNumbers the scan numbers to add
     */
    public void addScanNumbers(Collection<? extends ScanNumber> scanNumbers) {

        this.scanNumbers.addAll(scanNumbers);
    }

    /**
     * Add <code>scanNumber</code> to this MsnSpectrum
     *
     * @param scanNumber the scan number to add
     */
    public void addScanNumber(ScanNumber scanNumber) {

        scanNumbers.add(scanNumber);
    }

    /**
     * Add <code>scanNumber</code> to this MsnSpectrum
     *
     * @param scanNumber the scan number to add
     */
    public void addScanNumber(int scanNumber) {

        scanNumbers.add(scanNumber);
    }

    /**
     * Return this MsnSpectrum's retention time list.
     *
     * @return this MsnSpectrum's retention time list
     */
    public RetentionTimeList getRetentionTimes() {

        return retentionTimes;
    }

    /**
     * Add all retention times from the <code>retentionTimes</code> to this
     * MsnSpectrum's retention times.
     *
     * @param retentionTimes the retention times to add
     */
    public void addRetentionTimes(RetentionTimeList retentionTimes) {

        this.retentionTimes.addAll(retentionTimes);
    }

    /**
     * Add the <code>retentionTime</code> to this MsnSpectrum's retention time list
     *
     * @param retentionTime the retention time to add
     */
    public void addRetentionTime(RetentionTime retentionTime) {

        retentionTimes.add(retentionTime);
    }

    /**
     * Return the index of this MsnSpectrum.
     *
     * @return the index of this MsnSpectrum
     */
    public int getSpectrumIndex() {

        return spectrumIndex;
    }

    /**
     * Set the spectrum index
     *
     * @param spectrumIndex the new spectrum index
     */
    public void setSpectrumIndex(int spectrumIndex) {

        this.spectrumIndex = spectrumIndex;
    }

    /**
     * Return the comment associated with this MsnSpectrum.
     *
     * @return the comment associated with this MsnSpectrum
     */
    public String getComment() {

        return comment;
    }

    /**
     * Set the comment
     *
     * @param comment the new value of the comment
     */
    public void setComment(String comment) {

        this.comment = comment;
    }

    /**
     * Return the URI that describes the source of this MsnSpectrum.
     *
     * @return the URI that describes the source of this MsnSpectrum
     */
    public URI getSpectrumSource() {

        return spectrumSource;
    }

    /**
     * Set the URI that describes the source of this MsnSpectrum
     *
     * @param spectrumSource the new spectrum source URI
     */
    public void setSpectrumSource(URI spectrumSource) {

        this.spectrumSource = spectrumSource;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof MsnSpectrum)) return false;
        if (!super.equals(o)) return false;

        MsnSpectrum that = (MsnSpectrum) o;

        if (spectrumIndex != that.spectrumIndex) return false;
        if (!comment.equals(that.comment)) return false;
        if (!fragMethod.equals(that.fragMethod)) return false;
        if (!parentScanNumber.equals(that.parentScanNumber)) return false;
        if (!retentionTimes.equals(that.retentionTimes)) return false;
        if (!scanNumbers.equals(that.scanNumbers)) return false;
        if (!spectrumSource.equals(that.spectrumSource)) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = super.hashCode();
        result = 31 * result + fragMethod.hashCode();
        result = 31 * result + parentScanNumber.hashCode();
        result = 31 * result + scanNumbers.hashCode();
        result = 31 * result + retentionTimes.hashCode();
        result = 31 * result + spectrumIndex;
        result = 31 * result + comment.hashCode();
        result = 31 * result + spectrumSource.hashCode();
        return result;
    }
}
