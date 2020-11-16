package newance.mzjava.mol;


import newance.mzjava.mol.modification.Tolerance;
import newance.mzjava.ms.peaklist.PeakList;
import newance.mzjava.ms.spectrasim.peakpairprocessor.DefaultPeakListAligner;
import newance.mzjava.ms.spectrasim.peakpairprocessor.PeakListAligner;
import newance.mzjava.ms.spectrasim.peakpairprocessor.PeakPairSink;
import newance.mzjava.ms.spectrum.PepFragAnnotation;
import newance.mzjava.ms.spectrum.PeptideSpectrum;

import java.util.List;

import static com.google.common.base.Preconditions.checkNotNull;

/**
 * This facade highly simplifies and hides the complexity of annotated a <code>peakList</code>
 * given a <code>peptide</code> as this process implies peptide fragmentation, aligning peaks and reporting
 * annotations.
 *
 * @author fnikitin
 * Date: 5/31/13
 */
public class PeptideFragmentAnnotator {

    private final PeptideFragmenter fragmenter;
    private final PeakListAligner<PepFragAnnotation, PepFragAnnotation> aligner;

    public PeptideFragmentAnnotator(PeptideFragmenter fragmenter, Tolerance tolerance) {

        this(fragmenter, new DefaultPeakListAligner<PepFragAnnotation, PepFragAnnotation>(tolerance));
    }

    public PeptideFragmentAnnotator(PeptideFragmenter fragmenter, PeakListAligner<PepFragAnnotation, PepFragAnnotation> aligner) {

        checkNotNull(fragmenter);
        checkNotNull(aligner);

        this.fragmenter = fragmenter;
        this.aligner = aligner;
        aligner.setSink(new AnnotatedSpectrumSink());
    }

    /**
     * Annotate peaklist from peptide at precursor charge
     *
     * @param peaklist the peaklist to annotate
     * @param peptide the peptide to fragment and to annotate peaklist with
     */
    public <S extends PeakList<PepFragAnnotation>> S annotate(S peaklist, Peptide peptide) {

        checkNotNull(peaklist.getPrecursor());

        return annotate(peaklist, peptide, peaklist.getPrecursor().getCharge());
    }

    /**
     * Annotate peaklist from peptide at specific charge
     *
     * @param peaklist the peaklist to annotate
     * @param peptide the peptide to fragment and to annotate peaklist with
     * @param precursorCharge the precursor charge (all fragments will be generated from +1 to precursorCharge)
     */
    public <S extends PeakList<PepFragAnnotation>> S annotate(S peaklist, Peptide peptide, int precursorCharge) {

        int[] charges = new int[precursorCharge];
        for (int i = 0; i < charges.length; i++) {

            charges[i] = i + 1;
        }

        return annotate(peaklist, peptide, precursorCharge, charges);
    }

    /**
     * Annotate peaklist from peptide at specific charge
     *
     * @param peaklist the peaklist to annotate
     * @param peptide the peptide to fragment and to annotate peaklist with
     * @param precursorCharge the precursor charge
     * @param fragmentCharges the fragment charges
     */
    public <S extends PeakList<PepFragAnnotation>> S annotate(S peaklist, Peptide peptide, int precursorCharge, int... fragmentCharges) {

        // fragment peptide and get theoretical annotated spectrum
        PeptideSpectrum peptideSpectrum = fragment(peptide, precursorCharge, fragmentCharges);

        if (peptideSpectrum.isEmpty())
            throw new IllegalStateException("cannot generate fragments from peptide precursor "+peptide.toString()+" at given charge "+precursorCharge);

        // make alignment and report matched peak annots to peaklist
        aligner.align(peaklist, peptideSpectrum);

        return peaklist;
    }

    protected PeptideSpectrum fragment(Peptide peptidePrecursor, int precursorCharge) {

        PeptideSpectrum peptideSpectrum = fragmenter.fragment(peptidePrecursor, precursorCharge);
        peptideSpectrum.setMsLevel(2);

        return peptideSpectrum;
    }

    protected PeptideSpectrum fragment(Peptide peptidePrecursor, int precursorCharge, int... fragmentCharges) {

        PeptideSpectrum peptideSpectrum = fragmenter.fragment(peptidePrecursor, precursorCharge, fragmentCharges);
        peptideSpectrum.setMsLevel(2);

        return peptideSpectrum;
    }

    protected void reportAnnotations(PeakList<PepFragAnnotation> peakList, int index, List<PepFragAnnotation> fragAnnotations) {

        for (PepFragAnnotation fragAnnotation : fragAnnotations) {

            peakList.addAnnotation(index, fragAnnotation);
        }
    }

    private class AnnotatedSpectrumSink implements PeakPairSink<PepFragAnnotation, PepFragAnnotation> {

        private PeakList<PepFragAnnotation> annotatedPeakList;
        private int currentExpPeakIndex;

        @Override
        public void begin(PeakList<PepFragAnnotation> xPeakList, PeakList<PepFragAnnotation> theoreticalPeakList) {

            annotatedPeakList = xPeakList;
            currentExpPeakIndex = -1;
        }

        @Override
        public void processPeakPair(double centroid, double xIntensity, double yIntensity,
                                    List<PepFragAnnotation> xAnnotations, List<PepFragAnnotation> yAnnotations) {

            if (xIntensity > 0) {

                currentExpPeakIndex++;
            }

            if (xIntensity > 0 && yIntensity > 0) {

                reportAnnotations(annotatedPeakList, currentExpPeakIndex, yAnnotations);
            }
        }

        @Override
        public void end() {

            currentExpPeakIndex = -1;
        }
    }
}
