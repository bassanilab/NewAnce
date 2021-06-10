package newance.mzjava.mol;

import com.google.common.collect.Lists;
import newance.mzjava.mol.modification.AbsoluteTolerance;
import newance.mzjava.ms.consensus.PeptideFragmentCoverage;
import newance.mzjava.ms.peaklist.DoubleConstantPeakList;
import newance.mzjava.ms.peaklist.DoublePeakList;
import newance.mzjava.ms.peaklist.PeakList;
import newance.mzjava.ms.spectrum.PepFragAnnotation;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.Set;

/**
 * Created by markusmueller on 04.05.21.
 */
public class PeptideFragmentationTest {

    @Test
    public void test_psm_annot() {
// the fragmenter is responsible for generating a theoretical spectrum from the peptide
        Set<IonType> ionTypes = EnumSet.of(IonType.b, IonType.y, IonType.i, IonType.p);

        Mass waterLoss = Composition.parseComposition("H-2O-1");
        Mass ammoniumLoss = Composition.parseComposition("H-3N-1");

        List<PeptidePeakGenerator<PepFragAnnotation>> peakGenerators = Lists.newArrayList();

        peakGenerators.add(new BackbonePeakGenerator(ionTypes, 1));

        peakGenerators.add(new PeptideNeutralLossPeakGenerator(waterLoss,
                EnumSet.of(AminoAcid.S, AminoAcid.T, AminoAcid.D, AminoAcid.E), ionTypes, 1));

        peakGenerators.add(new PeptideNeutralLossPeakGenerator(ammoniumLoss,
                EnumSet.of(AminoAcid.Q, AminoAcid.K, AminoAcid.R, AminoAcid.N), ionTypes, 1));

        PeptideFragmenter fragmenter = new PeptideFragmenter(peakGenerators, PeakList.Precision.DOUBLE);

// the annotator needs to delegate to the fragmenter the fragmentation process and to
// the internal aligner the tolerance for aligning peaks
        PeptideFragmentAnnotator annotator = new PeptideFragmentAnnotator(fragmenter, new AbsoluteTolerance(0.1));

// the spectrum to annotate
        PeakList<PepFragAnnotation> spectrum = new DoubleConstantPeakList<>(1);

        double[] mzs = new double[] {234.1454, 321.1773, 322.1428, 330.1658, 353.1497, 359.1742, 365.1893, 365.219, 378.5135,
                402.6949, 402.718, 434.2426, 434.27, 434.2711, 439.2058, 456.2456, 471.2636, 491.2805, 508.2191, 553.2644,
                559.2601, 559.2625, 560.2307, 560.2513, 560.2559, 567.301, 575.7813, 576.2722, 576.3156, 592.3148, 592.3464,
                643.2743, 661.292, 674.309, 689.3341, 690.3297, 700.2772, 700.3068, 700.3408, 707.3326, 707.3625, 717.3297,
                718.2906, 718.3193, 735.3185, 735.3597, 780.1489, 786.3971, 804.3954, 804.4386};

        for (double mz : mzs) {

            spectrum.add(mz, 1);
        }

// the peptide as a source for theoretical fragments with annotations
        Peptide peptide = Peptide.parse("QVHPDTGISSK");

// start process: fragmenting peptide, aligning to the spectrum and reporting matched ions back to the spectrum
        annotator.annotate(spectrum, peptide, 3);


        for (int i=0;i<spectrum.size();i++) {
            String annotStr = "\t";
            if (spectrum.hasAnnotationsAt(i)) {
                List<PepFragAnnotation> annots = spectrum.getAnnotations(i);
                for (PepFragAnnotation annot : annots) annotStr += annot.toString()+",";
            }
            System.out.println(spectrum.getMz(i)+"\t"+spectrum.getIntensity(i)+annotStr);
        }
    }

    @Test
    public void calc_pept_cov_test() {
        Set<IonType> ionTypes = EnumSet.of(IonType.b, IonType.y, IonType.i, IonType.p);

        Mass waterLoss = Composition.parseComposition("H-2O-1");
        Mass ammoniumLoss = Composition.parseComposition("H-3N-1");

        List<PeptidePeakGenerator<PepFragAnnotation>> peakGenerators = Lists.newArrayList();

        peakGenerators.add(new BackbonePeakGenerator(ionTypes, 1));

        peakGenerators.add(new PeptideNeutralLossPeakGenerator(waterLoss,
                EnumSet.of(AminoAcid.S, AminoAcid.T, AminoAcid.D, AminoAcid.E), ionTypes, 1));

        peakGenerators.add(new PeptideNeutralLossPeakGenerator(ammoniumLoss,
                EnumSet.of(AminoAcid.Q, AminoAcid.K, AminoAcid.R, AminoAcid.N), ionTypes, 1));

        PeptideFragmenter fragmenter = new PeptideFragmenter(peakGenerators, PeakList.Precision.DOUBLE);

        PeptideFragmentAnnotator annotator = new PeptideFragmentAnnotator(fragmenter, new AbsoluteTolerance(0.1));

        PeakList<PepFragAnnotation> spectrum = new DoublePeakList<>(1);

        double[] mzs = new double[] {234.1454, 321.1773, 322.1428, 330.1658, 353.1497, 359.1742, 365.1893, 365.219, 378.5135,
                402.6949, 402.718, 434.2426, 434.27, 434.2711, 439.2058, 456.2456, 471.2636, 491.2805, 508.2191, 553.2644,
                559.2601, 559.2625, 560.2307, 560.2513, 560.2559, 567.301, 575.7813, 576.2722, 576.3156, 592.3148, 592.3464,
                643.2743, 661.292, 674.309, 689.3341, 690.3297, 700.2772, 700.3068, 700.3408, 707.3326, 707.3625, 717.3297,
                718.2906, 718.3193, 735.3185, 735.3597, 780.1489, 786.3971, 804.3954, 804.4386};


        for (double mz : mzs) {

            spectrum.add(mz, Math.random()*100.0);
        }

// the peptide as a source for theoretical fragments with annotations
        Peptide peptide = Peptide.parse("QVHPDTGISSK");

// start process: fragmenting peptide, aligning to the spectrum and reporting matched ions back to the spectrum
        annotator.annotate(spectrum, peptide, 3);

        PeptideFragmentCoverage peptideFragmentCoverage  = new PeptideFragmentCoverage(spectrum, peptide);

        System.out.println(peptideFragmentCoverage.getPeptideFragInfo());
        System.out.println(peptideFragmentCoverage.getPeptideFragInfoShort());
        System.out.println(peptideFragmentCoverage.getPeptideFragInfoLong());
        System.out.println(peptideFragmentCoverage.getSequenceCoverage());
        System.out.println(peptideFragmentCoverage.getSpectrumCoverage());
    }
}
