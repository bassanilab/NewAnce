package newance.mzjava.mol;

import com.google.common.collect.Lists;
import newance.mzjava.ms.peaklist.PeakList;
import newance.mzjava.ms.spectrum.PepFragAnnotation;

import java.util.EnumSet;
import java.util.List;
import java.util.Set;

/**
 * Creates custom made PeptideFragmenter instances
 *
 * @author Markus Muller
 * @version 0.0
 */
public class PeptideFragmenterFactory {

    private PeptideFragmenterFactory() {

    }

    public static PeptideFragmenter newSequestPeptideFragmenter(PeakList.Precision precision) {

        Set<IonType> mainIonTypes = EnumSet.of(IonType.b, IonType.y);
        Set<IonType> bIonTypes = EnumSet.of(IonType.b);
        List<PeptidePeakGenerator<PepFragAnnotation>> peakGenerators = Lists.newArrayList(

                //Generate peaks of intensity 50 for b and y ions
                new BackbonePeakGenerator(mainIonTypes, 50),

                //Generate peaks of intensity 10 for water losses from b and y ions.
                //The peaks are only generated if the peptide fragment contains a S, T, D or E
                new PeptideNeutralLossPeakGenerator(Composition.parseComposition("H-2O-1"),
                        EnumSet.of(AminoAcid.S, AminoAcid.T, AminoAcid.D, AminoAcid.E), mainIonTypes, 10.0),

                //Generate peaks of intensity 10 for ammonium losses from b and y ions.
                //The peaks are only generated if the peptide fragment contains a Q, K, R or N
                new PeptideNeutralLossPeakGenerator(Composition.parseComposition("H-3N-1"),
                        EnumSet.of(AminoAcid.Q, AminoAcid.K, AminoAcid.R, AminoAcid.N), mainIonTypes, 10.0),

                //Generate that generates peaks of intensity 10 for CO loss from b ions
                new PeptideNeutralLossPeakGenerator(Composition.parseComposition("C-1O-1"),
                        EnumSet.allOf(AminoAcid.class), bIonTypes, 10.0)
        );

        return new PeptideFragmenter(peakGenerators, precision);
    }
}
