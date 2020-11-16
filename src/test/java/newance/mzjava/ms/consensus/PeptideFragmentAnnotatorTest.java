/**
 * Copyright (c) 2010, SIB. All rights reserved.
 *
 * SIB (Swiss Institute of Bioinformatics) - http://www.isb-sib.ch Host -
 * http://mzjava.expasy.org
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer. Redistributions in binary
 * form must reproduce the above copyright notice, this list of conditions and
 * the following disclaimer in the documentation and/or other materials provided
 * with the distribution. Neither the name of the SIB/GENEBIO nor the names of
 * its contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL SIB/GENEBIO BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package newance.mzjava.ms.consensus;

import newance.mzjava.mol.*;
import newance.mzjava.mol.modification.AbsoluteTolerance;
import newance.mzjava.ms.peaklist.PeakCollectorSink;
import newance.mzjava.ms.peaklist.PeakList;
import newance.mzjava.ms.spectrum.PepFragAnnotation;
import newance.mzjava.ms.spectrum.PepLibPeakAnnotation;
import newance.mzjava.ms.spectrum.PeptideSpectrum;
import org.junit.Assert;
import org.junit.Test;

import java.util.*;

/**
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class PeptideFragmentAnnotatorTest {


    @Test
    public void test() throws Exception {

        List<PeptidePeakGenerator<PepFragAnnotation>> peakGeneratorList = new ArrayList<>();
        peakGeneratorList.add(new BackbonePeakGenerator(EnumSet.of(IonType.a, IonType.b, IonType.y), 10));
        Mass waterLoss = Composition.parseComposition("H-2O-1");
        Set<AminoAcid> aa = EnumSet.of(AminoAcid.S, AminoAcid.T, AminoAcid.D, AminoAcid.E);
        Set<IonType> ions = EnumSet.of(IonType.b, IonType.y);
        peakGeneratorList.add(new PeptideNeutralLossPeakGenerator(waterLoss, aa, ions, 5));

        Peptide peptide = Peptide.parse("CERV");
        PeptideSpectrum theoreticalSpectrum = new PeptideFragmenter(peakGeneratorList, PeakList.Precision.DOUBLE).fragment(peptide, 2, new int[]{1});

        Set<UUID> memberIds = new HashSet<>();
        memberIds.add(UUID.randomUUID());
        PeptideConsensusSpectrum consensus = new PeptideConsensusSpectrum(peptide,memberIds);

        PeptideFragmentAnnotator annotator = new PeptideFragmentAnnotator(new AbsoluteTolerance(0.5), theoreticalSpectrum);

        PeakCollectorSink<PepLibPeakAnnotation> sink = new PeakCollectorSink<>();
        sink.setPeakList(consensus);
        annotator.setSink(sink);

        annotator.start(20);

        annotator.processPeak(104.02, 10.0, Collections.singletonList(new PepLibPeakAnnotation(10, 0, 0)));	//b1
        annotator.processPeak(118.09, 10.0, Collections.singletonList(new PepLibPeakAnnotation(10, 0, 0)));	//y1
        annotator.processPeak(205.06, 10.0, Collections.singletonList(new PepLibPeakAnnotation(10, 0, 0)));	//a2
        annotator.processPeak(215.05, 5.0, Collections.singletonList(new PepLibPeakAnnotation(10, 0, 0)));	//b2-18.011

        annotator.processPeak(233.06, 10.0, Collections.singletonList(new PepLibPeakAnnotation(10, 0, 0)));	//b2
        annotator.processPeak(274.19, 10.0, Collections.singletonList(new PepLibPeakAnnotation(10, 0, 0)));	//y2
        annotator.processPeak(371.15, 5.0, Collections.singletonList(new PepLibPeakAnnotation(10, 0, 0)));	    //b3-18.011
        annotator.processPeak(385.22, 5.0, Collections.singletonList(new PepLibPeakAnnotation(10, 0, 0)));	    //y3-18.011
        annotator.processPeak(389.16, 10.0, Collections.singletonList(new PepLibPeakAnnotation(10, 0, 0)));	//b3
        annotator.processPeak(403.23, 10.0, Collections.singletonList(new PepLibPeakAnnotation(10, 0, 0)));	//y3

        annotator.end();

        Assert.assertEquals(10, consensus.size());

        checkAnnotation(IonType.b, 1, 1, Mass.ZERO, consensus.getAnnotations(0), consensus);
        checkAnnotation(IonType.y, 1, 1, Mass.ZERO, consensus.getAnnotations(1), consensus);
        checkAnnotation(IonType.a, 2, 1, Mass.ZERO, consensus.getAnnotations(2), consensus);
        checkAnnotation(IonType.b, 2, 1, waterLoss, consensus.getAnnotations(3), consensus);
        checkAnnotation(IonType.b, 2, 1, Mass.ZERO, consensus.getAnnotations(4), consensus);
        checkAnnotation(IonType.y, 2, 1, Mass.ZERO, consensus.getAnnotations(5), consensus);
        checkAnnotation(IonType.b, 3, 1, waterLoss, consensus.getAnnotations(6), consensus);
        checkAnnotation(IonType.y, 3, 1, waterLoss, consensus.getAnnotations(7), consensus);
        checkAnnotation(IonType.b, 3, 1, Mass.ZERO, consensus.getAnnotations(8), consensus);
        checkAnnotation(IonType.y, 3, 1, Mass.ZERO, consensus.getAnnotations(9), consensus);
    }

    @Test(timeout = 2000)
    public void test2() throws Exception {

        List<PeptidePeakGenerator<PepFragAnnotation>> peakGeneratorList = new ArrayList<>();
        peakGeneratorList.add(new BackbonePeakGenerator(EnumSet.of(IonType.a, IonType.b, IonType.y), 10));
        Mass waterLoss = Composition.parseComposition("H-2O-1");
        Set<AminoAcid> aa = EnumSet.of(AminoAcid.S, AminoAcid.T, AminoAcid.D, AminoAcid.E);
        Set<IonType> ions = EnumSet.of(IonType.b, IonType.y);
        peakGeneratorList.add(new PeptideNeutralLossPeakGenerator(waterLoss, aa, ions, 5));

        Peptide peptide = Peptide.parse("CERV");
        PeptideSpectrum theoreticalSpectrum = new PeptideFragmenter(peakGeneratorList, PeakList.Precision.DOUBLE).fragment(peptide, 2, new int[]{1});

        Set<UUID> memberIds = new HashSet<>();
        memberIds.add(UUID.randomUUID());
        memberIds.add(UUID.randomUUID());
        memberIds.add(UUID.randomUUID());
        memberIds.add(UUID.randomUUID());
        memberIds.add(UUID.randomUUID());
        PeptideConsensusSpectrum consensus = new PeptideConsensusSpectrum(peptide,memberIds);

        PeptideFragmentAnnotator annotator = new PeptideFragmentAnnotator(new AbsoluteTolerance(0.5),theoreticalSpectrum);

        PeakCollectorSink<PepLibPeakAnnotation> sink = new PeakCollectorSink<>();
        sink.setPeakList(consensus);
        annotator.setSink(sink);

        annotator.processPeak(103.50, 10.0, Collections.singletonList(new PepLibPeakAnnotation(10, 0, 0)));	//unknown
        annotator.processPeak(104.02, 10.0, Collections.singletonList(new PepLibPeakAnnotation(1, 0, 0)));	//b1
        annotator.processPeak(104.60, 10.0, Collections.singletonList(new PepLibPeakAnnotation(10, 0, 0)));	//unknown
        annotator.processPeak(371.15, 5.0, Collections.singletonList(new PepLibPeakAnnotation(10, 0, 0)));	//b3-18.011
        annotator.processPeak(371.65, 10.0, Collections.singletonList(new PepLibPeakAnnotation(10, 0, 0)));	//unknown
        annotator.processPeak(389.16, 10.0, Collections.singletonList(new PepLibPeakAnnotation(1, 0, 0)));	//b3
        annotator.processPeak(402.63, 10.0, Collections.singletonList(new PepLibPeakAnnotation(10, 0, 0)));	//unknown
        annotator.processPeak(403.23, 10.0, Collections.singletonList(new PepLibPeakAnnotation(10, 0, 0)));	//y3
        annotator.processPeak(457.23, 10.0, Collections.singletonList(new PepLibPeakAnnotation(1, 0, 0)));	//unknown

        annotator.end();

        Assert.assertEquals(9, consensus.size());

        Assert.assertEquals(false, consensus.getAnnotations(0).get(0).getOptFragmentAnnotation().isPresent());
        checkAnnotation(IonType.b, 1, 1, Mass.ZERO, consensus.getAnnotations(1), consensus);
        Assert.assertEquals(false, consensus.getAnnotations(2).get(0).getOptFragmentAnnotation().isPresent());
        checkAnnotation(IonType.b, 3, 1, waterLoss, consensus.getAnnotations(3), consensus);
        Assert.assertEquals(false, consensus.getAnnotations(4).get(0).getOptFragmentAnnotation().isPresent());
        checkAnnotation(IonType.b, 3, 1, Mass.ZERO, consensus.getAnnotations(5), consensus);
        Assert.assertEquals(false, consensus.getAnnotations(6).get(0).getOptFragmentAnnotation().isPresent());
        checkAnnotation(IonType.y, 3, 1, Mass.ZERO, consensus.getAnnotations(7), consensus);
    }

    @Test
    public void test3() throws Exception {

        List<PeptidePeakGenerator<PepFragAnnotation>> peakGeneratorList = new ArrayList<>();
        peakGeneratorList.add(new BackbonePeakGenerator(EnumSet.of(IonType.a, IonType.b, IonType.y), 10));
        Mass waterLoss = Composition.parseComposition("H-2O-1");
        Set<AminoAcid> aa = EnumSet.of(AminoAcid.S, AminoAcid.T, AminoAcid.D, AminoAcid.E);
        Set<IonType> ions = EnumSet.of(IonType.b, IonType.y);
        peakGeneratorList.add(new PeptideNeutralLossPeakGenerator(waterLoss, aa, ions, 5));

        Peptide peptide = Peptide.parse("CERV");
        PeptideSpectrum theoreticalSpectrum = new PeptideFragmenter(peakGeneratorList, PeakList.Precision.DOUBLE).fragment(peptide, 2, new int[]{1});

        Set<UUID> memberIds = new HashSet<>();
        memberIds.add(UUID.randomUUID());
        final PeptideConsensusSpectrum consensus = new PeptideConsensusSpectrum(peptide,memberIds);

        PeptideFragmentAnnotator annotator = new PeptideFragmentAnnotator(new AbsoluteTolerance(0.5),theoreticalSpectrum);

        PeakCollectorSink<PepLibPeakAnnotation> sink = new PeakCollectorSink<>();
        sink.setPeakList(consensus);
        annotator.setSink(sink);

        annotator.processPeak(103.50, 10.0, Collections.singletonList(new PepLibPeakAnnotation(10, 0, 0)));	//unknown
        annotator.processPeak(104.02, 10.0, Collections.singletonList(new PepLibPeakAnnotation(1, 0, 0)));	//b1
        annotator.processPeak(104.60, 10.0, Collections.singletonList(new PepLibPeakAnnotation(10, 0, 0)));	//unknown
        annotator.processPeak(371.15, 5.0, Collections.singletonList(new PepLibPeakAnnotation(10, 0, 0)));	//b3-18.011
        annotator.processPeak(371.65, 10.0, Collections.singletonList(new PepLibPeakAnnotation(10, 0, 0)));	//unknown
        annotator.processPeak(389.16, 10.0, Collections.singletonList(new PepLibPeakAnnotation(1, 0, 0)));	//b3
        annotator.processPeak(402.63, 10.0, Collections.singletonList(new PepLibPeakAnnotation(10, 0, 0)));	//unknown
        annotator.processPeak(403.23, 10.0, Collections.singletonList(new PepLibPeakAnnotation(10, 0, 0)));	//y3
        annotator.processPeak(457.23, 10.0, Collections.singletonList(new PepLibPeakAnnotation(1, 0, 0)));	//unknown

        annotator.end();

        Assert.assertEquals(9, consensus.size());

        Assert.assertEquals(false, consensus.getAnnotations(0).get(0).getOptFragmentAnnotation().isPresent());
        checkAnnotation(IonType.b, 1, 1, Mass.ZERO, consensus.getAnnotations(1), consensus);
        Assert.assertEquals(false, consensus.getAnnotations(2).get(0).getOptFragmentAnnotation().isPresent());
        checkAnnotation(IonType.b, 3, 1, waterLoss, consensus.getAnnotations(3), consensus);
        Assert.assertEquals(false, consensus.getAnnotations(4).get(0).getOptFragmentAnnotation().isPresent());
        checkAnnotation(IonType.b, 3, 1, Mass.ZERO, consensus.getAnnotations(5), consensus);
        Assert.assertEquals(false, consensus.getAnnotations(6).get(0).getOptFragmentAnnotation().isPresent());
        checkAnnotation(IonType.y, 3, 1, Mass.ZERO, consensus.getAnnotations(7), consensus);
    }

    private void checkAnnotation(IonType expectedIonType, int expectedResidueNumber, int expectedCharge, Mass expectedNeutralLoss, List<PepLibPeakAnnotation> actual, PeptideConsensusSpectrum peptideConsensusSpectrum) {

        Assert.assertEquals(1, actual.size());
        PepFragAnnotation actualAnnotation = actual.get(0).getOptFragmentAnnotation().get();

        Assert.assertEquals(expectedIonType, actualAnnotation.getIonType());
        Assert.assertEquals(peptideConsensusSpectrum.getPeptide().subSequence(expectedIonType, expectedResidueNumber), peptideConsensusSpectrum.getPeptide().subSequence(expectedIonType, expectedResidueNumber));
        Assert.assertEquals(expectedCharge, actualAnnotation.getCharge());
        Assert.assertEquals(expectedNeutralLoss, actualAnnotation.getNeutralLoss());
    }
}
