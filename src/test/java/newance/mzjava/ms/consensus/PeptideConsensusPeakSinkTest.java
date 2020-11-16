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

/**
 * @author Markus Muller
 * @version sqrt -1
 */
public class PeptideConsensusPeakSinkTest {

/*
    @Test
    public void testWithLosses() throws Exception {

        int minPeakCount = 2;
        double peakFraction = 0.5;

        MockPeakSink<PepLibPeakAnnotation> mockPeakSink = new MockPeakSink<>();

        PeptideConsensusPeakSink filter = new PeptideConsensusPeakSink(peakFraction, true);
        filter.setSink(mockPeakSink);

        PeptideSpectrum theoreticalSpectrum = makeTheoreticalSpectrum();
        filter.setParameters(theoreticalSpectrum, 10);

        filter.start(5);

        //Within tolerance of a2, keep
        filter.processPeak(205.2, 45, Collections.singletonList(new PepLibPeakAnnotation(1, 0, 0)));
        //Not within tolerance of theoretical peak, discard
        filter.processPeak(205.7,  5, Collections.singletonList(new PepLibPeakAnnotation(1, 0, 0)));
        //Not within tolerance of theoretical peak but mergedPeakCount >= 10*0.5, keep
        filter.processPeak(208.3, 65, Collections.singletonList(new PepLibPeakAnnotation(5, 0, 0)));
        //Within tolerance of b3-H2 and mergedPeakCount >= 10*0.5, keep
        filter.processPeak(371.1, 77, Collections.singletonList(new PepLibPeakAnnotation(5, 0, 0)));
        //Within tolerance of y3-H2O and mergedPeakCount < 10*0.5, keep
        filter.processPeak(385.1, 73, Collections.singletonList(new PepLibPeakAnnotation(1, 0, 0)));

        filter.end();

        Assert.assertEquals(1, mockPeakSink.getStartCalled());
        Assert.assertEquals(1, mockPeakSink.getEndCalled());

        double[] mzList = mockPeakSink.getMzList();
        double[] intensityList = mockPeakSink.getIntensityList();

        Assert.assertEquals(4, mzList.length);
        Assert.assertEquals(4, intensityList.length);

        Assert.assertEquals(205.2, mzList[0], 0);
        Assert.assertEquals(45, intensityList[0], 0);

        Assert.assertEquals(208.3, mzList[1], 0);
        Assert.assertEquals(65, intensityList[1], 0);

        Assert.assertEquals(371.1, mzList[2], 0);
        Assert.assertEquals(77, intensityList[2], 0);

        Assert.assertEquals(385.1, mzList[3], 0);
        Assert.assertEquals(73, intensityList[3], 0);
    }

    @Test
    public void testWithoutLosses() throws Exception {

        int absoluteMinPeakCount = 2;
        double peakFraction = 0.5;
        Tolerance tolerance = new AbsoluteTolerance(0.5);

        MockPeakSink<PepLibPeakAnnotation> mockPeakSink = new MockPeakSink<>();

        PeptideConsensusPeakSink filter = new PeptideConsensusPeakSink(tolerance, peakFraction, absoluteMinPeakCount, false);
        filter.setSink(mockPeakSink);

        PeptideSpectrum theoreticalSpectrum = makeTheoreticalSpectrum();
        filter.setParameters(theoreticalSpectrum, 10);

        filter.start(5);

        //Within tolerance of a2, keep
        filter.processPeak(205.2, 45, Collections.singletonList(new PepLibPeakAnnotation(1, 0, 0)));
        //Not within tolerance of theoretical peak, discard
        filter.processPeak(205.7,  5, Collections.singletonList(new PepLibPeakAnnotation(1, 0, 0)));
        //Not within tolerance of theoretical peak but mergedPeakCount >= 10*0.5, keep
        filter.processPeak(208.3, 65, Collections.singletonList(new PepLibPeakAnnotation(5, 0, 0)));
        //Within tolerance of b3-H2 and mergedPeakCount >= 10*0.5, keep
        filter.processPeak(371.1, 77, Collections.singletonList(new PepLibPeakAnnotation(5, 0, 0)));
        //Within tolerance of y3-H2O and mergedPeakCount < 10*0.5, discard
        filter.processPeak(385.1, 73, Collections.singletonList(new PepLibPeakAnnotation(1, 0, 0)));

        filter.end();

        Assert.assertEquals(1, mockPeakSink.getStartCalled());
        Assert.assertEquals(1, mockPeakSink.getEndCalled());

        double[] mzList = mockPeakSink.getMzList();
        double[] intensityList = mockPeakSink.getIntensityList();

        Assert.assertEquals(3, mzList.length);
        Assert.assertEquals(3, intensityList.length);

        Assert.assertEquals(205.2, mzList[0], 0);
        Assert.assertEquals(45, intensityList[0], 0);

        Assert.assertEquals(208.3, mzList[1], 0);
        Assert.assertEquals(65, intensityList[1], 0);

        Assert.assertEquals(371.1, mzList[2], 0);
        Assert.assertEquals(77, intensityList[2], 0);
    }

    @Test
    public void testSpectraPair() throws Exception {

        Tolerance tolerance = new AbsoluteTolerance(0.5);
        int absoluteMinPeakCount = 2;
        double peakFraction = 0.5;

        MockPeakSink<PepLibPeakAnnotation> mockPeakSink = new MockPeakSink<>();

        PeptideConsensusPeakSink filter = new PeptideConsensusPeakSink(tolerance, peakFraction, absoluteMinPeakCount, false);

        filter.setSink(mockPeakSink);

        PeptideSpectrum theoreticalSpectrum = makeTheoreticalSpectrum();
        filter.setParameters(theoreticalSpectrum, 2);

        filter.start(5);

        //Within tolerance of a2, keep
        filter.processPeak(205.2, 45, Collections.singletonList(new PepLibPeakAnnotation(1, 0, 0)));
        //Not within tolerance of theoretical peak, discard
        filter.processPeak(205.7,  5, Collections.singletonList(new PepLibPeakAnnotation(1, 0, 0)));
        //Not within tolerance of theoretical peak but mergedPeakCount >= 2*0.5, keep
        filter.processPeak(208.3, 65, Collections.singletonList(new PepLibPeakAnnotation(2, 0, 0)));
        //Within tolerance of b3-H2 and mergedPeakCount >= 2*0.5, keep
        filter.processPeak(371.1, 77, Collections.singletonList(new PepLibPeakAnnotation(2, 0, 0)));
        //Within tolerance of y3-H2O and mergedPeakCount < 2*0.5, discard
        filter.processPeak(385.1, 73, Collections.singletonList(new PepLibPeakAnnotation(1, 0, 0)));

        filter.end();

        Assert.assertEquals(1, mockPeakSink.getStartCalled());
        Assert.assertEquals(1, mockPeakSink.getEndCalled());

        double[] mzList = mockPeakSink.getMzList();
        double[] intensityList = mockPeakSink.getIntensityList();

        Assert.assertEquals(3, mzList.length);
        Assert.assertEquals(3, intensityList.length);

        Assert.assertEquals(205.2, mzList[0], 0);
        Assert.assertEquals(45, intensityList[0], 0);

        Assert.assertEquals(208.3, mzList[1], 0);
        Assert.assertEquals(65, intensityList[1], 0);

        Assert.assertEquals(371.1, mzList[2], 0);
        Assert.assertEquals(77, intensityList[2], 0);
    }

    private PeptideSpectrum makeTheoreticalSpectrum() {

        List<PeptidePeakGenerator<PepFragAnnotation>> peakGeneratorList = new ArrayList<>();
        peakGeneratorList.add(new BackbonePeakGenerator(EnumSet.of(IonType.a, IonType.b, IonType.y), 10));
        Mass waterLoss = Composition.parseComposition("H-2O-1");
        Set<AminoAcid> aa = EnumSet.of(AminoAcid.S, AminoAcid.T, AminoAcid.D, AminoAcid.E);
        Set<IonType> ions = EnumSet.of(IonType.b, IonType.y);
        peakGeneratorList.add(new PeptideNeutralLossPeakGenerator(waterLoss, aa, ions, 5));

        return new PeptideFragmenter(peakGeneratorList, PeakList.Precision.DOUBLE).fragment(Peptide.parse("CERV"), 2, new int[]{1});
    }
*/
}
