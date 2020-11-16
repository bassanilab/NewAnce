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

import com.google.common.base.Optional;
import com.google.common.collect.Sets;
import newance.mzjava.mol.*;
import newance.mzjava.mol.modification.AbsoluteTolerance;
import newance.mzjava.ms.peaklist.*;
import newance.mzjava.ms.peaklist.peakfilter.AbstractMergePeakFilter;
import newance.mzjava.ms.peaklist.peaktransformer.IdentityPeakProcessor;
import newance.mzjava.ms.spectrum.*;
import org.junit.Assert;
import org.junit.Test;

import java.lang.reflect.Field;
import java.util.*;

/**
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class PeptideConsensusSpectrumTest extends BasePeakListTest {

    public PeptideConsensusSpectrumTest() {

        super(PeakList.Precision.DOUBLE);
    }

    private double[] getMzField(DoublePeakList peakList) throws NoSuchFieldException, IllegalAccessException {

        Field field = peakList.getClass().getDeclaredField("mzList");
        field.setAccessible(true);
        return (double[]) field.get(peakList);
    }

    @Override
    protected void checkPeakList(double[] expectedMzs, double[] expectedIntensities, PeakList<? extends PeakAnnotation> spectrum) throws NoSuchFieldException, IllegalAccessException {

        DoublePeakList peakList = getPeakList(spectrum);

        peakList.trimToSize();
        double[] mzs = getMzField(peakList);
        Assert.assertEquals(expectedMzs.length, mzs.length);
        for (int i = 0; i < mzs.length; i++) {

            Assert.assertEquals(expectedMzs[i], mzs[i], 0.0001);
        }

        double[] intensities = getIntensityField(peakList);
        Assert.assertEquals(expectedIntensities.length, intensities.length);
        for (int i = 0; i < intensities.length; i++) {

            Assert.assertEquals(expectedIntensities[i], intensities[i], 0.0001);
        }
    }

    private double[] getIntensityField(DoublePeakList peakList) throws NoSuchFieldException, IllegalAccessException {

        Field field = peakList.getClass().getDeclaredField("intensityList");
        field.setAccessible(true);
        return (double[]) field.get(peakList);
    }

    private PeptideConsensusSpectrum newPeakList(int initialCapacity) {

        Set<UUID> memberIds = new HashSet<>();
        return new PeptideConsensusSpectrum(Peptide.parse("PEPTIDE"), initialCapacity, PeakList.Precision.DOUBLE,memberIds);
    }

    private PeptideConsensusSpectrum newPeakListRT() {

        Set<UUID> memberIds = new HashSet<>();
        memberIds.add(UUID.randomUUID());
        memberIds.add(UUID.randomUUID());
        Set<String> prots = new HashSet<>();
        prots.add("prot1");
        prots.add("prot2");
        return new PeptideConsensusSpectrum(Peptide.parse("PEPTIDE"), PeakList.Precision.DOUBLE,memberIds,prots,new RetentionTimeDiscrete(600.0, TimeUnit.SECOND));
    }

    @Override
    protected <A extends PeakAnnotation> int getIntensityArrayLength(PeakList<A> spectrum) throws NoSuchFieldException, IllegalAccessException {

        return getMzField(getPeakList(spectrum)).length;
    }

    @Override
    protected <A extends PeakAnnotation> int getMzArrayLength(PeakList<A> spectrum) throws NoSuchFieldException, IllegalAccessException {

        return getIntensityField(getPeakList(spectrum)).length;
    }

    private DoublePeakList getPeakList(PeakList spectrum) throws NoSuchFieldException, IllegalAccessException {

        Field field = Spectrum.class.getDeclaredField("peakList");
        field.setAccessible(true);
        return (DoublePeakList)field.get(spectrum);
    }

    protected PeptideConsensusSpectrum newPeakList() {

        return newPeakList(0);
    }

    @Override
    @Test
    public void testEquals() throws Exception {

        runTestEquals(newPeakList(), newPeakList());
    }

    @Override
    @Test
    public void testHashCode() throws Exception {

        runTestHashCode(newPeakList(), newPeakList());
    }

    @Override
    @Test
    public void testGetIntensitiesArr() throws Exception {

        runTestGetIntensitiesArr(newPeakList());
    }

    @Override
    @Test
    public void testGetIntensities() throws Exception {

        PeptideConsensusSpectrum peakList = newPeakList();

        runTestGetIntensities(peakList);
    }

    @Override
    @Test
    public void testGetIntensities2() throws Exception {

        runTestGetIntensities2(newPeakList());
    }

    @Override
    @Test
    public void testGetIntensities3() throws Exception {

        runTestGetIntensities3(newPeakList());
    }

    @Override
    @Test
    public void testGetIntensities4() throws Exception {

        runTestGetIntensities4(newPeakList());
    }

    @Override
    @Test
    public void testGetIntensities5() throws Exception {

        runTestGetIntensities5(newPeakList());
    }

    @Override
    @Test
    public void testGetMzs1() throws Exception {

        runTestGetMzs1(newPeakList());
    }

    @Override
    @Test
    public void testGetMzs2() throws Exception {

        runTestGetMzs2(newPeakList());
    }

    @Override
    @Test
    public void testGetMzs3() throws Exception {

        runTestGetMzs3(newPeakList());
    }

    @Override
    @Test
    public void testGetIntensityAt() throws Exception {

        PeptideConsensusSpectrum peakList = newPeakList();

        runTestGetIntensityAt(peakList);
    }

    @Test
    public void testAddSame() {

        runTestAddSame(newPeakList());
    }

    @Override
    @Test
    public void testEnsureCapacity() throws NoSuchFieldException, IllegalAccessException {

        runTestEnsureCapacity(newPeakList());
    }

    @Override
    @Test
    public void testSetLoadFactor() throws NoSuchFieldException, IllegalAccessException {

        //Only for PeakLists
    }

    @Override
    @Test
    public void testGetMostIntenseIndex() {

        runTestGetMostIntenseIndex(newPeakList(15));
    }

    @Override
    @Test
    public void testGetMostIntenseIndex2() {

        runTestGetMostIntenseIndex2(newPeakList(0));
    }

    @Override
    @Test
    public void testGetMostIntenseIndex3() {

        runTestGetMostIntenseIndex3(newPeakList(0));
    }

    @Override
    @Test
    public void testGetMostIntenseIndex4() {

        runTestGetMostIntenseIndex4(newPeakList(0));
    }

    @Override
    @Test
    public void testGetMostIntenseIndex5() {

        runTestGetMostIntenseIndex5(newPeakList(0));
    }

    @Override
    @Test
    public void testTotalIonCurrent() {

        runTestTotalIonCurrent(newPeakList());
    }

    @Override
    @Test
    public void testCopyMzs() {

        runTestCopyMzs(newPeakList(0));
    }

    @Override
    @Test
    public void testSetIntensity() {

        int initialCapacity = 0;
        runTestSetIntensity(newPeakList(initialCapacity));
    }

    @Override
    @Test
    public void testValueOf() {

        runTestValueOf(precision);
    }

    @Override
    @Test
    public void testValueOfWithIntensities() {

        runTestValueOfWithIntensities(precision);
    }

    @Override
    @Test
    public void testAddPeakWithAnnotation() throws Exception {

        // only for PeakList
    }

    @Override
    @Test
    public void testGetMostIntenseIndexWholeSpectra() throws Exception {

        runGetMostIntenseIndexWholeSpectra(newPeakList());
    }

    @Override
    @Test(expected = IndexOutOfBoundsException.class)
    public void testEmptyPeakListGetMz() throws Exception {

        PeptideConsensusSpectrum peakList = newPeakList();
        runTestEmptyPeakListGetMz(peakList);
    }

    @Override
    @Test
    public void testEmptyPeakListGetMzs() throws Exception {

        runTestEmptyPeakListGetMzs(newPeakList());
    }

    @Override
    @Test(expected = IndexOutOfBoundsException.class)
    public void testEmptyPeakListGetIntensityAt() throws Exception {

        runTestEmptyPeakListGetIntensityAt(newPeakList());
    }

    @Override
    @Test
    public void testEmptyPeakListGetIntensities() throws Exception {

        runTestEmptyPeakListGetIntensities(newPeakList());
    }

    @Override
    @Test
    public void testClear() throws Exception {

        PeptideConsensusSpectrum peakList = newPeakList();
        runTestClear(peakList);
    }

    @Override
    @Test
    public void testPrecision() throws Exception {

        PeptideConsensusSpectrum peakList = newPeakList();
        PeakList.Precision precision = this.precision;
        runTestPrecision(peakList, precision);
    }

    @Test
    public void testBulkAdd() throws Exception {

        runTestBulkAdd(newPeakList());
    }

    @Override
    @Test
    public void testBulkAdd2() throws Exception {

        runTestBulkAdd2(newPeakList());
    }

    @Override
    @Test
    public void testBulkAdd3() throws Exception {

        runTestBulkAdd3(newPeakList());
    }

    @Override
    @Test
    public void testBulkAddException() throws Exception {

        runTestAddException(newPeakList());
    }

    @Override
    @Test
    public void testBulkAddException2() throws Exception {

        runTestBulkAddException2(newPeakList());
    }

    @Override
    @Test
    public void testBulkAddException3() throws Exception {

        runTestBulkAddException3(newPeakList());
    }

    @Override
    @Test
    public void testBulkAddException4() throws Exception {

        runTestBulkAddException4(newPeakList());
    }

    @Test
    public void testInsert() throws Exception {

        runTestInsert(newPeakList());
    }

    @Test
    public void testMerge() throws Exception {

        runTestMerge(newPeakList());
    }

    @Test
    public void testMerge2() {

        runTestMerge2(newPeakList());
    }

    @Test
    public void testMerge3() {

        runTestMerge3(newPeakList());
    }

    @Test
    public void testMerge4() {

        runTestMerge4(newPeakList());
    }

    @Test
    public void testMerge5() throws Exception {

        //Only for PeakList
    }

    @Test
    public void testMerge6() throws Exception {

        //Only for PeakList
    }

    @Test
    public void testDoInsert() throws Exception {

        //Only for PeakList
    }

    @Override
    public void testIndexOf() throws Exception {

        runTestIndexOf(newPeakList());
    }

    @Override
    public void testGetClosestIndex() {

        runTestGetClosestIndex(newPeakList());
    }

    @Test
    public void testAddUnsortedPeaks() throws Exception {

        //Only for PeakList
    }

    @Test
    public void testGetPeptide() throws Exception {

        Peptide peptide = Peptide.parse("PEPTIDE");
        PeptideConsensusSpectrum spectrum = new PeptideConsensusSpectrum(peptide,new HashSet<UUID>());

        Assert.assertEquals(peptide, spectrum.getPeptide());
        Assert.assertSame(peptide, spectrum.getPeptide());
    }

    @Test
    public void testAddProteinAccessionNumber() throws Exception {

        PeptideConsensusSpectrum spectrum = new PeptideConsensusSpectrum(Peptide.parse("PEPTIDE"),new HashSet<UUID>());

        Assert.assertEquals(Collections.<String>emptySet(), spectrum.getProteinAccessionNumbers());

        spectrum.addProteinAccessionNumbers("ACC1");
        Assert.assertEquals(Collections.singleton("ACC1"), spectrum.getProteinAccessionNumbers());

        spectrum.addProteinAccessionNumbers(Collections.singleton("ACC2"));
        Assert.assertEquals(Sets.newLinkedHashSet(Arrays.asList("ACC1", "ACC2")), spectrum.getProteinAccessionNumbers());

        spectrum.clearProteinAccessionNumber();
        Assert.assertEquals(Collections.<String>emptySet(), spectrum.getProteinAccessionNumbers());
    }

    @Test(expected = UnsupportedOperationException.class)
    public void testUnmodifiableProteinAccessionNumber() throws Exception {

        PeptideConsensusSpectrum spectrum = new PeptideConsensusSpectrum(Peptide.parse("PEPTIDE"),new HashSet<UUID>());

        Assert.assertEquals(Collections.<String>emptySet(), spectrum.getProteinAccessionNumbers());

        spectrum.addProteinAccessionNumbers("ACC1");
        Assert.assertEquals(Collections.singleton("ACC1"), spectrum.getProteinAccessionNumbers());

        spectrum.addProteinAccessionNumbers(Collections.singleton("ACC2"));
        Assert.assertEquals(Sets.newLinkedHashSet(Arrays.asList("ACC1", "ACC2")), spectrum.getProteinAccessionNumbers());

        spectrum.getProteinAccessionNumbers().clear();
    }

    @Test
    public void testStatus() throws Exception {

        PeptideConsensusSpectrum spectrum = new PeptideConsensusSpectrum(Peptide.parse("PEPTIDE"),new HashSet<UUID>());

        Assert.assertEquals(PeptideConsensusSpectrum.Status.UNKNOWN, spectrum.getStatus());
        spectrum.setStatus(PeptideConsensusSpectrum.Status.INQUORATE_UNCONFIRMED);
        Assert.assertEquals(PeptideConsensusSpectrum.Status.INQUORATE_UNCONFIRMED, spectrum.getStatus());
    }

    @Test
    public void testGetSpectrumSource() throws Exception {

        PeptideConsensusSpectrum spectrum = new PeptideConsensusSpectrum(Peptide.parse("PEPTIDE"),new HashSet<UUID>());

        Assert.assertEquals(URIBuilder.UNDEFINED_URI, spectrum.getSpectrumSource());

        spectrum.setSpectrumSource(new URIBuilder("org.expasy", "test").build());
        Assert.assertEquals(new URIBuilder("org.expasy", "test").build(), spectrum.getSpectrumSource());
    }

    @Test
    public void testComment() throws Exception {

        PeptideConsensusSpectrum spectrum = new PeptideConsensusSpectrum(Peptide.parse("PEPTIDE"), PeakList.Precision.FLOAT,new HashSet<UUID>());

        Assert.assertEquals("", spectrum.getComment());

        spectrum.setComment("Some comment");
        Assert.assertEquals("Some comment", spectrum.getComment());
    }

    @Test
    public void testConsensusSpectrumCreation() {

        PepFragAnnotation pepFragAnnotation1 = new PepFragAnnotation(IonType.y,1, PeptideFragment.parse("PEPT", FragmentType.REVERSE));
        PepFragAnnotation pepFragAnnotation2 = new PepFragAnnotation(IonType.y,1, PeptideFragment.parse("PEPTI", FragmentType.REVERSE));
        PepFragAnnotation pepFragAnnotation3 = new PepFragAnnotation(IonType.y,1, PeptideFragment.parse("PEPTID", FragmentType.REVERSE));
        PepLibPeakAnnotation libPeakAnnotation1 = new PepLibPeakAnnotation(1,0.0,0.0, Optional.of(pepFragAnnotation1));
        PepLibPeakAnnotation libPeakAnnotation2 = new PepLibPeakAnnotation(1,0.0,0.0, Optional.of(pepFragAnnotation2));
        PepLibPeakAnnotation libPeakAnnotation3 = new PepLibPeakAnnotation(1,0.0,0.0, Optional.of(pepFragAnnotation3));

        PeptideConsensusSpectrum consensus = new PeptideConsensusSpectrum(Peptide.parse("PEPTIDE"),100, PeakList.Precision.DOUBLE,new HashSet<UUID>());
        consensus.add(pepFragAnnotation1.getTheoreticalMz(), 10.0, libPeakAnnotation1);
        consensus.add(pepFragAnnotation2.getTheoreticalMz(), 20.0, libPeakAnnotation2);
        consensus.add(pepFragAnnotation3.getTheoreticalMz(), 30.0, libPeakAnnotation3);

        PeptideConsensusSpectrum copy = consensus.copy(new IdentityPeakProcessor<PepLibPeakAnnotation>());

        Assert.assertTrue(consensus.equals(copy));

        Assert.assertEquals(copy.getAnnotationIndexes().length,3);
        Assert.assertEquals(copy.getAnnotations(0).get(0),libPeakAnnotation1);
        Assert.assertEquals(copy.getAnnotations(1).get(0),libPeakAnnotation2);
        Assert.assertEquals(copy.getAnnotations(2).get(0),libPeakAnnotation3);

        PepLibPeakAnnotation libPeakAnnotation4 = new PepLibPeakAnnotation(1,0.0,0.0, Optional.of(pepFragAnnotation1));

        PeptideConsensusSpectrum consensus2 = new PeptideConsensusSpectrum(Peptide.parse("PEPTIDE"),100, PeakList.Precision.DOUBLE,new HashSet<UUID>());

        consensus.add(pepFragAnnotation1.getTheoreticalMz(), 10.0, libPeakAnnotation1);
        consensus.add(pepFragAnnotation2.getTheoreticalMz(), 20.0, libPeakAnnotation2);
        consensus.add(pepFragAnnotation3.getTheoreticalMz(), 30.0, libPeakAnnotation4);

        Assert.assertFalse(consensus.equals(consensus2));
    }

    @Test
    public void test() {
        Peptide peptide = Peptide.parse("CERV");

        Set<MsnSpectrum> spectra = new HashSet<>();

        MsnSpectrum spectrum = new MsnSpectrum();
        spectrum.add(103.50, 10.0);    //unknown
        spectrum.add(104.02, 12.0);    //b1
        spectrum.add(371.15, 23.0);    //b3-18.011
        spectrum.add(371.64, 23.0);    //unknown
        spectrum.add(389.16, 23.0);    //b3
        spectrum.add(402.63, 23.0);    //unknown
        spectrum.add(403.23, 23.0);    //y3
        spectrum.add(457.23, 23.0);    //unknown
        spectrum.setPrecursor(new Peak(peptide.calculateMz(2) + 0.001, 10.0, 2));
        spectrum.setId(UUID.randomUUID());
        spectra.add(spectrum);

        MsnSpectrum spectrum2 = new MsnSpectrum();
        spectrum2.add(103.50, 11.0);    //unknown
        spectrum2.add(104.02, 12.0);    //b1
        spectrum2.add(371.65, 23.0);    //unknown
        spectrum2.add(389.16, 23.0);    //b3
        spectrum2.add(403.23, 23.0);    //y3
        spectrum2.add(457.23, 23.0);    //unknown
        spectrum2.setPrecursor(new Peak(peptide.calculateMz(2) + 0.001, 11.0, 2));
        spectrum2.setId(UUID.randomUUID());
        spectra.add(spectrum2);

        MsnSpectrum spectrum3 = new MsnSpectrum();
        spectrum3.add(103.50, 12.0);    //unknown
        spectrum3.add(104.02, 12.0);    //b1
        spectrum3.add(104.60, 23.0);    //unknown
        spectrum3.add(371.66, 23.0);    //unknown
        spectrum3.add(389.16, 23.0);    //b3
        spectrum3.add(402.63, 23.0);    //unknown
        spectrum3.add(403.23, 23.0);    //y3
        spectrum3.setPrecursor(new Peak(peptide.calculateMz(2) + 0.001, 12.0, 2));
        spectrum3.setId(UUID.randomUUID());
        spectra.add(spectrum3);

        List<PeptidePeakGenerator<PepFragAnnotation>> peakGeneratorList = new ArrayList<>();
        peakGeneratorList.add(new BackbonePeakGenerator(EnumSet.of(IonType.a, IonType.b, IonType.y), 10));
        Mass waterLoss = Composition.parseComposition("H-2O-1");
        Set<AminoAcid> aa = EnumSet.of(AminoAcid.S, AminoAcid.T, AminoAcid.D, AminoAcid.E);
        Set<IonType> ions = EnumSet.of(IonType.b, IonType.y);
        peakGeneratorList.add(new PeptideNeutralLossPeakGenerator(waterLoss, aa, ions, 5));

        PeptideConsensusSpectrum consensusSpectrum =
                PeptideConsensusSpectrum.builder().
                        peptide(peptide).
                        fragmenter(new PeptideFragmenter(peakGeneratorList, PeakList.Precision.DOUBLE)).
                        spectra(spectra).
                        setPeakFilterParams(0.2, 2).
                        fragMzTolerance(0.2).
                        intensityCombMethod(AbstractMergePeakFilter.IntensityMode.MEAN_ALL_INTENSITY).
                        build();

        Assert.assertEquals(consensusSpectrum.size(), 8);

        Assert.assertEquals(consensusSpectrum.getAnnotationIndexes().length,8);

        Assert.assertEquals(consensusSpectrum.getAnnotations(0).size(),1);
        Assert.assertEquals(consensusSpectrum.getAnnotations(0).get(0),new PepLibPeakAnnotation(3,0.0,Math.sqrt(2.0/3.0)));
        Assert.assertEquals(consensusSpectrum.getAnnotations(1).size(),1);

        PepFragAnnotation fragAnnotation = new PepFragAnnotation(IonType.b,1, Peptide.parse("C"));
        Assert.assertEquals(consensusSpectrum.getAnnotations(1).size(),1);
        Assert.assertEquals(consensusSpectrum.getAnnotations(1).get(0).getMergedPeakCount(),3);
        Assert.assertEquals(consensusSpectrum.getAnnotations(1).get(0).getMzStd(),0.0,0.00001);
        Assert.assertEquals(consensusSpectrum.getAnnotations(1).get(0).getIntensityStd(),0.0,0.00001);
        Assert.assertEquals(consensusSpectrum.getAnnotations(1).get(0).getOptFragmentAnnotation().get(),fragAnnotation);

        PepFragAnnotation.Builder builder = new PepFragAnnotation.Builder(IonType.b,1, Peptide.parse("CER"));
        fragAnnotation = builder.setNeutralLoss(Composition.parseComposition("H-2O-1")).build();
        Assert.assertEquals(consensusSpectrum.getAnnotations(2).size(),1);
        Assert.assertEquals(consensusSpectrum.getAnnotations(2).get(0).getMergedPeakCount(),1);
        Assert.assertEquals(consensusSpectrum.getAnnotations(2).get(0).getMzStd(),0.0,0.00001);
        Assert.assertEquals(consensusSpectrum.getAnnotations(2).get(0).getIntensityStd(),0.0,0.00001);
        Assert.assertEquals(consensusSpectrum.getAnnotations(2).get(0).getOptFragmentAnnotation().get(),fragAnnotation);

        Assert.assertEquals(consensusSpectrum.getAnnotations(3).size(),1);
        Assert.assertEquals(consensusSpectrum.getAnnotations(3).get(0).getMergedPeakCount(),3);
        double sd = Math.sqrt(2.0*0.0001/3);
        Assert.assertEquals(consensusSpectrum.getAnnotations(3).get(0).getMzStd(),sd,0.00001);
        Assert.assertEquals(consensusSpectrum.getAnnotations(3).get(0).getIntensityStd(),0.0,0.00001);
        Assert.assertFalse(consensusSpectrum.getAnnotations(3).get(0).getOptFragmentAnnotation().isPresent());

        fragAnnotation = new PepFragAnnotation(IonType.b,1, Peptide.parse("CER"));
        Assert.assertEquals(consensusSpectrum.getAnnotations(4).size(),1);
        Assert.assertEquals(consensusSpectrum.getAnnotations(4).get(0).getMergedPeakCount(),3);
        Assert.assertEquals(consensusSpectrum.getAnnotations(4).get(0).getMzStd(),0.0,0.00001);
        Assert.assertEquals(consensusSpectrum.getAnnotations(4).get(0).getIntensityStd(),0.0,0.00001);
        Assert.assertEquals(consensusSpectrum.getAnnotations(4).get(0).getOptFragmentAnnotation().get(),fragAnnotation);

        Assert.assertEquals(consensusSpectrum.getAnnotations(5).size(),1);
        Assert.assertEquals(consensusSpectrum.getAnnotations(5).get(0).getMergedPeakCount(),2);
        Assert.assertEquals(consensusSpectrum.getAnnotations(5).get(0).getMzStd(),0.0,0.00001);
        Assert.assertEquals(consensusSpectrum.getAnnotations(5).get(0).getIntensityStd(),0.0,0.00001);
        Assert.assertFalse(consensusSpectrum.getAnnotations(5).get(0).getOptFragmentAnnotation().isPresent());

        fragAnnotation = new PepFragAnnotation(IonType.y,1, Peptide.parse("ERV"));
        Assert.assertEquals(consensusSpectrum.getAnnotations(6).size(),1);
        Assert.assertEquals(consensusSpectrum.getAnnotations(6).get(0).getMergedPeakCount(),3);
        Assert.assertEquals(consensusSpectrum.getAnnotations(6).get(0).getMzStd(),0.0,0.00001);
        Assert.assertEquals(consensusSpectrum.getAnnotations(6).get(0).getIntensityStd(),0.0,0.00001);
        Assert.assertEquals(consensusSpectrum.getAnnotations(6).get(0).getOptFragmentAnnotation().get(),fragAnnotation);

        Assert.assertEquals(consensusSpectrum.getAnnotations(7).size(),1);
        Assert.assertEquals(consensusSpectrum.getAnnotations(7).get(0).getMergedPeakCount(),2);
        Assert.assertEquals(consensusSpectrum.getAnnotations(7).get(0).getMzStd(),0.0,0.00001);
        Assert.assertEquals(consensusSpectrum.getAnnotations(7).get(0).getIntensityStd(),0.0,0.00001);
        Assert.assertFalse(consensusSpectrum.getAnnotations(7).get(0).getOptFragmentAnnotation().isPresent());
    }

    @Test
    public void testBuilder2() {

        Peptide peptide = Peptide.parse("CERV");

        Set<MsnSpectrum> spectra = new HashSet<>();

        MsnSpectrum spectrum = new MsnSpectrum();
        spectrum.add(103.50, 10.0);    //unknown
        spectrum.add(104.02, 12.0);    //b1
        spectrum.add(371.15, 23.0);    //b3-18.011
        spectrum.add(371.64, 23.0);    //unknown
        spectrum.add(389.16, 23.0);    //b3
        spectrum.add(402.63, 23.0);    //unknown
        spectrum.add(403.23, 23.0);    //y3
        spectrum.add(457.23, 23.0);    //unknown
        spectrum.setPrecursor(new Peak(peptide.calculateMz(2) + 0.001, 10.0, 2));
        spectrum.setId(UUID.randomUUID());
        spectra.add(spectrum);

        MsnSpectrum spectrum2 = new MsnSpectrum();
        spectrum2.add(103.50, 11.0);    //unknown
        spectrum2.add(104.02, 12.0);    //b1
        spectrum2.add(371.65, 23.0);    //unknown
        spectrum2.add(389.16, 23.0);    //b3
        spectrum2.add(403.23, 23.0);    //y3
        spectrum2.add(457.23, 23.0);    //unknown
        spectrum2.setPrecursor(new Peak(peptide.calculateMz(2) + 0.001, 11.0, 2));
        spectrum2.setId(UUID.randomUUID());
        spectra.add(spectrum2);

        MsnSpectrum spectrum3 = new MsnSpectrum();
        spectrum3.add(103.50, 12.0);    //unknown
        spectrum3.add(104.02, 12.0);    //b1
        spectrum3.add(104.60, 23.0);    //unknown
        spectrum3.add(371.66, 23.0);    //unknown
        spectrum3.add(389.16, 23.0);    //b3
        spectrum3.add(402.63, 23.0);    //unknown
        spectrum3.add(403.23, 23.0);    //y3
        spectrum3.setPrecursor(new Peak(peptide.calculateMz(2) + 0.001, 12.0, 2));
        spectrum3.setId(UUID.randomUUID());
        spectra.add(spectrum3);

        List<PeptidePeakGenerator<PepFragAnnotation>> peakGeneratorList = new ArrayList<>();
        peakGeneratorList.add(new BackbonePeakGenerator(EnumSet.of(IonType.a, IonType.b, IonType.y), 10));
        Mass waterLoss = Composition.parseComposition("H-2O-1");
        Set<AminoAcid> aa = EnumSet.of(AminoAcid.S, AminoAcid.T, AminoAcid.D, AminoAcid.E);
        Set<IonType> ions = EnumSet.of(IonType.b, IonType.y);
        peakGeneratorList.add(new PeptideNeutralLossPeakGenerator(waterLoss, aa, ions, 5));

        PeptideConsensusSpectrum consensusSpectrum =
                PeptideConsensusSpectrum.builder(PeakList.Precision.DOUBLE, new URIBuilder("org.expasy.mzjava", "test").build())
                        .setConsensusParameters(0.2, 0.2, AbstractMergePeakFilter.IntensityMode.MEAN_ALL_INTENSITY)
                .setAnnotationParameters(new AbsoluteTolerance(0.2), new PeptideFragmenter(peakGeneratorList, PeakList.Precision.DOUBLE))
                .setFilterParameters(0.2, 2)
                .buildConsensus(2, peptide, spectra, Collections.<String>emptySet());

        Assert.assertEquals(consensusSpectrum.size(), 8);

        Assert.assertEquals(consensusSpectrum.getAnnotationIndexes().length,8);

        Assert.assertEquals(consensusSpectrum.getAnnotations(0).size(),1);
        Assert.assertEquals(consensusSpectrum.getAnnotations(0).get(0),new PepLibPeakAnnotation(3,0.0,Math.sqrt(2.0/3.0)));
        Assert.assertEquals(consensusSpectrum.getAnnotations(1).size(),1);

        PepFragAnnotation fragAnnotation = new PepFragAnnotation(IonType.b,1, Peptide.parse("C"));
        Assert.assertEquals(consensusSpectrum.getAnnotations(1).size(),1);
        Assert.assertEquals(consensusSpectrum.getAnnotations(1).get(0).getMergedPeakCount(),3);
        Assert.assertEquals(consensusSpectrum.getAnnotations(1).get(0).getMzStd(),0.0,0.00001);
        Assert.assertEquals(consensusSpectrum.getAnnotations(1).get(0).getIntensityStd(),0.0,0.00001);
        Assert.assertEquals(consensusSpectrum.getAnnotations(1).get(0).getOptFragmentAnnotation().get(),fragAnnotation);

        PepFragAnnotation.Builder builder = new PepFragAnnotation.Builder(IonType.b,1, Peptide.parse("CER"));
        fragAnnotation = builder.setNeutralLoss(Composition.parseComposition("H-2O-1")).build();
        Assert.assertEquals(consensusSpectrum.getAnnotations(2).size(),1);
        Assert.assertEquals(consensusSpectrum.getAnnotations(2).get(0).getMergedPeakCount(),1);
        Assert.assertEquals(consensusSpectrum.getAnnotations(2).get(0).getMzStd(),0.0,0.00001);
        Assert.assertEquals(consensusSpectrum.getAnnotations(2).get(0).getIntensityStd(),0.0,0.00001);
        Assert.assertEquals(consensusSpectrum.getAnnotations(2).get(0).getOptFragmentAnnotation().get(),fragAnnotation);

        Assert.assertEquals(consensusSpectrum.getAnnotations(3).size(),1);
        Assert.assertEquals(consensusSpectrum.getAnnotations(3).get(0).getMergedPeakCount(),3);
        double sd = Math.sqrt(2.0*0.0001/3);
        Assert.assertEquals(consensusSpectrum.getAnnotations(3).get(0).getMzStd(),sd,0.00001);
        Assert.assertEquals(consensusSpectrum.getAnnotations(3).get(0).getIntensityStd(),0.0,0.00001);
        Assert.assertFalse(consensusSpectrum.getAnnotations(3).get(0).getOptFragmentAnnotation().isPresent());

        fragAnnotation = new PepFragAnnotation(IonType.b,1, Peptide.parse("CER"));
        Assert.assertEquals(consensusSpectrum.getAnnotations(4).size(),1);
        Assert.assertEquals(consensusSpectrum.getAnnotations(4).get(0).getMergedPeakCount(),3);
        Assert.assertEquals(consensusSpectrum.getAnnotations(4).get(0).getMzStd(),0.0,0.00001);
        Assert.assertEquals(consensusSpectrum.getAnnotations(4).get(0).getIntensityStd(),0.0,0.00001);
        Assert.assertEquals(consensusSpectrum.getAnnotations(4).get(0).getOptFragmentAnnotation().get(),fragAnnotation);

        Assert.assertEquals(consensusSpectrum.getAnnotations(5).size(),1);
        Assert.assertEquals(consensusSpectrum.getAnnotations(5).get(0).getMergedPeakCount(),2);
        Assert.assertEquals(consensusSpectrum.getAnnotations(5).get(0).getMzStd(),0.0,0.00001);
        Assert.assertEquals(consensusSpectrum.getAnnotations(5).get(0).getIntensityStd(),0.0,0.00001);
        Assert.assertFalse(consensusSpectrum.getAnnotations(5).get(0).getOptFragmentAnnotation().isPresent());

        fragAnnotation = new PepFragAnnotation(IonType.y,1, Peptide.parse("ERV"));
        Assert.assertEquals(consensusSpectrum.getAnnotations(6).size(),1);
        Assert.assertEquals(consensusSpectrum.getAnnotations(6).get(0).getMergedPeakCount(),3);
        Assert.assertEquals(consensusSpectrum.getAnnotations(6).get(0).getMzStd(),0.0,0.00001);
        Assert.assertEquals(consensusSpectrum.getAnnotations(6).get(0).getIntensityStd(),0.0,0.00001);
        Assert.assertEquals(consensusSpectrum.getAnnotations(6).get(0).getOptFragmentAnnotation().get(),fragAnnotation);

        Assert.assertEquals(consensusSpectrum.getAnnotations(7).size(),1);
        Assert.assertEquals(consensusSpectrum.getAnnotations(7).get(0).getMergedPeakCount(),2);
        Assert.assertEquals(consensusSpectrum.getAnnotations(7).get(0).getMzStd(),0.0,0.00001);
        Assert.assertEquals(consensusSpectrum.getAnnotations(7).get(0).getIntensityStd(),0.0,0.00001);
        Assert.assertFalse(consensusSpectrum.getAnnotations(7).get(0).getOptFragmentAnnotation().isPresent());
    }

    @Test
    public void test_retentionTime() {
        PeptideConsensusSpectrum consensusSpectrum = newPeakList();
        Assert.assertFalse(consensusSpectrum.getRetentionTime().isPresent());

        PeptideConsensusSpectrum copy = consensusSpectrum.copy(new PeakProcessorChain<PepLibPeakAnnotation>());
        Assert.assertFalse(copy.getRetentionTime().isPresent());

        consensusSpectrum = newPeakListRT();

        Assert.assertTrue(consensusSpectrum.getRetentionTime().isPresent());
        Assert.assertEquals(600.0, consensusSpectrum.getRetentionTime().get().getTime(), 0.000001);

        copy = consensusSpectrum.copy(new PeakProcessorChain<PepLibPeakAnnotation>());

        Assert.assertTrue(copy.getRetentionTime().isPresent());
        Assert.assertEquals(600.0, copy.getRetentionTime().get().getTime(), 0.000001);

        copy = new PeptideConsensusSpectrum(consensusSpectrum,new PeakProcessorChain<PepLibPeakAnnotation>());

        Assert.assertTrue(copy.getRetentionTime().isPresent());
        Assert.assertEquals(600.0,copy.getRetentionTime().get().getTime(),0.000001);
    }
}
