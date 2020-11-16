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
package newance.mzjava.ms.peaklist;

import cern.colt.GenericPermuting;
import com.google.common.base.Function;
import com.google.common.base.Optional;
import newance.mzjava.ms.peaklist.peaktransformer.SqrtTransformer;
import org.junit.Assert;
import org.junit.Test;

import java.util.*;


/**
 * @author Oliver Horlacher
 * @version 1.0
 */
public abstract class BasePeakListTest {

    protected final double deltaMz;
    private final double deltaIntensity;
    protected final PeakList.Precision precision;

    private final boolean checkTotalIonCurrent;

    protected BasePeakListTest(PeakList.Precision precision) {

        switch (precision) {

            case DOUBLE:
            case DOUBLE_FLOAT:
            case DOUBLE_CONSTANT:
                this.deltaMz = 0.00000000001;
                break;
            case FLOAT:
            case FLOAT_CONSTANT:
                this.deltaMz = 0.0001;
                break;
            default:
                throw new IllegalArgumentException("Cannot set delta mz for " + precision);
        }

        switch (precision) {

            case DOUBLE:
                this.deltaIntensity = 0.000000000001;
                checkTotalIonCurrent = true;
                break;
            case DOUBLE_FLOAT:
            case FLOAT:
                this.deltaIntensity = 0.1;
                checkTotalIonCurrent = true;
                break;
            case DOUBLE_CONSTANT:
            case FLOAT_CONSTANT:
                this.deltaIntensity = 100000000;
                checkTotalIonCurrent = false;
                break;
            default:
                throw new IllegalArgumentException("Cannot set delta mz for " + precision);
        }

        this.precision = precision;
    }

    protected abstract <A extends PeakAnnotation> int getIntensityArrayLength(PeakList<A> peakList) throws NoSuchFieldException, IllegalAccessException;

    protected abstract <A extends PeakAnnotation> int getMzArrayLength(PeakList<A> peakList) throws NoSuchFieldException, IllegalAccessException;

    protected abstract void checkPeakList(double[] expectedMzs, double[] expectedIntensities, PeakList<? extends PeakAnnotation> peakList) throws NoSuchFieldException, IllegalAccessException;

    @Test
    protected abstract void testInsert() throws Exception;

    protected <A extends PeakAnnotation> void runTestInsert(PeakList<A> peakList) throws NoSuchFieldException, IllegalAccessException {

        int index;
        //Setup
        index = peakList.add(10, 1);
        Assert.assertEquals(0, index);
        index = peakList.add(15, 3);
        Assert.assertEquals(1, index);

        //Test insert in middle
        index = peakList.add(11, 2);
        Assert.assertEquals(1, index);
        checkPeakList(new double[]{10, 11, 15}, new double[]{1, 2, 3}, peakList);
        if (checkTotalIonCurrent) Assert.assertEquals(6, peakList.getTotalIonCurrent(), 0.0000001);

        //Test insert duplicate m/z
        index = peakList.add(11, 2);
        Assert.assertEquals(1, index);
        checkPeakList(new double[]{10, 11, 15}, new double[]{1, 4, 3}, peakList);
        if (checkTotalIonCurrent) Assert.assertEquals(8, peakList.getTotalIonCurrent(), 0.0000001);

        //Test insert at start
        index = peakList.add(9.99, 8);
        Assert.assertEquals(0, index);
        checkPeakList(new double[]{9.99, 10, 11, 15}, new double[]{8, 1, 4, 3}, peakList);
        if (checkTotalIonCurrent) Assert.assertEquals(16, peakList.getTotalIonCurrent(), 0.0000001);
    }

    @Test
    protected abstract void testAddSame() throws Exception;

    protected <A extends PeakAnnotation> void runTestAddSame(PeakList<A> peakList) {

        peakList.add(12.35823, 1);
        peakList.add(12.35823, 1);

        Assert.assertEquals(1, peakList.size());
        Assert.assertEquals(12.35823, peakList.getMz(0), deltaMz);
        Assert.assertEquals(2.0, peakList.getIntensity(0), deltaIntensity);
    }

    @Test
    public abstract void testBulkAdd() throws Exception;

    protected <A extends PeakAnnotation> void runTestBulkAdd(PeakList<A> peakList) throws NoSuchFieldException, IllegalAccessException {

        peakList = build(peakList, new double[]{87.6f, 98.65f, 123.54f, 169.54f});
        peakList.trimToSize();

        Assert.assertEquals(4, peakList.size());

        Assert.assertEquals(4, getMzArrayLength(peakList));
        Assert.assertEquals(4, getIntensityArrayLength(peakList));

        peakList.addSorted(
		        new double[]{201, 202, 203, 204, 205, 206},
		        new double[]{9, 8, 7, 6, 5, 4}, 3);

        checkPeakList(new double[]{87.6f, 98.65f, 123.54f, 169.54f, 201, 202, 203}, new double[]{1, 1, 1, 1, 9, 8, 7}, peakList);
    }

    <A extends PeakAnnotation> PeakList<A> build(PeakList<A> peakList, double[] mzs) {

        for (double mz : mzs) {

            peakList.add(mz, 1);
        }

        peakList.trimToSize();
        return peakList;
    }

    protected abstract <A extends PeakAnnotation> PeakList<A> newPeakList();

    protected <A extends PeakAnnotation> PeakList<A> build(PeakList<A> peakList, double[] mzs, double[] intensities) {

        Assert.assertEquals(mzs.length, intensities.length);
        peakList.ensureCapacity(mzs.length);

        for (int i = 0; i < mzs.length; i++) {
            peakList.add(mzs[i], intensities[i]);
        }

        peakList.trimToSize();
        return peakList;
    }

    <A extends PeakAnnotation> void build(PeakList<A> peakList, double[] mzs, double[] intensities, Map<Integer, A> annotationMap) {

        Assert.assertEquals(mzs.length, intensities.length);
        peakList.ensureCapacity(mzs.length);

        Set<Integer> annotSites = annotationMap.keySet();

        for (int i = 0; i < mzs.length; i++) {
            if (annotSites.contains(i)) {
                peakList.add(mzs[i], intensities[i], annotationMap.get(i));
            } else {
                peakList.add(mzs[i], intensities[i]);
            }
        }

        peakList.trimToSize();
    }

    @Test
    public abstract void testBulkAdd2() throws Exception;

    protected <A extends PeakAnnotation> void runTestBulkAdd2(PeakList<A> peakList) throws NoSuchFieldException, IllegalAccessException {

        build(peakList, new double[]{87.6f, 98.65f, 123.54f, 169.54f});
        peakList.trimToSize();

        Assert.assertEquals(4, peakList.size());

        Assert.assertEquals(4, getMzArrayLength(peakList));
        Assert.assertEquals(4, getIntensityArrayLength(peakList));

        peakList.addSorted(
		        new double[]{1, 2, 3, 4, 5, 6},
		        new double[]{9, 8, 7, 6, 5, 4}, 0);

        Assert.assertEquals(4, peakList.size());

        Assert.assertEquals(4, getMzArrayLength(peakList));
        Assert.assertEquals(4, getIntensityArrayLength(peakList));

        checkPeakList(new double[]{87.6f, 98.65f, 123.54f, 169.54f}, new double[]{1, 1, 1, 1}, peakList);
    }

    @Test
    public abstract void testBulkAdd3() throws Exception;

    protected <A extends PeakAnnotation> void runTestBulkAdd3(PeakList<A> peakList) throws NoSuchFieldException, IllegalAccessException {

        build(peakList, new double[]{87.6f, 98.65f, 123.54f, 169.54f});
        peakList.trimToSize();

        Assert.assertEquals(4, peakList.size());

        Assert.assertEquals(4, getMzArrayLength(peakList));
        Assert.assertEquals(4, getIntensityArrayLength(peakList));

        peakList.addSorted(
		        new double[]{201, 202, 203, 204, 205, 206},
		        new double[]{9, 8, 7, 6, 5, 4}, 6);

        Assert.assertEquals(10, peakList.size());

        Assert.assertEquals(10, getMzArrayLength(peakList));
        Assert.assertEquals(10, getIntensityArrayLength(peakList));

        checkPeakList(
                new double[]{87.6f, 98.65f, 123.54f, 169.54f, 201f, 202f, 203f, 204f, 205f, 206f},
                new double[]{1, 1, 1, 1, 9, 8, 7, 6, 5, 4},
                peakList);
    }

    @Test
    public abstract void testBulkAddException() throws Exception;

    protected <A extends PeakAnnotation> void runTestAddException(PeakList<A> peakList) throws NoSuchFieldException, IllegalAccessException {

        build(peakList, new double[]{87.6f, 98.65f, 123.54f, 169.54f});
        peakList.trimToSize();

        Assert.assertEquals(4, peakList.size());

        Assert.assertEquals(4, getMzArrayLength(peakList));
        Assert.assertEquals(4, getIntensityArrayLength(peakList));

        IndexOutOfBoundsException exception = null;
        try {
            peakList.addSorted(
		            new double[]{1, 2, 3, 4, 5, 6},
		            new double[]{9, 8, 7, 6, 5, 4}, -3);
        } catch (IndexOutOfBoundsException e) {

            exception = e;
        }
        Assert.assertNotNull(exception);

        Assert.assertEquals(4, peakList.size());

        Assert.assertEquals(4, getMzArrayLength(peakList));
        Assert.assertEquals(4, getIntensityArrayLength(peakList));

        checkPeakList(
                new double[]{87.6f, 98.65f, 123.54f, 169.54f},
                new double[]{1, 1, 1, 1},
                peakList);
    }

    @Test
    public abstract void testBulkAddException2() throws Exception;

    protected <A extends PeakAnnotation> void runTestBulkAddException2(PeakList<A> peakList) throws NoSuchFieldException, IllegalAccessException {

        build(peakList, new double[]{87.6f, 98.65f, 123.54f, 169.54f});
        peakList.trimToSize();

        Assert.assertEquals(4, peakList.size());

        Assert.assertEquals(4, getMzArrayLength(peakList));
        Assert.assertEquals(4, getIntensityArrayLength(peakList));

        IndexOutOfBoundsException exception = null;
        try {
            peakList.addSorted(
		            new double[]{1, 2, 3, 4, 5, 6},
		            new double[]{9, 8, 7, 6, 5, 4}, 25);
        } catch (IndexOutOfBoundsException e) {

            exception = e;
        }
        Assert.assertNotNull(exception);

        Assert.assertEquals(4, peakList.size());

        Assert.assertEquals(4, getMzArrayLength(peakList));
        Assert.assertEquals(4, getIntensityArrayLength(peakList));

        checkPeakList(
                new double[]{87.6f, 98.65f, 123.54f, 169.54f},
                new double[]{1, 1, 1, 1},
                peakList);
    }

    @Test
    public abstract void testBulkAddException3() throws Exception;

    protected <A extends PeakAnnotation> void runTestBulkAddException3(PeakList<A> peakList) throws NoSuchFieldException, IllegalAccessException {

        build(peakList, new double[]{87.6f, 98.65f, 123.54f, 169.54f});
        peakList.trimToSize();

        Assert.assertEquals(4, peakList.size());

        Assert.assertEquals(4, getMzArrayLength(peakList));
        Assert.assertEquals(4, getIntensityArrayLength(peakList));

        IndexOutOfBoundsException exception = null;
        try {
            peakList.addSorted(
		            new double[]{1, 2, 3, 4, 5, 6},
		            new double[]{9, 8, 7}, 6);
        } catch (IndexOutOfBoundsException e) {

            exception = e;
        }
        Assert.assertNotNull(exception);

        Assert.assertEquals(4, peakList.size());

        Assert.assertEquals(4, getMzArrayLength(peakList));
        Assert.assertEquals(4, getIntensityArrayLength(peakList));

        checkPeakList(
                new double[]{87.6f, 98.65f, 123.54f, 169.54f},
                new double[]{1, 1, 1, 1},
                peakList);
    }

    @Test
    public abstract void testBulkAddException4() throws Exception;

    protected <A extends PeakAnnotation> void runTestBulkAddException4(PeakList<A> peakList) throws NoSuchFieldException, IllegalAccessException {

        build(peakList, new double[]{87.6f, 98.65f, 123.54f, 169.54f});
        peakList.trimToSize();

        Assert.assertEquals(4, peakList.size());

        Assert.assertEquals(4, getMzArrayLength(peakList));
        Assert.assertEquals(4, getIntensityArrayLength(peakList));

        IndexOutOfBoundsException exception = null;
        try {
            peakList.addSorted(
		            new double[]{1, 2, 3},
		            new double[]{9, 8, 7, 6, 5, 4}, 6);

        } catch (IndexOutOfBoundsException e) {

            exception = e;
        }
        Assert.assertNotNull(exception);

        Assert.assertEquals(4, peakList.size());

        Assert.assertEquals(4, getMzArrayLength(peakList));
        Assert.assertEquals(4, getIntensityArrayLength(peakList));

        checkPeakList(
                new double[]{87.6f, 98.65f, 123.54f, 169.54f},
                new double[]{1, 1, 1, 1},
                peakList);
    }

    @Test
    public abstract void testClear() throws Exception;

    protected <A extends PeakAnnotation> void runTestClear(PeakList<A> peakList) throws NoSuchFieldException, IllegalAccessException {

        build(peakList, new double[]{87.6f, 98.65f, 123.54f, 169.54f});

        Assert.assertEquals(4, peakList.size());
        Assert.assertEquals(4, getMzArrayLength(peakList));
        Assert.assertEquals(4, getIntensityArrayLength(peakList));

        peakList.clear();
        Assert.assertEquals(0, getMzArrayLength(peakList));
        Assert.assertEquals(0, getIntensityArrayLength(peakList));
    }

    @Test
    public abstract void testCopyMzs();

    protected <A extends PeakAnnotation> void runTestCopyMzs(PeakList<A> peakList) {

        peakList.add(106.05, 1000);
        peakList.add(177.09, 10);
        peakList.add(233.06, 100);
        peakList.add(290.17, 10);
        peakList.add(389.16, 1000);
        peakList.trimToSize();

        double[] copy = peakList.getMzs(2, new double[3], 0, 3);

        Assert.assertArrayEquals(new double[]{233.06, 290.17, 389.16}, copy, deltaMz);
    }

    @Test
    public abstract void testEmptyPeakListGetIntensities() throws Exception;

    protected <A extends PeakAnnotation> void runTestEmptyPeakListGetIntensities(PeakList<A> peakList) {

        build(peakList, new double[0]);

        Assert.assertEquals(0, peakList.getIntensities(null).length);
    }

    @Test(expected = IndexOutOfBoundsException.class)
    public abstract void testEmptyPeakListGetMz() throws Exception;

    protected <A extends PeakAnnotation> void runTestEmptyPeakListGetMz(PeakList<A> peakList) {

        build(peakList, new double[0]);

        Assert.assertEquals(0.0, peakList.getMz(0), deltaMz);
    }

    @Test(expected = IndexOutOfBoundsException.class)
    public abstract void testEmptyPeakListGetIntensityAt() throws Exception;

    protected <A extends PeakAnnotation> void runTestEmptyPeakListGetIntensityAt(PeakList<A> peakList) {

        build(peakList, new double[0]);

        Assert.assertEquals(0.0, peakList.getIntensity(0), deltaIntensity);
    }

    @Test
    public abstract void testEmptyPeakListGetMzs() throws Exception;

    protected <A extends PeakAnnotation> void runTestEmptyPeakListGetMzs(PeakList<A> peakList) {

        build(peakList, new double[0]);

        Assert.assertEquals(0, peakList.getMzs(null).length);
    }

    @Test
    public abstract void testEnsureCapacity() throws NoSuchFieldException, IllegalAccessException;

    protected <A extends PeakAnnotation> void runTestEnsureCapacity(PeakList<A> peakList) throws NoSuchFieldException, IllegalAccessException {

        peakList.ensureCapacity(200);

        Assert.assertEquals(200, getMzArrayLength(peakList));
    }

    @Test
    public abstract void testEquals() throws Exception;

    protected <A extends PeakAnnotation> void runTestEquals(PeakList<A> peakList1, PeakList<A> peakList2) {

        build(peakList1, new double[]{4, 5, 6}, new double[]{1, 1, 1});
        build(peakList2, new double[]{4, 5, 6}, new double[]{2, 2, 2});

        Assert.assertEquals(false, peakList1.equals(peakList2));
    }

    @Test
    public abstract void testGetIntensities() throws Exception;

    protected <A extends PeakAnnotation> void runTestGetIntensities(PeakList<A> peakList) {

        build(peakList, new double[]{56, 87, 125, 986});

        double[] intensities = peakList.getIntensities(null);
        double[] expectedIntensities = new double[]{1, 1, 1, 1};
        Assert.assertArrayEquals(expectedIntensities, intensities, deltaIntensity);
    }

    @Test
    public abstract void testGetIntensities2() throws Exception;

    protected <A extends PeakAnnotation> void runTestGetIntensities2(PeakList<A> peakList) {

        build(peakList, new double[]{56, 87, 125, 986});

        double[] intensities = peakList.getIntensities(new double[4]);
        double[] expectedIntensities = new double[]{1, 1, 1, 1};
        Assert.assertArrayEquals(expectedIntensities, intensities, deltaIntensity);
    }

    @Test
    public abstract void testGetIntensities3() throws Exception;

    protected <A extends PeakAnnotation> void runTestGetIntensities3(PeakList<A> peakList) {

        build(peakList, new double[]{56, 87, 125, 986}, new double[]{100, 202, 354, 1});

        double[] dest = new double[]{-1, -1, -1, -1};
        double[] intensities = peakList.getIntensities(1, dest, 2, 2);
        double[] expectedIntensities = new double[]{-1, -1, 202, 354};
        Assert.assertArrayEquals(expectedIntensities, intensities, deltaIntensity);
    }

    @Test
    public abstract void testGetIntensities4() throws Exception;

    protected <A extends PeakAnnotation> void runTestGetIntensities4(PeakList<A> peakList) {

        build(peakList, new double[]{56, 87, 125, 986}, new double[]{100, 202, 354, 1});

        double[] intensities = peakList.getIntensities(1, null, 2, 2);
        double[] expectedIntensities = new double[]{0, 0, 202, 354};
        Assert.assertArrayEquals(expectedIntensities, intensities, deltaIntensity);
    }

    @Test
    public abstract void testGetIntensities5() throws Exception;

    protected <A extends PeakAnnotation> void runTestGetIntensities5(PeakList<A> peakList) {

        build(peakList, new double[]{56, 87, 125, 986}, new double[]{100, 202, 354, 1});

        double[] intensities = peakList.getIntensities(null, 3);
        double[] expectedIntensities = new double[]{0, 0, 0, 100, 202, 354, 1};
        Assert.assertArrayEquals(expectedIntensities, intensities, deltaIntensity);
    }

    @Test
    public abstract void testGetMzs1() throws Exception;

    protected <A extends PeakAnnotation> void runTestGetMzs1(PeakList<A> peakList) {

        build(peakList, new double[]{56, 87, 125, 986}, new double[]{100, 202, 354, 1});

        double[] dest = new double[]{-1, -1, -1, -1};
        double[] mzs = peakList.getMzs(1, dest, 2, 2);
        double[] expectedMzs = new double[]{-1, -1, 87, 125};
        Assert.assertArrayEquals(expectedMzs, mzs, deltaMz);
    }

    @Test
    public abstract void testGetMzs2() throws Exception;

    protected <A extends PeakAnnotation> void runTestGetMzs2(PeakList<A> peakList) {

        build(peakList, new double[]{56, 87, 125, 986}, new double[]{100, 202, 354, 1});

        double[] mzs = peakList.getMzs(1, null, 3, 2);
        double[] expectedMzs = new double[]{0, 0, 0, 87, 125};
        Assert.assertArrayEquals(expectedMzs, mzs, deltaMz);
    }

    @Test
    public abstract void testGetIntensitiesArr() throws Exception;

    protected <A extends PeakAnnotation> void runTestGetIntensitiesArr(PeakList<A> peakList) {

        build(peakList, new double[]{56, 87, 125, 986}, new double[]{1, 2, 3, 8});

        double[] actualMzs = peakList.getIntensities(new double[4]);
        Assert.assertArrayEquals(new double[]{1, 2, 3, 8}, actualMzs, deltaIntensity);
    }

    @Test
    public abstract void testGetIntensityAt() throws Exception;

    protected <A extends PeakAnnotation> void runTestGetIntensityAt(PeakList<A> peakList) {

        build(peakList, new double[]{56, 87, 125, 986}, new double[]{1, 2, 3, 4});

        double intensity = peakList.getIntensity(2);

        Assert.assertEquals(3.0, intensity, deltaIntensity);
    }

    @Test
    public abstract void testGetMostIntenseIndex();

    protected <A extends PeakAnnotation> void runTestGetMostIntenseIndex(PeakList<A> peakList) {

        if (!checkTotalIonCurrent) return;

        peakList.add(106.05, 1000);
        peakList.add(177.09, 10);
        peakList.add(233.06, 100);
        peakList.add(290.17, 10);
        peakList.add(389.16, 1000);
        peakList.trimToSize();

        Assert.assertEquals(2, peakList.getMostIntenseIndex(177.09, 300.17));
        Assert.assertEquals(2, peakList.getMostIntenseIndex(160.09, 300.17));
        Assert.assertEquals(4, peakList.getMostIntenseIndex(177.09, 389.16));
        Assert.assertEquals(0, peakList.getMostIntenseIndex(106.05, 210.39));
        Assert.assertEquals(-1, peakList.getMostIntenseIndex(390, 987.64));
    }

    @Test
    public abstract void testGetMostIntenseIndex2();

    protected <A extends PeakAnnotation> void runTestGetMostIntenseIndex2(PeakList<A> peakList) {

        peakList.trimToSize();

        Assert.assertEquals(-1, peakList.getMostIntenseIndex(177.09, 300.17));
    }

    @Test
    public abstract void testGetMostIntenseIndex3();

    protected <A extends PeakAnnotation> void runTestGetMostIntenseIndex3(PeakList<A> peakList) {

        peakList.add(106.05, 1000);
        peakList.add(177.09, 10);
        peakList.add(233.06, 100);
        peakList.add(290.17, 10);
        peakList.add(389.16, 1000);
        peakList.trimToSize();

        Assert.assertEquals(-1, peakList.getMostIntenseIndex(458.09, 500.17));
    }

    @Test
    public abstract void testGetMostIntenseIndex4();

    protected <A extends PeakAnnotation> void runTestGetMostIntenseIndex4(PeakList<A> peakList) {

        peakList.add(106.05, 1000);
        peakList.add(177.09, 10);
        peakList.add(233.06, 100);
        peakList.add(290.17, 10);
        peakList.add(389.16, 1000);
        peakList.trimToSize();

        Assert.assertEquals(2, peakList.getMostIntenseIndex(233.06, 300.0));
    }

    @Test
    public abstract void testGetMostIntenseIndex5();

    protected <A extends PeakAnnotation> void runTestGetMostIntenseIndex5(PeakList<A> peakList) {

        if (!checkTotalIonCurrent) return;

        peakList.add(106.05, 1000);
        peakList.add(177.09, 10);
        peakList.add(233.06, 100);
        peakList.add(290.17, 10);
        peakList.add(389.16, 1000);
        peakList.trimToSize();

        Assert.assertEquals(4, peakList.getMostIntenseIndex(233.06, 389.16));
    }

    @Test
    public void testGetMostIntenseIndex6() {

        if (!checkTotalIonCurrent) return;

        PeakList<PeakAnnotation> peakList = newPeakList();
        peakList.trimToSize();

        Assert.assertEquals(-1, peakList.getMostIntenseIndex(177.09, 300.17));
        Assert.assertEquals(-1, peakList.getMostIntenseIndex(160.09, 300.17));
        Assert.assertEquals(-1, peakList.getMostIntenseIndex(177.09, 389.16));
        Assert.assertEquals(-1, peakList.getMostIntenseIndex(106.05, 210.39));
        Assert.assertEquals(-1, peakList.getMostIntenseIndex(390, 987.64));
    }

    @Test
    public void testGetMostIntenseIndex7() {

        if (!checkTotalIonCurrent) return;

        PeakList<PeakAnnotation> peakList = newPeakList();
        peakList.add(100, 10);
        peakList.trimToSize();

        Assert.assertEquals(0, peakList.getMostIntenseIndex(90, 300.17));
    }

    @Test
    public void testGetMostIntenseIndex8() {

        if (!checkTotalIonCurrent) return;

        PeakList<PeakAnnotation> peakList = newPeakList();
        peakList.add(100, 10);
        peakList.trimToSize();

        Assert.assertEquals(-1, peakList.getMostIntenseIndex(100.1, 300.17));
    }

    @Test
    public abstract void testGetMzs3() throws Exception;

    protected <A extends PeakAnnotation>  void runTestGetMzs3(PeakList<A> peakList) {

        build(peakList, new double[]{56, 87, 125, 986}, new double[]{100, 202, 354, 1});

        double[] mzs = peakList.getMzs(null, 1);
        double[] expectedMzs = new double[]{0, 56, 87, 125, 986};
        Assert.assertArrayEquals(expectedMzs, mzs, deltaMz);
    }

    @Test
    public abstract void testHashCode() throws Exception;

    protected <A extends PeakAnnotation> void runTestHashCode(PeakList<A> peakList1, PeakList<A> peakList2) {

        build(peakList1, new double[]{4, 5, 6}, new double[]{1, 1, 1});
        build(peakList2, new double[]{4, 5, 6}, new double[]{2, 2, 2});

        Assert.assertTrue(peakList1.hashCode() != peakList2.hashCode());
    }

    @Test
    public abstract void testPrecision() throws Exception;

    protected <A extends PeakAnnotation> void runTestPrecision(PeakList<A> peakList, PeakList.Precision precision) {

        build(peakList, new double[]{87.6f, 98.65f, 123.54f, 169.54f});
        Assert.assertEquals(precision, peakList.getPrecision());
    }

    @Test
    public abstract void testSetIntensity();

    protected <A extends PeakAnnotation> void runTestSetIntensity(PeakList<A> peakList) {

        peakList.add(106.05, 1000);
        peakList.add(177.09, 10);
        peakList.add(233.06, 100);
        peakList.add(290.17, 10);
        peakList.add(389.16, 1000);
        peakList.trimToSize();

        if (checkTotalIonCurrent) Assert.assertEquals(2120.0, peakList.getTotalIonCurrent(), 0.00001);
        peakList.setIntensityAt(20, 1);

        if (checkTotalIonCurrent) Assert.assertEquals(2130.0, peakList.getTotalIonCurrent(), 0.00001);
    }

    @Test
    public abstract void testSetLoadFactor() throws NoSuchFieldException, IllegalAccessException;

    protected <A extends PeakAnnotation> void runTestSetLoadFactor(AbstractPeakList<A> peakList) throws NoSuchFieldException, IllegalAccessException {

        peakList.setLoadFactor(20);
        peakList.add(12.89, 785.3);

        Assert.assertEquals(20, getMzArrayLength(peakList));
        Assert.assertEquals(20, peakList.getLoadFactor());
    }

    @Test
    public abstract void testTotalIonCurrent();

    protected <A extends PeakAnnotation> void runTestTotalIonCurrent(PeakList<A> peakList) {

        build(peakList, new double[]{58.3, 87.6, 894.6}, new double[]{1, 2, 3});

        Assert.assertEquals(6.0, peakList.getTotalIonCurrent(), 0.00001);
    }

    @Test
    public abstract void testValueOf();

    protected void runTestValueOf(PeakList.Precision precision) {

        PeakList<PeakAnnotation> pl =
                AbstractPeakList
                        .valueOf("120.07913 Da, " +
                                "147.14932 Da, " +
                                "158.09187 Da, " +
                                "165.10173 Da", precision);

        Assert.assertEquals(4, pl.size());

        Assert.assertEquals(120.07913f, (float) pl.getMz(0), deltaMz);
        Assert.assertEquals(147.14932f, (float) pl.getMz(1), deltaMz);
        Assert.assertEquals(158.09187f, (float) pl.getMz(2), deltaMz);
        Assert.assertEquals(165.10173f, (float) pl.getMz(3), deltaMz);
    }

    @Test
    public abstract void testValueOfWithIntensities();

    protected void runTestValueOfWithIntensities(PeakList.Precision precision) {

        PeakList<PeakAnnotation> pl =
                AbstractPeakList
                        .valueOf("120.07913 Da (85002.6)," +
                                "147.14932 Da (20311.2)," +
                                "158.09187 Da (33321.5), " +
                                "165.10173 Da (487886.5), " +
                                "166.10524 Da (26983.5), " +
                                "175.11842 Da (428709.2), " +
                                "189.95845 Da (19962.7), " +
                                "191.11719 Da (33407.9), " +
                                "194.08276 Da (24121.1), " +
                                "214.08656 Da (55177.2), " +
                                "219.11382 Da (25720.4), " +
                                "231.11353 Da (44814.6), " +
                                "242.08058 Da (456825.2), " +
                                "243.08276 Da (46099.4), " +
                                "243.10849 Da (52023.8), " +
                                "248.13817 Da (37212.6)", precision);

        Assert.assertEquals(16, pl.size());

        Assert.assertEquals(120.07913, pl.getMz(0), deltaMz);
        Assert.assertEquals(85002.6, pl.getIntensity(0), deltaIntensity);
        Assert.assertEquals(147.14932, pl.getMz(1), deltaMz);
        Assert.assertEquals(20311.2, pl.getIntensity(1), deltaIntensity);
        Assert.assertEquals(158.09187, pl.getMz(2), deltaMz);
        Assert.assertEquals(33321.5, pl.getIntensity(2), deltaIntensity);
        Assert.assertEquals(165.10173, pl.getMz(3), deltaMz);
        Assert.assertEquals(487886.5, pl.getIntensity(3), deltaIntensity);
        Assert.assertEquals(166.10524, pl.getMz(4), deltaMz);
        Assert.assertEquals(26983.5, pl.getIntensity(4), deltaIntensity);
        Assert.assertEquals(175.11842, pl.getMz(5), deltaMz);
        Assert.assertEquals(428709.2, pl.getIntensity(5), deltaIntensity);
        Assert.assertEquals(189.95845, pl.getMz(6), deltaMz);
        Assert.assertEquals(19962.7, pl.getIntensity(6), deltaIntensity);
        Assert.assertEquals(191.11719, pl.getMz(7), deltaMz);
        Assert.assertEquals(33407.9, pl.getIntensity(7), deltaIntensity);
        Assert.assertEquals(194.08276, pl.getMz(8), deltaMz);
        Assert.assertEquals(24121.1, pl.getIntensity(8), deltaIntensity);
        Assert.assertEquals(214.08656, pl.getMz(9), deltaMz);
        Assert.assertEquals(55177.2, pl.getIntensity(9), deltaIntensity);
        Assert.assertEquals(219.11382, pl.getMz(10), deltaMz);
        Assert.assertEquals(25720.4, pl.getIntensity(10), deltaIntensity);
        Assert.assertEquals(231.11353, pl.getMz(11), deltaMz);
        Assert.assertEquals(44814.6, pl.getIntensity(11), deltaIntensity);
        Assert.assertEquals(242.08058, pl.getMz(12), deltaMz);
        Assert.assertEquals(456825.2, pl.getIntensity(12), deltaIntensity);
        Assert.assertEquals(243.08276, pl.getMz(13), deltaMz);
        Assert.assertEquals(46099.4, pl.getIntensity(13), deltaIntensity);
        Assert.assertEquals(243.10849, pl.getMz(14), deltaMz);
        Assert.assertEquals(52023.8, pl.getIntensity(14), deltaIntensity);
        Assert.assertEquals(248.13817, pl.getMz(15), deltaMz);
        Assert.assertEquals(37212.6, pl.getIntensity(15), deltaIntensity);
    }

    @Test
    public abstract void testAddPeakWithAnnotation() throws Exception;


    @Test
    public abstract void testGetMostIntenseIndexWholeSpectra() throws Exception;

    protected void runGetMostIntenseIndexWholeSpectra(PeakList peakList) {

        Assert.assertEquals(-1, peakList.getMostIntenseIndex());

        peakList.addSorted(new double[]{1, 2, 3, 4, 5, 6}, new double[]{1, 2, 5, 5, 4, 2});
        Assert.assertEquals(2, peakList.getMostIntenseIndex());
    }

    @Test
    public abstract void testMerge() throws Exception;

    protected <A extends PeakAnnotation> void runTestMerge(PeakList<A> peakList) {

        double[] masses =
                new double[]{995.465, 997.549, 1002.213, 1009.445, 1012.809, 1015.269, 1023.958, 1024.746, 1032.158, 1038.543,
                        1049.680, 1055.237, 1056.236, 1061.347, 1063.269, 1069.886, 1074.268, 1075.122, 1080.235, 1088.610,
                        1094.540, 1096.541, 1102.804, 1107.303, 1112.495, 1117.861, 1120.826, 1125.643, 1129.048, 1136.221,
                        1138.760, 1141.752, 1149.159, 1154.643, 1167.231, 1168.390, 1172.392, 1177.803, 1185.206, 1191.259,
                        1192.455, 1199.366, 1209.212, 1210.209, 1216.670};

        double[] masses2 =
                new double[]{995.833, 997.966, 1002.981, 1010.051, 1013.222, 1015.804, 1024.514, 1025.329, 1032.446, 1039.026,
                        1049.747, 1055.329, 1056.479, 1061.676, 1064.147, 1070.748, 1074.629, 1075.850, 1080.466, 1089.234,
                        1094.547, 1096.997, 1103.555, 1108.065, 1112.761, 1118.130, 1121.452, 1126.181, 1129.544, 1137.025,
                        1138.924, 1142.534, 1149.749, 1155.341, 1168.126, 1169.230, 1173.107, 1178.550, 1186.049, 1191.336,
                        1192.786, 1199.839, 1209.575, 1210.388, 1216.734};

        double[] intensities =
                new double[]{18.288, 10.924, 6.442, 17.300, 17.788, 18.947, 21.898, 11.378, 22.440, 10.302,
                        50.827, 17.054, 22.280, 11.892, 7.995, 10.341, 215.683, 88.473, 18.228, 13.208,
                        3.964, 5.049, 2.980, 1.780, 13.860, 5.808, 5.583, 2.364, 1.033, 13.762,
                        21.549, 5.541, 14.566, 19.965, 239.495, 40.558, 21.821, 9.397, 1820.518, 630.056,
                        23.832, 3.843, 6.159, 383.467, 58.726};

        peakList.addSorted(masses, intensities, masses.length);
        peakList.addSorted(masses2, intensities, masses.length);

        for (int i = 0; i < masses.length; i++) {

            Assert.assertEquals(masses[i], peakList.getMz(2 * i), deltaMz);
            Assert.assertEquals(masses2[i], peakList.getMz(2 * i + 1), deltaMz);
            Assert.assertEquals(intensities[i], peakList.getIntensity(2 * i), deltaIntensity);
            Assert.assertEquals(intensities[i], peakList.getIntensity(2 * i + 1), deltaIntensity);
        }
    }

    @Test
    public abstract void testMerge2();

    protected <A extends PeakAnnotation> void runTestMerge2(PeakList<A> consPeakList) {

        double[] masses =
                new double[] { 995.465, 997.549,1002.213,1009.445,1012.809,1015.269,1023.958,1024.746,1032.158,1038.543,
                        1049.680,1055.237,1056.236,1061.347,1063.269,1069.886,1074.268,1075.122,1080.235,1088.610,
                        1094.540,1096.541,1102.804,1107.303,1112.495,1117.861,1120.826,1125.643,1129.048,1136.221,
                        1138.760,1141.752,1149.159,1154.643,1167.231,1168.390,1172.392,1177.803,1185.206,1191.259,
                        1192.455,1199.366,1209.212,1210.209,1216.670};
        double[] intensities =
                new double[] {  18.288,  10.924,   6.442,  17.300,  17.788,  18.947,  21.898,  11.378,  22.440,  10.302,
                        50.827,  17.054,  22.280,  11.892,   7.995,  10.341, 215.683,  88.473,  18.228,  13.208,
                        3.964,   5.049,   2.980,   1.780,  13.860,   5.808,   5.583,   2.364,   1.033,  13.762,
                        21.549,   5.541,  14.566,  19.965, 239.495,  40.558,  21.821,   9.397,1820.518, 630.056,
                        23.832,   3.843,   6.159, 383.467,  58.726};

        int[] permutation = GenericPermuting.permutation(20, masses.length);
        double[] pMasses = new double[masses.length];
        double[] pIntensities = new double[masses.length];
        for (int i = 0; i < masses.length; i++) {
            int j = permutation[i];
            pMasses[i] = masses[j];
            pIntensities[i] = intensities[j];
        }


        double[] pSubSetMasses = new double[masses.length];
        double[] pSubSetIntensities = new double[masses.length];
        final Comparator<double[]> comparator = new Comparator<double[]>() {
            @Override
            public int compare(double[] o1, double[] o2) {

                return Double.compare(o1[0], o2[0]);
            }
        };

        double[][] data = new double[9][2];
        for (int i = 0; i < 5; i++) {
            System.arraycopy(pMasses, i * 9, pSubSetMasses, 0, 9);
            System.arraycopy(pIntensities, i * 9, pSubSetIntensities, 0, 9);

            for (int j = 0; j < 9; j++) {
                data[j][0] = pSubSetMasses[j];
                data[j][1] = pSubSetIntensities[j];
            }

            // Sort the IntDoublePairs
            Arrays.sort(data, 0, 9, comparator);

            for (int j = 0; j < 9; j++) {
                pSubSetMasses[j] = data[j][0];
                pSubSetIntensities[j] = data[j][1];
            }

            consPeakList.addSorted(pSubSetMasses, pSubSetIntensities, 9);
        }

        for (int i = 0; i < masses.length; i++) {
            Assert.assertEquals(masses[i], consPeakList.getMz(i), deltaMz);
            Assert.assertEquals(intensities[i], consPeakList.getIntensity(i), deltaIntensity);
        }
    }

    @Test
    public abstract void testMerge3();

    protected <A extends PeakAnnotation> void runTestMerge3(PeakList<A> peakList) {

        double[] masses =
                new double[] { 995.465, 997.549,1002.213,1009.445,1012.809,1015.269,1023.958,1024.746,1032.158,1038.543,
                        1049.680,1055.237,1056.236,1061.347,1063.269,1069.886,1074.268,1075.122,1080.235,1088.610,
                        1094.540,1096.541,1102.804,1107.303,1112.495,1117.861,1120.826,1125.643,1129.048,1136.221,
                        1138.760,1141.752,1149.159,1154.643,1167.231,1168.390,1172.392,1177.803,1185.206,1191.259,
                        1192.455,1199.366,1209.212,1210.209,1216.670};
        double[] intensities =
                new double[] {  18.288,  10.924,   6.442,  17.300,  17.788,  18.947,  21.898,  11.378,  22.440,  10.302,
                        50.827,  17.054,  22.280,  11.892,   7.995,  10.341, 215.683,  88.473,  18.228,  13.208,
                        3.964,   5.049,   2.980,   1.780,  13.860,   5.808,   5.583,   2.364,   1.033,  13.762,
                        21.549,   5.541,  14.566,  19.965, 239.495,  40.558,  21.821,   9.397,1820.518, 630.056,
                        23.832,   3.843,   6.159, 383.467,  58.726};

        for (int i = masses.length - 1; i >= 0; i--) {
            peakList.add(masses[i], intensities[i]);
        }
        Assert.assertArrayEquals(peakList.getMzs(null), masses, deltaMz);
        Assert.assertArrayEquals(peakList.getIntensities(null), intensities, deltaIntensity);
    }

    @Test
    public abstract void testMerge4();

    protected <A extends PeakAnnotation> void runTestMerge4(PeakList<A> consPeakList) {

        double[] masses =
                new double[] { 995.465, 997.549,1002.213,1009.445,1012.809,1015.269,1023.958,1024.746,1032.158,1038.543,
                        1049.680,1055.237,1056.236,1061.347,1063.269,1069.886,1074.268,1075.122,1080.235,1088.610,
                        1094.540,1096.541,1102.804,1107.303,1112.495,1117.861,1120.826,1125.643,1129.048,1136.221,
                        1138.760,1141.752,1149.159,1154.643,1167.231,1168.390,1172.392,1177.803,1185.206,1191.259,
                        1192.455,1199.366,1209.212,1210.209,1216.670};
        double[] intensities =
                new double[] {  18.288,  10.924,   6.442,  17.300,  17.788,  18.947,  21.898,  11.378,  22.440,  10.302,
                        50.827,  17.054,  22.280,  11.892,   7.995,  10.341, 215.683,  88.473,  18.228,  13.208,
                        3.964,   5.049,   2.980,   1.780,  13.860,   5.808,   5.583,   2.364,   1.033,  13.762,
                        21.549,   5.541,  14.566,  19.965, 239.495,  40.558,  21.821,   9.397,1820.518, 630.056,
                        23.832,   3.843,   6.159, 383.467,  58.726};

        int[] permutation = GenericPermuting.permutation(20, masses.length);
        double[] pMasses = new double[masses.length];
        double[] pIntensities = new double[masses.length];
        for (int i = 0; i < masses.length; i++) {
            int j = permutation[i];
            pMasses[i] = masses[j];
            pIntensities[i] = intensities[j];
        }

        for (int i = masses.length - 1; i >= 0; i--) {
            consPeakList.add(pMasses[i], pIntensities[i]);
        }
        Assert.assertArrayEquals(consPeakList.getMzs(null), masses, deltaMz);
        Assert.assertArrayEquals(consPeakList.getIntensities(null), intensities, deltaIntensity);
    }

    @Test
    public abstract void testMerge5() throws Exception;


    @Test
    public abstract void testDoInsert() throws Exception;

    protected <A extends PeakAnnotation> void runTestDoInsert(AbstractPeakList<A> peakList) {

        double vd = 58.0287398307;
        peakList.grow();
        peakList.doInsert(vd, 1);
        peakList.grow();
        peakList.doInsert(vd, 1);

        Assert.assertEquals(1, peakList.size());
    }

    @Test
    public abstract void testIndexOf() throws Exception;

    protected void runTestIndexOf(PeakList peakList) {

        peakList.add(12, 1);
        peakList.add(12.5, 1);
        peakList.add(124, 1);

        Assert.assertEquals(-1, peakList.indexOf(10));
        Assert.assertEquals(0, peakList.indexOf(12));
        Assert.assertEquals(-2, peakList.indexOf(12.4));
        Assert.assertEquals(1, peakList.indexOf(12.5));
        Assert.assertEquals(-3, peakList.indexOf(12.51));
        Assert.assertEquals(2, peakList.indexOf(124));
        Assert.assertEquals(-4, peakList.indexOf(124.1));
    }

    @Test
    public abstract void testGetClosestIndex();

    protected void runTestGetClosestIndex(PeakList peakList) {

        peakList.add(12, 1);
        peakList.add(12.5, 1);
        peakList.add(124, 1);

        Assert.assertEquals(0, peakList.getClosestIndex(10));
        Assert.assertEquals(0, peakList.getClosestIndex(12));
        Assert.assertEquals(0, peakList.getClosestIndex(12.249));

        Assert.assertEquals(0, peakList.getClosestIndex(12.25));
        Assert.assertEquals(1, peakList.getClosestIndex(12.5));
        Assert.assertEquals(1, peakList.getClosestIndex(20));

        Assert.assertEquals(2, peakList.getClosestIndex(124));
        Assert.assertEquals(2, peakList.getClosestIndex(1000));
    }

    @Test
    public void testIsEmpty() throws Exception {

        PeakList peakList = newPeakList();

        Assert.assertEquals(true, peakList.isEmpty());

        peakList.add(12, 3);

        Assert.assertEquals(false, peakList.isEmpty());
    }
    @Test
    public void testSetPrecursor() throws Exception {

        PeakList peakList = newPeakList();

        Peak peak = Peak.noIntensity(12);
        Assert.assertNotSame(peak, peakList.getPrecursor());

        peakList.setPrecursor(peak);
        Assert.assertSame(peak, peakList.getPrecursor());
    }
    @Test
    public void testApplyPeakProcessorChain() throws Exception {

        if (!checkTotalIonCurrent) return;

        PeakList<PeakAnnotation> peakList = newPeakList();

        peakList.addSorted(
                new double[]{1.2, 1.3, 1.4, 1.5, 1.6},
                new double[]{4, 9, 16, 25, 36}
        );

        peakList.apply(new PeakProcessorChain<>(new SqrtTransformer<>()));

        Assert.assertArrayEquals(new double[]{2, 3, 4, 5, 6}, peakList.getIntensities(new double[5]), 0.00001);
        Assert.assertEquals(2 + 3 + 4 + 5 + 6, peakList.getTotalIonCurrent(), deltaIntensity);
    }

    @Test (expected = UnsortedPeakListException.class)
    public void testAddUnsortedPeaksNotExpected() throws Exception {

        PeakList pl = newPeakList();

        pl.addSorted(new double[]{1, 20, 3}, new double[]{1, 1, 1});
    }
    @Test
    public void testIonCurrent() throws Exception {

        if (!checkTotalIonCurrent) return;

        PeakList<PeakAnnotation> peakList = newPeakList();

        peakList.add(100.0, 10.0);
        Assert.assertEquals(10, peakList.getTotalIonCurrent(), deltaIntensity);
        peakList.add(100.0, 10.0);
        Assert.assertEquals(20, peakList.getTotalIonCurrent(), deltaIntensity);
        peakList.add(200.0, 12.0);
        Assert.assertEquals(32, peakList.getTotalIonCurrent(), deltaIntensity);
        peakList.add(300.0, 13.0);
        Assert.assertEquals(45, peakList.getTotalIonCurrent(), deltaIntensity);
    }
    @Test
    public void testBasePeakMzAndIntensity() {

        if (!checkTotalIonCurrent) return;

        PeakList peakList = newPeakList();

        peakList.add(106.05, 1000);
        peakList.add(177.09, 10);
        peakList.add(233.06, 100);
        peakList.add(290.17, 10);
        peakList.add(389.16, 1000);
        peakList.trimToSize();

        Assert.assertEquals(0, peakList.getMostIntenseIndex());
    }

    @Test
    public void testGetBasePeak() throws Exception {

        PeakList peakList = newPeakList();

        if(peakList.getPrecision() == PeakList.Precision.FLOAT_CONSTANT || peakList.getPrecision() == PeakList.Precision.DOUBLE_CONSTANT)
            return;

        peakList.add(1, 1);
        peakList.add(2, 2);
        peakList.add(3, 3);
        peakList.add(4, 4);
        peakList.add(5, 10);
        peakList.add(6, 9);
        peakList.add(7, 8);
        peakList.add(8, 7);
        peakList.add(9, 6);
        peakList.add(10, 5);

        Assert.assertEquals(5.0, peakList.getBasePeakMz(), 0.00000001);
        Assert.assertEquals(10.0, peakList.getBasePeakIntensity(), 0.00000001);
    }
    @Test
    public void testCursor() throws Exception {

        PeakList<PeakAnnotation> peakList = newPeakList();

        peakList.add(1, 1);
        peakList.add(2, 2);
        peakList.add(3, 3);
        peakList.add(4, 4);
        peakList.add(5, 5);
        peakList.add(6, 6);

        PeakCursor<PeakAnnotation> cursor = peakList.cursor();

        Assert.assertEquals(3, cursor.getClosestIndex(3.9));
        Assert.assertEquals(2, cursor.getClosestIndex(3.5));

        cursor.movePast(4.3);
        Assert.assertEquals(5.0, cursor.currMz(), deltaMz);

        cursor.moveBefore(2.5);
        Assert.assertEquals(2.0, cursor.currMz(), deltaMz);
        cursor.moveBefore(2.9);
        Assert.assertEquals(2.0, cursor.currMz(), deltaMz);

        cursor.moveToClosest(5.5);
        Assert.assertEquals(5.0, cursor.currMz(), deltaMz);

        cursor.moveToClosest(5.51);
        Assert.assertEquals(6.0, cursor.currMz(), deltaMz);

        Assert.assertEquals(false, cursor.canPeek(1));
        Assert.assertEquals(false, cursor.canPeek(2));

        Assert.assertEquals(true, cursor.canPeek(-3));
        Assert.assertEquals(3.0, cursor.peekMz(-3), deltaMz);

        boolean previous = cursor.previous();
        Assert.assertEquals(true, previous);
        Assert.assertEquals(5.0, cursor.currMz(), deltaMz);

        cursor.resetCursor();
        boolean next = cursor.next();
        Assert.assertEquals(true, next);
        Assert.assertEquals(1.0, cursor.currMz(), deltaMz);
    }
}
