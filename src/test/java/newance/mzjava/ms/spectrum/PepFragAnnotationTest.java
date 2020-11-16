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
package newance.mzjava.ms.spectrum;

import newance.mzjava.mol.*;
import org.junit.Assert;
import org.junit.Test;

/**
 * @author Oliver Horlacher
 * @version 1.0
 */
public class PepFragAnnotationTest {

    @Test
    public void testEquals() throws Exception {

        PeptideFragment fragment = PeptideFragment.parse("CER", newance.mzjava.mol.FragmentType.INTERNAL);
        PepFragAnnotation annotation1 = new PepFragAnnotation.Builder(IonType.ay_internal, -2, fragment).addC13(2).setNeutralLoss(new NumericMass(0.001)).build();
        PepFragAnnotation annotation2 = new PepFragAnnotation.Builder(IonType.ay_internal, -2, fragment).addC13(2).setNeutralLoss(new NumericMass(0.001)).build();
        PepFragAnnotation annotation3 = new PepFragAnnotation.Builder(IonType.ay_internal, -1, fragment).addC13(2).setNeutralLoss(new NumericMass(0.001)).build();

        Assert.assertNotSame(annotation1, annotation2);

        Assert.assertTrue(annotation1.equals(annotation2));
        Assert.assertTrue(annotation2.equals(annotation1));

        Assert.assertFalse(annotation1.equals(annotation3));
        Assert.assertFalse(annotation3.equals(annotation1));
    }

    @Test
    public void testHashCode() throws Exception {

        PeptideFragment fragment = PeptideFragment.parse("CER", FragmentType.INTERNAL);
        PepFragAnnotation annotation1 = new PepFragAnnotation.Builder(IonType.ay_internal, -2, fragment).addC13(2).setNeutralLoss(new NumericMass(0.001)).build();
        PepFragAnnotation annotation2 = new PepFragAnnotation.Builder(IonType.ay_internal, -2, fragment).addC13(2).setNeutralLoss(new NumericMass(0.001)).build();
        PepFragAnnotation annotation3 = new PepFragAnnotation.Builder(IonType.ay_internal, -1, fragment).addC13(2).setNeutralLoss(new NumericMass(0.001)).build();


        int hash1 = annotation1.hashCode();
        int hash2 = annotation2.hashCode();
        int hash3 = annotation3.hashCode();

        Assert.assertEquals(hash1, hash2);
        Assert.assertTrue(hash1 != hash3);
    }

    @Test
    public void testCopy() throws Exception {

        PeptideFragment fragment = PeptideFragment.parse("CER", FragmentType.INTERNAL);
        PepFragAnnotation annotation = new PepFragAnnotation.Builder(IonType.ay_internal, -2, fragment).addC13(2).setNeutralLoss(new NumericMass(0.001)).build();

        PepFragAnnotation copyAnnotation = annotation.copy();
        PepFragAnnotation copyConstructedAnnotation = new PepFragAnnotation(annotation);

        Assert.assertNotSame(copyAnnotation, copyConstructedAnnotation);
        Assert.assertTrue(copyAnnotation.equals(copyConstructedAnnotation));
        Assert.assertTrue(copyConstructedAnnotation.equals(copyAnnotation));
    }

    @Test
    public void testHasNeutralLoss() throws Exception {

        PepFragAnnotation annotation = new PepFragAnnotation(IonType.c, 1, Peptide.parse("CER"));
        Assert.assertEquals(false, annotation.hasNeutralLoss());

        annotation = new PepFragAnnotation.Builder(IonType.z, 1, Peptide.parse("PEP")).setNeutralLoss(new NumericMass(0.0000000001)).build();
        Assert.assertEquals(true, annotation.hasNeutralLoss());
    }

    @Test
    public void testGetters() throws Exception {

        PepFragAnnotation annotation = new PepFragAnnotation(IonType.c, 1, Peptide.parse("CS(79.9)R"));

        Assert.assertEquals(IonType.c, annotation.getIonType());
        Assert.assertEquals(1, annotation.getCharge());
        Assert.assertEquals(PeptideFragment.parse("CS(79.9)R", FragmentType.FORWARD), annotation.getFragment());
        Assert.assertEquals(0, annotation.getIsotopeCount());
        Assert.assertEquals(Mass.ZERO, annotation.getNeutralLoss());
    }

    @Test
    public void testGetTheoreticalMZ() throws Exception {

        Assert.assertEquals(115.0865890597,
                new PepFragAnnotation(IonType.c, 1, Peptide.parse("P")).getTheoreticalMz(), 0.0000000001);

        Assert.assertEquals(116.08994389770001,
                new PepFragAnnotation.Builder(IonType.c, 1, Peptide.parse("P")).addC13(1).build().getTheoreticalMz(), 0.0000000001);

        Assert.assertEquals(116.0928658057,
                new PepFragAnnotation.Builder(IonType.c, 1, Peptide.parse("P")).addH2(1).build().getTheoreticalMz(), 0.0000000001);

        Atom n15 = PeriodicTable.getInstance().getAtom(AtomicSymbol.N, 15);
        Assert.assertEquals(116.08362395270001,
                new PepFragAnnotation.Builder(IonType.c, 1, Peptide.parse("P")).addIsotope(n15).build().getTheoreticalMz(), 0.0000000001);
    }

    @Test
    public void testBuilderCopyConstructor() throws Exception {

        PepFragAnnotation annotation = new PepFragAnnotation(IonType.b, 3, Peptide.parse("CERV"));

        PepFragAnnotation copy = new PepFragAnnotation.Builder(annotation).build();

        Assert.assertNotSame(annotation, copy);
        Assert.assertEquals(annotation, copy);
    }

    @Test
    public void testReuseBuilder() throws Exception {

        PepFragAnnotation.Builder builder = new PepFragAnnotation.Builder(IonType.y, 1, Peptide.parse("P"));

        PepFragAnnotation y1 = builder.build();
        Assert.assertEquals(IonType.y, y1.getIonType());
        Assert.assertEquals(1, y1.getCharge());
        Assert.assertEquals(0, y1.getIsotopeCount());
        Assert.assertEquals(PeptideFragment.parse("P", FragmentType.REVERSE), y1.getFragment());

        builder.setFragment(Peptide.parse("PE"));
        PepFragAnnotation y2 = builder.build();
        Assert.assertEquals(IonType.y, y2.getIonType());
        Assert.assertEquals(1, y2.getCharge());
        Assert.assertEquals(0, y2.getIsotopeCount());
        Assert.assertEquals(PeptideFragment.parse("PE", FragmentType.REVERSE), y2.getFragment());
    }

    @Test
    public void testAsString() throws Exception {

        Assert.assertEquals("y2^2", new PepFragAnnotation(IonType.y, 2, Peptide.parse("PE")).toSptxtString());
        Assert.assertEquals("y2i^2", new PepFragAnnotation.Builder(IonType.y, 2, Peptide.parse("PE")).addC13(1).build().toSptxtString());
        Assert.assertEquals("y2ii^2", new PepFragAnnotation.Builder(IonType.y, 2, Peptide.parse("PE")).addC13(2).build().toSptxtString());

        Assert.assertEquals("y2-18^2", new PepFragAnnotation.Builder(IonType.y, 2, Peptide.parse("PE")).setNeutralLoss(new NumericMass(-18)).build().toSptxtString());
    }
}
