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
package newance.mzjava.mol.modification.unimod;

import newance.mzjava.mol.AminoAcid;
import newance.mzjava.mol.Composition;
import newance.mzjava.mol.modification.AbsoluteTolerance;
import newance.mzjava.mol.modification.Modification;
import newance.mzjava.mol.modification.NeutralLoss;
import newance.mzjava.mol.modification.Tolerance;
import org.junit.Assert;
import org.junit.Test;

import java.text.ParseException;
import java.util.*;

/**
 * @author Oliver Horlacher
 * @version 1.0
 */
public class UnimodManagerTest {

    @Test
    public void getModificationCount() throws Exception{

        Assert.assertEquals(971, UnimodManager.size());
    }

    @Test
    public void testGetModifications() throws ParseException {

        List<Modification> mods = UnimodManager.getModifications(29, new AbsoluteTolerance(0.5));

        Assert.assertEquals(mods.toString(), 5, mods.size());
        Assert.assertEquals(Modification.parseModification("Nitrosyl:H-1NO"), mods.get(0));
        Assert.assertEquals(Modification.parseModification("Val->Gln:H-1NO"), mods.get(1));
        Assert.assertEquals(Modification.parseModification("Ethyl+Deamidated:H3C2N-1O"), mods.get(2));
        Assert.assertEquals(Modification.parseModification("Val->Lys:CH3N"), mods.get(3));
        Assert.assertEquals(Modification.parseModification("Delta:H(5)C(2):C2H5"), mods.get(4));
    }

    @Test
    public void testPhospho() throws Exception {

        UnimodMod mod = UnimodManager.getModification("Phospho").get();

        Assert.assertEquals("Phospho", mod.getLabel());
        Assert.assertEquals("HO3P", mod.getMass().getFormula());

        List<NeutralLoss> neutralLosses = mod.getNeutralLosses();
        Assert.assertEquals(2, neutralLosses.size());

        Assert.assertEquals(0.0, neutralLosses.get(0).getMolecularMass(), 0.000001);
        Assert.assertEquals(new HashSet<String>(Arrays.asList("S", "T")), neutralLosses.get(0).getSites());

        Assert.assertEquals(97.97689508399999, neutralLosses.get(1).getMolecularMass(), 0.00000001);
        Assert.assertEquals(new HashSet<String>(Arrays.asList("S", "T")), neutralLosses.get(1).getSites());
    }

    @Test
    public void testIssue41() throws Exception {

        for(UnimodMod mod : UnimodManager.getInstance().getModificationList()) {

            Modification fromString = Modification.parseModification(mod.toString());

            Assert.assertEquals(mod, fromString);
        }
    }

    @Test(expected = UnsupportedOperationException.class)
    public void testGetModificationList() throws Exception {

        UnimodManager manager = UnimodManager.getInstance();

        List<UnimodMod> mods = manager.getModificationList();
        mods.clear();
    }

    @Test
    public void test_getModifications_withAA() throws Exception {
        double mass = Composition.parseComposition("CH2").getMolecularMass();
        Tolerance tol = new AbsoluteTolerance(0.01);
        Set<String> aas = new HashSet<>(Arrays.asList(new String[]{AminoAcid.S.getSymbol()}));

        List<Modification> modifs = UnimodManager.getModifications(mass,tol,aas);

        Assert.assertEquals(modifs.size(),2);
        Assert.assertEquals(modifs.get(0).toString(),"Methyl:CH2");
        Assert.assertEquals(modifs.get(1).toString(),"Ser->Thr:CH2");

        aas = new HashSet<>(Arrays.asList(new String[]{AminoAcid.K.getSymbol()}));

        modifs = UnimodManager.getModifications(mass,tol,aas);

        Assert.assertEquals(modifs.size(),1);
        Assert.assertEquals(modifs.get(0).toString(),"Methyl:CH2");

        mass = Composition.parseComposition("O").getMolecularMass();
        aas = new HashSet<>();
        for (AminoAcid aa : AminoAcid.values()) {
            aas.add(aa.getSymbol());
        }

        modifs = UnimodManager.getModifications(mass,tol,aas);

        Assert.assertEquals(modifs.size(),3);

        mass = Composition.parseComposition("C2H2O").getMolecularMass();
        aas = new HashSet<>(Arrays.asList(new String[]{"N-term"}));

        modifs = UnimodManager.getModifications(mass,tol,aas);

        Assert.assertEquals(modifs.size(),1);
    }

}