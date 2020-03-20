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
package newance.mzjava.mol;

import com.google.common.base.Preconditions;

import java.util.HashMap;
import java.util.Map;

/**
 * @author Oliver Horlacher
 * @version 1.0
 */
public class AAMassCalculator extends MassCalculator {

    private static AAMassCalculator instance;

    public static synchronized AAMassCalculator getInstance(){

        if(instance == null) instance = new AAMassCalculator();

        return instance;
    }

    private final Map<IonType, Composition> ionTypeCompositionMassMap;

    private AAMassCalculator() {

        // Source: Computational Methods for Mass Spectrometry Proteomics p.134
        ionTypeCompositionMassMap = new HashMap<>();
                                                                                                 //N-term    C-term
        ionTypeCompositionMassMap.put(IonType.a, Composition.parseComposition("C-1O-1H(+)"));    //[H]       [C-1O-1]
        ionTypeCompositionMassMap.put(IonType.b, Composition.parseComposition("H(+)"));         //[H]
        ionTypeCompositionMassMap.put(IonType.c, Composition.parseComposition("NH4(+)"));       //[H]       [NH3]

        ionTypeCompositionMassMap.put(IonType.x, Composition.parseComposition("CO2H(+)"));       //[CO]      [OH]
        ionTypeCompositionMassMap.put(IonType.y, Composition.parseComposition("H3O(+)"));       //[H2]      [OH]
        ionTypeCompositionMassMap.put(IonType.z, Composition.parseComposition("N-1O(+)"));   //[N-1H-1]  [OH]

        ionTypeCompositionMassMap.put(IonType.i, Composition.parseComposition("C-1O-1H(+)"));

        ionTypeCompositionMassMap.put(IonType.p, Composition.parseComposition("H2O"));       //[H]       [OH]
    }

    /**
     * Calculate the neutral molecular weight given a m/z, charge and mass of the charge carrying particle.
     *
     *
     * @param mz the m/z
     * @param charge the charge
     * @return the neutral molecular weight
     */
    public double calculateNeutralMolecularMass(double mz, int charge) {

        Preconditions.checkArgument(charge > 0, "Charge must be greater than 0");

        return (mz*charge) - PROTON_MASS*Math.abs(charge);
    }

    public double calculateMz(double molecularWeight, int charge, IonType ionType) {

        Preconditions.checkArgument(charge > 0, "Charge must be greater than 0");

        int chargeCount = ionType == IonType.p ? charge : charge - 1;

        return (molecularWeight + (chargeCount * PROTON_MASS)) /charge;
    }

    public final double getDeltaMass(IonType ionType) {

        if (ionTypeCompositionMassMap.containsKey(ionType)) {

            return ionTypeCompositionMassMap.get(ionType).getMolecularMass();
        } else {

            throw new IllegalStateException("Requesting an ion for which there is no data in this calculator.  The requested ion type was " + ionType);
        }
    }

    public final Composition getDeltaComposition(IonType ionType) {

        if (ionTypeCompositionMassMap.containsKey(ionType)) {

            return ionTypeCompositionMassMap.get(ionType);
        } else {

            throw new IllegalStateException("Requesting an ion for which there is no data in this calculator.  The requested ion type was " + ionType);
        }
    }
}
