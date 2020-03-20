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

/**
 * @author Oliver Horlacher
 *
 * @version 1.0
 */
public abstract class MassCalculator {

    public static final double PROTON_MASS = PeriodicTable.H_MASS - PeriodicTable.ELECTRON_MASS;

    public abstract double getDeltaMass(IonType ionType);

    /**
     * Calculate the neutral molecular weight given a m/z, charge and mass of the charge carrying particle.
     *
     *
     * @param mz the m/z
     * @param charge the charge
     * @return the neutral molecular weight
     */
    public abstract double calculateNeutralMolecularMass(double mz, int charge);

    /**
     * Calculate the m/z for the given neutral molecular weight.
     *
     * @param molecularWeight the neutral molecular weight
     * @param charge the charge
     * @return the corresponding m/z
     */
    public abstract double calculateMz(double molecularWeight, int charge, IonType ionType);

    public static double calcIsotopeDelta(Composition composition) {

        PeriodicTable periodicTable = PeriodicTable.getInstance();
        double nonIsotopeMass = 0;
        for(Atom atom : composition.getAtoms()) {

            int count = composition.getCount(periodicTable.getAtom(atom.getSymbol(), atom.getMassNumber()));
            nonIsotopeMass += periodicTable.getAtom(atom.getSymbol()).getMass() * count;
        }

        return composition.getMolecularMass() - nonIsotopeMass;
    }

    /**
     * Calculates average mass of Atom
     *
     * @return the average mass of the given atom.
     */
    public static double calculateAverageMass(Atom atom) {

        PeriodicTable periodicTable = PeriodicTable.getInstance();

        if (!atom.isDefaultIsotope())
            return atom.getMass();

        double avgAtomMass = 0;

        for (Atom isotope : periodicTable.getAtoms(atom.getSymbol())) {

            avgAtomMass += isotope.getMass()*isotope.getAbundance();
        }

        return avgAtomMass;
    }

    /**
     * Calculates the sum of the average mass of each Atom in this
     * Composition.
     *
     * @return the average mass for this composition
     */
    public static double calculateAverageMass(Composition composition) {

        double avgMass = 0;

        for (Atom atom : composition.getAtoms()) {

            avgMass += calculateAverageMass(atom) * composition.getCount(atom);
        }

        avgMass -= composition.getCharge() * PeriodicTable.ELECTRON_MASS;

        return avgMass;
    }
}
