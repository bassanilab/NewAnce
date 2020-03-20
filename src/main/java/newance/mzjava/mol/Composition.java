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

import com.google.common.collect.ImmutableSortedSet;
import gnu.trove.impl.unmodifiable.TUnmodifiableObjectIntMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.procedure.TObjectIntProcedure;

import java.io.Serializable;
import java.util.*;

import static com.google.common.base.Preconditions.checkNotNull;

/**
 * A collection of atoms and charge
 *
 * @author Oliver Horlacher
 * @author nikitin
 * @version 1.0
 */
public class Composition extends Mass {

    private double mass = 0;
    private final String formula;
    protected final int charge;
    protected final TObjectIntMap<Atom> atomCounterMap;

    /**
     * Construct a new composition that contains the <code>charge</code> and the atoms in the atomCounterMap.
     *
     * @param atomCounterMap the atoms to include in the composition
     * @param charge the charge of the composition
     */
    public Composition(TObjectIntMap<Atom> atomCounterMap, int charge) {

        this(Builder.makeFormula(atomCounterMap, charge),
                new TUnmodifiableObjectIntMap<Atom>(new TObjectIntHashMap<Atom>(atomCounterMap)), charge);
    }

    protected Composition(String formula, TObjectIntMap<Atom> atomCountMap, int charge) {

        this.charge = charge;
        this.atomCounterMap = atomCountMap;
        this.formula = formula;

        for(Atom atom : atomCountMap.keySet()) {

            mass += (atom.getMass() * atomCountMap.get(atom));
        }

        mass -= (this.charge * PeriodicTable.ELECTRON_MASS);
    }

    /**
     * Constructs a composition that is the sum of the <code>compositions</code>
     *
     * @param compositions the Compositions to sum
     */
    public Composition(Composition... compositions) {

        this.charge = sumCharges(compositions);

        atomCounterMap = new TObjectIntHashMap<Atom>();
        for(Composition composition : compositions) {

            for(Atom atom : composition.atomCounterMap.keySet()) {

                addAtoms(atom, composition.atomCounterMap.get(atom));
            }
        }

        if (compositions.length != 1) {
            formula = Builder.makeFormula(atomCounterMap, charge);
        } else {
            formula = compositions[0].formula;
        }

        mass -= (charge * PeriodicTable.ELECTRON_MASS);
    }

    private int sumCharges(Composition[] compositions) {

        int charge = 0;
        for (Composition composition : compositions) {

            charge += composition.getCharge();
        }

        return charge;
    }

    /**
     * Returns a set containing  all the Atoms in this composition.
     *
     * @return a set containing  all the Atoms in this composition.
     */
    public Set<Atom> getAtoms(){

        return ImmutableSortedSet.copyOf(atomCounterMap.keySet());
    }

    /**
     * Returns the monoisotopic molecular mass of this composition.
     *
     * <h2>Roundoff error and error propagation</h2>
     * At rare occasion, masses of different instances with exact same composition may slightly differ.
     * <p>
     * This is due to the fact that as our composition is a Map of unsorted atoms
     * the same atoms traversing order is not garantied.
     * <p>
     * As each instance mass is calculated by traversing atomic masses and that mass
     * are unexactly represented in the standard binary floating point, these roundoff errors<br>
     * can propagate through the calculation in non-intuitive ways and may lead to teany tiny different results.
     *
     * @return the monoisotopic molecular mass of this composition
     */

    @Override
    public double getMolecularMass() {

        return mass;
    }

    /**
     * Return the formula for this composition.
     *
     * @see AtomicCompositionParser
     *
     * @return the formula for this composition
     */
    public String getFormula() {

        return formula;
    }

    @Override
    public double getMassDefect() {

        final Counter counter = new Counter();

        atomCounterMap.forEachEntry(new TObjectIntProcedure<Atom>() {
            @Override
            public boolean execute(Atom atom, int n) {

                counter.increment(atom.getMassNumber() * n);
                return true;
            }
        });

        return getMolecularMass() - counter.getCount();
    }

    /**
     * Return the charge for this composition.
     *
     * @return the charge for this composition
     */
    public int getCharge() {

        return charge;
    }

    /**
     * Counts the number of times <code>atom</code> occurs in this Composition
     *
     * @param atom the atom
     * @return the number of times <code>atom</code> occurs in this Composition
     */
    public int getCount(Atom atom) {

        return atomCounterMap.containsKey(atom) ? atomCounterMap.get(atom) : 0;
    }

    /**
     * Returns the number of atoms in this composition
     *
     * @return the number of atoms in this composition
     */


    public int size() {

        int count = 0;

        for (int counter : atomCounterMap.values()) {

            count += counter;
        }

        return count;
    }

    /**
     * Returns true if this composition contains no atoms or charge, false otherwise.
     *
     * @return true if this composition contains no atoms or charge, false otherwise
     */
    public boolean isEmpty(){

        return atomCounterMap.isEmpty() && charge == 0;
    }

    private void addAtoms(Atom atom, int count) {

        if(count == 0) return;

        checkNotNull(atom);

        mass += atom.getMass() * count;

        atomCounterMap.adjustOrPutValue(atom, count, count);
    }

    @Override
    public boolean equals(Object o) {

        if (this == o) return true;
        if (!(o instanceof Composition)) return false;

        Composition that = (Composition) o;

        return charge == that.charge &&
                atomCounterMap.equals(that.atomCounterMap);

    }

    @Override
    public int hashCode() {

        int result = charge;
        result = 31 * result + atomCounterMap.hashCode();
        return result;
    }

    @Override
    public String toString() {

        return getFormula();
    }

    /**
     * Creates a Composition by parsing the given string.
     * For the format that is allowed see AtomicCompositionParser
     *
     * @see AtomicCompositionParser
     *
     * @param s the string to parse
     * @return the new Composition
     */
    public static Composition parseComposition(String s) {

        TObjectIntMap<Atom> atomCounterMap = new TObjectIntHashMap<Atom>();
        int charge = AtomicCompositionParser.getInstance().parse(s, atomCounterMap);

        return new Composition(s, atomCounterMap, charge);
    }

    /**
     * Class to build a Composition
     *
     * @author Oliver Horlacher
     * @version sqrt -1
     */
    public static final class Builder {

        private final PeriodicTable periodicTable = PeriodicTable.getInstance();
        private final TObjectIntHashMap<Atom> atomCounterMap = new TObjectIntHashMap<Atom>();
        private int charge = 0;

        private static final MzJavaIsotopeComparatorBySymbolAndA atomComparator = new MzJavaIsotopeComparatorBySymbolAndA();

        /**
         * Construct a new Builder
         */
        public Builder() {

        }

        /**
         * Construct a new Builder and add the <code>atom</code>
         *
         * @param atom the first atom to add to the builder
         */
        public Builder(AtomicSymbol atom) {

            add(atom, 1);
        }

        /**
         * Construct a new Builder and add the <code>atom</code> <code>count</code> times
         *
         * @param atom the first atom to add to the builder
         * @param count the number of times to add the atom
         */
        public Builder(AtomicSymbol atom, int count) {

            add(atom, count);
        }

        /**
         * Set the charge.
         *
         * @return this builder
         */
        public Builder charge(int charge) {

            this.charge = charge;
            return this;
        }

        /**
         * Add an atom for the specified AtomicSymbol
         *
         * @param atom the symbol of the atom to add
         * @return this builder
         */
        public Builder add(AtomicSymbol atom) {

            checkNotNull(atom);
            add(periodicTable.getAtom(atom), 1);

            return this;
        }

        /**
         * Add <code>count</code> atoms for the specified AtomicSymbol
         *
         * @param atom the symbol of the atom to add
         * @param count the number of times that the atom is to be added
         * @return this builder
         */
        public Builder add(AtomicSymbol atom, int count) {

            checkNotNull(atom);
            add(periodicTable.getAtom(atom), count);

            return this;
        }

        /**
         * Add an atom for the specified symbol and isotope
         *
         * @param symbol the symbol of the atom to add
         * @param isotope the isotope number
         * @return this builder
         */
        public Builder addIsotope(AtomicSymbol symbol, int isotope) {

            checkNotNull(symbol);
            add(periodicTable.getAtom(symbol, isotope), 1);

            return this;
        }

        /**
         * Add an atom <code>count<code/> times for the specified symbol and isotope
         *
         * @param symbol the symbol of the atom to add
         * @param isotope the isotope number
         * @param count the number of times to add the atom
         * @return this builder
         */
        public Builder addIsotope(AtomicSymbol symbol, int isotope, int count) {

            checkNotNull(symbol);
            add(periodicTable.getAtom(symbol, isotope), count);

            return this;
        }

        /**
         * Add the <code>atom</code> <code>count</code> times.
         *
         * @param atom the atom to add
         * @param count the number of atoms to add
         * @return this builder
         */
        public Builder add(Atom atom, int count) {

            if(count == 0) return this;
            checkNotNull(atom);
            atomCounterMap.adjustOrPutValue(atom, count, count);

            return this;
        }

        /**
         * Build the composition
         *
         * @return the new Composition
         */
        public Composition build(){

            return new Composition(atomCounterMap, charge);
        }

        /**
         * Return the atomic composition string defined in AtomicCompositionParser
         *
         * @see AtomicCompositionParser
         *
         * @return the formula
         */
        @Override
        public String toString() {

            return makeFormula(atomCounterMap, charge);
        }

        /**
         * Make the atomic formula given an <code>atomCounterMap</code> and <code>charge</code>
         *
         * @param atomCounterMap map containing the counts for the atoms
         * @param charge the charge
         * @return the formula
         */
        public static String makeFormula(TObjectIntMap<Atom> atomCounterMap, int charge) {

            // sort by symbol + A (number of mass)
            List<Atom> sortedIsotopes = new ArrayList<Atom>(atomCounterMap.keySet());

            Collections.sort(sortedIsotopes, atomComparator);

            StringBuilder sbFormula = new StringBuilder();

            /* add isotopic mass */
            for (Atom atom : sortedIsotopes) {

                int count = atomCounterMap.get(atom);

                sbFormula.append(atom.toString());

                if (count != 1) {
                    sbFormula.append(count);
                }
            }

            if (charge != 0) {
                sbFormula.append('(');

                if (charge != 1 && charge != -1) {
                    sbFormula.append((charge > 0) ? charge + "+" : -charge + "-");
                } else {
                    sbFormula.append((charge > 0) ? "+" : "-");
                }
                sbFormula.append(')');

            }

            return sbFormula.toString();
        }

        /**
         * Add all of the atoms in <code>composition</code> to this
         *
         * @param composition the composition
         * @return this builder
         */
        public Builder addAll(Composition composition) {

            for (Atom atom : composition.atomCounterMap.keySet()){

                int count = composition.atomCounterMap.get(atom);
                atomCounterMap.adjustOrPutValue(atom, count, count);
            }

            return this;
        }

        /**
         * Returns true if there are no atoms in this builder and the charge is 0
         *
         * @return true if there are no atoms in this builder and the charge is 0
         */
        public boolean isEmpty() {

            return atomCounterMap.isEmpty() && charge == 0;
        }

        /**
         * Sorts Atoms alphabetically and by neutron count
         */
        private static class MzJavaIsotopeComparatorBySymbolAndA implements Comparator<Atom>, Serializable {

            public int compare(Atom o1, Atom o2) {

                int comp = o1.getSymbol().getSymbol().compareTo(o2.getSymbol().getSymbol());
                if (comp == 0) {
                    comp = o1.getMassNumber() - o2.getMassNumber();
                }
                return comp;
            }
        }
    }


    public static Composition subtractCompositions(Composition composition1, Composition composition2){

        checkNotNull(composition1);
        checkNotNull(composition2);

        int charge = composition1.getCharge()-composition2.getCharge();
        TObjectIntMap<Atom> atomCounterMap=composition1.atomCounterMap;

        for(Atom atom : composition2.atomCounterMap.keySet()) {

            if (! composition1.atomCounterMap.containsKey(atom)){
                throw new IllegalStateException("Invalid subtraction!");
            }

            int subtraction =composition1.atomCounterMap.get(atom)- composition2.atomCounterMap.get(atom);

            if (subtraction == 0){
                atomCounterMap.remove(atom);
            }else{
                atomCounterMap.put(atom, subtraction);
            }

        }

        return new Composition(atomCounterMap, charge);

    }



}
