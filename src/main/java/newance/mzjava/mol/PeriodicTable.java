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

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Periodic table of the elements. Is used to get Atom instances.
 *
 * @author nikitin
 * @author Oliver Horlacher
 * @version 1.0
 */
public class PeriodicTable {

    private static final PeriodicTable instance;

    public static final double ELECTRON_MASS = 0.0005489243;

    public static final double H_MASS;
    public static final double O_MASS;
    public static final double C_MASS;
    public static final double C13_MASS;
    public static final double N_MASS;
    public static final double P_MASS;
    public static final double S_MASS;
    public static final double K_MASS;
    public static final double Ca_MASS;
    public static final double Fe_MASS;
    public static final double Mg_MASS;
    public static final double Na_MASS;
    public static final double Zn_MASS;
    public static final double Mn_MASS;
    public static final double Cu_MASS;
    public static final double I_MASS;

    public static final Atom H;
    public static final Atom O;
    public static final Atom C;
    public static final Atom C13;
    public static final Atom N;
    public static final Atom P;
    public static final Atom S;
    public static final Atom K;
    public static final Atom Ca;
    public static final Atom Fe;
    public static final Atom Mg;
    public static final Atom Na;
    public static final Atom Zn;
    public static final Atom Mn;
    public static final Atom Cu;
    public static final Atom I;

    public static PeriodicTable getInstance() {

        return instance;
    }

    static {

        instance = new PeriodicTable();

        H = instance.getAtom(AtomicSymbol.H);
        O = instance.getAtom(AtomicSymbol.O);
        C = instance.getAtom(AtomicSymbol.C);
        C13 = instance.getAtom(AtomicSymbol.C, 13);
        N = instance.getAtom(AtomicSymbol.N);
        P = instance.getAtom(AtomicSymbol.P);
        S = instance.getAtom(AtomicSymbol.S);
        K = instance.getAtom(AtomicSymbol.K);
        Ca = instance.getAtom(AtomicSymbol.Ca);
        Na = instance.getAtom(AtomicSymbol.Na);
        Fe = instance.getAtom(AtomicSymbol.Fe);
        Mg = instance.getAtom(AtomicSymbol.Mg);
        Zn = instance.getAtom(AtomicSymbol.Zn);
        Mn = instance.getAtom(AtomicSymbol.Mn);
        Cu = instance.getAtom(AtomicSymbol.Cu);
        I = instance.getAtom(AtomicSymbol.I);

        H_MASS = H.getMass();
        O_MASS = O.getMass();
        C_MASS = C.getMass();
        C13_MASS = C13.getMass();
        N_MASS = N.getMass();
        P_MASS = P.getMass();
        S_MASS = S.getMass();
        K_MASS = K.getMass();
        Ca_MASS = Ca.getMass();
        Na_MASS = Na.getMass();
        Fe_MASS = Fe.getMass();
        Mg_MASS = Mg.getMass();
        Zn_MASS = Zn.getMass();
        Mn_MASS = Mn.getMass();
        Cu_MASS = Cu.getMass();
        I_MASS = I.getMass();
    }

    private final Map<Integer, Atom> atomMap = new HashMap<Integer, Atom>();
    private final Map<AtomicSymbol, Atom> defaultIsotopeMap = new HashMap<AtomicSymbol, Atom>();

    /**
     * Construct a PeriodicTable by reading the periodic-table.dat file
     */
    private PeriodicTable() {


        InputStream in = getClass().getResourceAsStream("periodic-table.dat");
        InputStreamReader isr;
        try {

            isr = new InputStreamReader(in, "UTF8");
        } catch (UnsupportedEncodingException e) {

            throw new IllegalStateException(e);
        }

        /* regular expression definitions for atoms and isotopes */
        final Pattern atomPattern =
                Pattern.compile("^([A-Z]+)\\s+([A-Z][a-z]?)");

        final Pattern isotopePattern =
                Pattern.compile("[A-Z][a-z]?\\((\\d+)\\)\\s+(\\S+)\\s+(\\S+)\\s*");

        BufferedReader br = new BufferedReader(isr);
        try {

            String line = br.readLine();

            // parse each line
            while (line != null) {
                // ignore comment lines
                if (!line.matches("//.+")) {

                    // parse atom line

                    /* match atom pattern */
                    Matcher matcher = atomPattern.matcher(line);
                    if (!matcher.find()) {
                        throw new IllegalStateException("bad atom format on " + line);
                    }

                    AtomicSymbol symbol = AtomicSymbol.valueOf(matcher.group(2));

                    /* match isotope pattern */
                    matcher = isotopePattern.matcher(line);

                    extractIsotopes(matcher, symbol, true);
                }
                // next line
                line = br.readLine();
            }
            br.close();
        } catch (IOException e) {

            try {

                br.close();
            } catch (IOException e1) {

                throw new IllegalStateException(e1);
            }

            throw new IllegalStateException("IOException caught while reading periodic-table.dat", e);
        }
    }

    private void extractIsotopes(Matcher matcher, AtomicSymbol symbol, boolean defaultIsotope) {

        while (matcher.find()) {

            int massNumber = Integer.parseInt(matcher.group(1));
            double mass = Double.parseDouble(matcher.group(2));
            double abundance = Double.parseDouble(matcher.group(3));

            int hash = Atom.staticHash(symbol, massNumber);

            Atom atom = new Atom(symbol, massNumber, mass, abundance / 100., defaultIsotope);
            atomMap.put(hash, atom);

            if (defaultIsotope) {

                defaultIsotopeMap.put(symbol, atom);
            }
            defaultIsotope = false;
        }
    }

    /**
     * Return the atom for the atomicSymbol with the number of neutrons given by massNumber
     *
     * @param atomicSymbol the atomic symbol
     * @param massNumber   the number of neutrons
     * @return the atom for the atomicSymbol with the number of neutrons given by massNumber
     */
    public Atom getAtom(AtomicSymbol atomicSymbol, int massNumber) {

        Atom atom = atomMap.get(Atom.staticHash(atomicSymbol, massNumber));
        if (atom != null) {

            return atom;
        }

        throw new IllegalArgumentException(massNumber + " is not a valid argument for getting the atom for the symbol " + atomicSymbol);
    }

    /**
     * Return the default isotope for the <code>atomicSymbol</code>
     *
     * @param atomicSymbol the atomic symbol
     * @return the default isotope for the <code>atomicSymbol</code>
     */
    public Atom getAtom(AtomicSymbol atomicSymbol) {

        Preconditions.checkNotNull(atomicSymbol);

        return defaultIsotopeMap.get(atomicSymbol);
    }

    /**
     * Parse the symbol and return the corresponding Atom. For example,
     * <ul>
     * <li>if symbol {@code C} given then an instance of {@code Atom
     * C[12]} returned</li>
     * <li>if symbol {@code C[13]} given then an instance of {@code
     * Atom C[13]} instance returned</li>
     * </ul>
     * <p/>
     * <h2>Atom definition</h2>
     * An atom is defined by a symbol with an optional explicit mass number (nucleon number)<br>
     * Example: hydrogen atom could be described like "H", "H[1]"
     *
     * @param symbol the symbol
     * @return the atom for the symbol
     * @throws IllegalArgumentException if bad element format
     */
    public Atom getAtom(String symbol) {

        int index = symbol.indexOf('[');

        try {
            if (index == -1) {

                return getAtom(AtomicSymbol.valueOf(symbol));
            } else {

                return getAtom(AtomicSymbol.valueOf(symbol.substring(0, index)), Integer.parseInt(symbol.substring(index + 1, symbol.length() - 1)));
            }
        } catch (NumberFormatException e) {

            throw new IllegalArgumentException(symbol + ": invalid symbol", e);
        } catch (IllegalArgumentException e) {

            throw new IllegalArgumentException(symbol + ": invalid symbol", e);
        }
    }

    /**
     * Return all the atoms for the given symbol sorted by mass
     * <p/>
     * For example AtomicSymbol.C will return {C12, C13, C14}
     *
     * @param atomicSymbol the symbol
     * @return all the atoms for the given symbol sorted by mass
     */
    public List<Atom> getAtoms(AtomicSymbol atomicSymbol) {

        List<Atom> atoms = new ArrayList<Atom>(4);

        for (Atom atom : atomMap.values()) {

            if (atom.getSymbol() == atomicSymbol) atoms.add(atom);
        }

        Collections.sort(atoms);

        return atoms;
    }
}
