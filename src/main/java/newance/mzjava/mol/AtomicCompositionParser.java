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

import gnu.trove.map.TObjectIntMap;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Parsing atomic composition.
 *
 * <h2>Atom definition</h2>
 * An atom is defined by a symbol with an optional explicit mass number (nucleon number)<br>
 * Example: hydrogen atom could be described like "H" or "H[1]"
 *
 * <h2>Atom occurrence definition</h2>
 * The number of atoms directly follows the atom definition else there is only one occurrence.
 * This number can be positive or negative, signed or not.<br>
 * Example 1: "O2" or "O+2" means 2 atoms of oxygen (main isotope O[16])<br>
 * Example 2: "H-1" defines a loss of H
 *
 * <h2>Composition definition</h2>
 * A composition is an assembling of atoms (or groups of atom) and is simply described as a sequence
 * of previously defined atoms.<br>
 * Example: "CH4", "H2O" or "C2H6OH"
 * <p>
 * A composition may also contains groups of atoms simplifying some composition definition<br>
 * Example: "CH3(CH2)2OH"<br>
 * Observe the two occurrences of CH2 group in the previous propanol molecule.
 *
 * <h2>Composition charge definition</h2>
 * A composition can be neutral or charged. The charge has to close the composition.
 * The charge, surrounded with parenthesis, may be positive or negative with the sign placed at the end<br>
 * Example: "H3O(+)", "OH(-)", "SO4(2-)"
 *
 * <h2>Grammar</h2>
 *
 * <pre>
 * Composition      :   AtomicExpression+ ('(' Charge ')') ?
 * AtomicExpression :   HomoAtoms | AtomGroup
 * HomoAtoms        :   AtomicSymbol ('[' MassNumber ']')? Count?
 * AtomGroup        :   '(' HomoAtoms ')' Count
 * AtomicSymbol     :   [A-Z][a-z]?
 * MassNumber       :   [0-9]+
 * Count            :   [+-]?[0-9]+
 * Charge           :   [0-9]* [+-]
 * </pre>
 *
 * <h2>Remark</h2>
 * Sometimes things can be a bit tricky with the overall charge state when atoms are removed.
 * Let's take an example by describing an expression for proton removal.
 * <p>
 * An intuitive way of writing it would be "H-1(+)" but it is wrong as it leads to subtracting an hydrogen
 * and an electron (-H-e) whereas it should be subtracting an hydrogen and adding an electron (H+ = H-1e => -H+ = -H+1e).
 * <p>
 * As a result the right way of defining it is "H-1(-)" with the correct overall deduced charge given to the parser.
 *
 * @author nikitin
 * @version 1.0
 */
public final class AtomicCompositionParser {

    private static final AtomicCompositionParser INSTANCE = new AtomicCompositionParser();

    private static final Pattern CHARGE_PATTERN = Pattern.compile("\\((\\d*)([-+])\\)$");
    private static final Pattern GROUP_PATTERN = Pattern.compile("\\(([^)]+)\\)([-+]?\\d+)?");
    private static final Pattern ELEMENT_PATTERN = Pattern.compile("([A-Z][a-z]?)(?:\\[(\\d{1,3})\\])?([+-]?\\d+)?\\s*");

    private final PeriodicTable periodicTable = PeriodicTable.getInstance();

    private int charge;

    private AtomicCompositionParser() {
    }

    /**
     * @return the singleton
     */
    public static AtomicCompositionParser getInstance() {
        return INSTANCE;
    }

    /**
     * Parses the given string content and transmit data to builder.
     *
     * @param content the data string to parse.
     */
    public int parse(String content, TObjectIntMap<Atom> composition) {

        composition.clear();
        charge = 0;

        /* parse charge */
        // i.e. content = CH3(CH2)2O(CH2)3OH(2-)
        content = parseCharge(content);

        /* parse molecule */
        // i.e. content = CH3(CH2)2O(CH2)3OH
        content = parseGroup(content, composition);

        /* parse last atoms */
        // i.e. content = CH3OOH
        int offset = parseMolecule(content, 1, composition);

        if (offset != 0) {
            throw new IllegalArgumentException("molecule " + content + " not found");
        }

        return charge;
    }

    private String parseCharge(String content) {

        /* match molecular charge pattern */
        Matcher matcher = CHARGE_PATTERN.matcher(content);

        if (matcher.find()) {
            int chargeNb = 1;
            String chargeType;

            if (!matcher.group(1).isEmpty()) {
                chargeNb = Integer.parseInt(matcher.group(1));
            }
            chargeType = matcher.group(2);

            if ("+".equals(chargeType)) {

                this.charge = chargeNb;
            } else {

                this.charge = -chargeNb;
            }

            // remove the charge part
            return content.substring(0, matcher.start());
        }
        return content;
    }

    private String parseGroup(String content, TObjectIntMap<Atom> composition) {

        StringBuilder sb = new StringBuilder();

        Matcher groupMatcher = GROUP_PATTERN.matcher(content);

        int from = 0;
        int offset;
        while (groupMatcher.find()) {

            int moleculeStartPos = groupMatcher.start(1);

            String group =
                    content.substring(moleculeStartPos, groupMatcher.end(1));
            int groupNumber = 1;

            String occNumberString = groupMatcher.group(2);

            if (occNumberString != null) {
                groupNumber = Integer.parseInt(occNumberString);
            }

            offset = parseMolecule(group, groupNumber, composition);
            if (offset != 0) {
                throw new IllegalArgumentException("molecule " + content + " not found");
            }

            sb.append(content.substring(from, moleculeStartPos - 1));

            if (groupMatcher.group(2) != null)
                from = groupMatcher.end(2);
            else
                from = groupMatcher.end(1)+1;
        }

        if (from < content.length())
            sb.append(content.substring(from));

        return sb.toString();
    }

    /**
     * Parse formula: if at least one explicit isotope, every atoms are isotopes
     * else all element are atoms.
     *
     * @param formula the formula
     * @param repetition the repetition
     * @param composition the composition
     * @return an int
     */
    private int parseMolecule(String formula, int repetition, TObjectIntMap<Atom> composition) {

        // 3 groups (atom, (mass number), (quantity))
        Matcher eltMatcher = ELEMENT_PATTERN.matcher(formula);

        int lastMatchEnd = 0;
        while (eltMatcher.find()) {

            AtomicSymbol atomicSymbol = AtomicSymbol.valueOf(eltMatcher.group(1));
            String massNumberGroup = eltMatcher.group(2);
            String occurrencesGroup = eltMatcher.group(3);

            int occNumber = occurrencesGroup != null ? Integer.parseInt(occurrencesGroup) : 1;

            occNumber *= repetition;

            if (lastMatchEnd != eltMatcher.start()) {
                throw new IllegalArgumentException("cannot recognize "
                        + formula.substring(lastMatchEnd, eltMatcher.start())
                        + " in " + formula);
            }

            Atom atom;
            if (massNumberGroup == null) {

                atom = periodicTable.getAtom(atomicSymbol);
            } else {

                atom = periodicTable.getAtom(atomicSymbol, Integer.parseInt(massNumberGroup));
            }
            composition.adjustOrPutValue(atom, occNumber, occNumber);

            lastMatchEnd = eltMatcher.end();

        }

        return formula.length() - lastMatchEnd;
    }
}
