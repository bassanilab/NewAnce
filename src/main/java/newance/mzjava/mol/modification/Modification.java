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
package newance.mzjava.mol.modification;

import newance.mzjava.mol.Composition;
import newance.mzjava.mol.Mass;
import newance.mzjava.mol.Weighable;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static com.google.common.base.Preconditions.checkNotNull;

/**
 * @author Oliver Horlacher
 * @version 1.0
 */
public class Modification implements Weighable {

    public static final Pattern parsePattern = Pattern.compile("(?:(.*):)?((?:[A-Z][a-z]?(?:\\[\\d+\\])?[+-]?\\d*)+(?:\\(\\d*[+-]\\))?)");

    private final String label;
    private final Mass mass;

    public Modification(String label, Mass mass) {

        checkNotNull(label);
        this.label = label;

        checkNotNull(mass);
        this.mass = mass;
    }

    public String getLabel() {

        return label;
    }

    public Mass getMass() {

        return mass;
    }

    @Override
    public double getMolecularMass() {

        return mass.getMolecularMass();
    }

    @Override
    public boolean equals(Object o) {

        if (this == o) return true;
        if (!(o instanceof Modification)) return false;

        Modification that = (Modification) o;

        return label.equals(that.label) && mass.equals(that.mass);
    }

    @Override
    public int hashCode() {

        int result = label.hashCode();
        result = 31 * result + mass.hashCode();
        return result;
    }

    @Override
    public String toString() {

        String formula = mass.toString();
        if (formula.equals(label)) {

            return formula;
        } else {
            return String.format("%s:%s", label, formula);
        }
    }

    /**
     * Parse a modification.  The format is (label:)?formula.  An optional label followed by a formula
     * delimited by :
     *
     * @param s the string to parse
     * @return the parsed Modification
     * @throws IllegalArgumentException - if s is not a valid modification
     */
    public static Modification parseModification(String s) {

        Matcher matcher = parsePattern.matcher(s);

        if(!matcher.matches()) throw new IllegalArgumentException(s + " is not a valid format for modifications");

        String label = matcher.group(1);
        String formula = matcher.group(2);

        if(label == null) label = formula;

        return new Modification(label, Composition.parseComposition(formula));
    }
}
