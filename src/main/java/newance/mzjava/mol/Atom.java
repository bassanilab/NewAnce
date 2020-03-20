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

/**
 * @author Oliver Horlacher
 * @version 1.0
 */
public class Atom implements Comparable<Atom>{

    private final AtomicSymbol symbol;
    private final double mass;
    private final int massNumber;
    private final double abundance;
    private final boolean defaultIsotope;

    Atom(AtomicSymbol symbol, int massNumber, double mass, double abundance, boolean defaultIsotope) {

        Preconditions.checkNotNull(symbol);
        Preconditions.checkArgument(mass > 0);
        Preconditions.checkArgument(abundance>=0);
        Preconditions.checkArgument(abundance<=1);

        this.abundance = abundance;
        this.defaultIsotope = defaultIsotope;
        this.symbol = symbol;
        this.massNumber = massNumber;
        this.mass = mass;
    }

    public AtomicSymbol getSymbol() {

        return symbol;
    }

    public double getMass() {

        return mass;
    }

    public int getMassNumber() {

        return massNumber;
    }

    public double getAbundance() {

        return abundance;
    }

    public boolean isDefaultIsotope() {

        return defaultIsotope;
    }

    @Override
    public boolean equals(Object o) {

        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Atom atom = (Atom) o;

        return massNumber == atom.massNumber && symbol == atom.symbol;

    }

    @Override
    public int hashCode() {

        return staticHash(symbol, massNumber);
    }

    public static int staticHash(AtomicSymbol symbol, int neutronCount) {

        int result = symbol.hashCode();
        result = 31 * result + neutronCount;
        return result;
    }

    @Override
    public String toString() {

        if (!isDefaultIsotope())
            return symbol.toString() + "[" + massNumber + "]";
        else
            return symbol.toString();
    }

    @Override
    public int compareTo(Atom o) {

        return Double.compare(mass, o.mass);
    }
}
