/*
Copyright (C) SIB - Swiss Institute of Bioinformatics, Lausanne, Switzerland
Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
*/

package newance.proteinmatch;

/**
 * @author Markus Müller
 */

public class PeptideUniProtSequenceMatch {

    protected final int start;
    protected final int end;
    protected final FastaProtein protein;
    protected final char[] peptide;


    public PeptideUniProtSequenceMatch(int start, int end, char[] peptide, FastaProtein protein) {
        this.start = start;
        this.end = end;

        this.protein = protein;
        this.peptide = peptide;
    }

    public String toString() {
        return peptide+"\t"+protein.getProteinID()+"\t"+start+"\t"+end;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public FastaProtein getProtein() {
        return protein;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        PeptideUniProtSequenceMatch that = (PeptideUniProtSequenceMatch) o;

        if (start != that.start) return false;
        if (end != that.end) return false;
        if (!protein.equals(that.protein)) return false;
        return peptide.equals(that.peptide);
    }

    @Override
    public int hashCode() {
        int result = start;
        result = 31 * result + end;
        result = 31 * result + protein.hashCode();
        result = 31 * result + peptide.hashCode();
        return result;
    }
}

