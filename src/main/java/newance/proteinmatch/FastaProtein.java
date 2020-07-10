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

public class FastaProtein {
    final protected String proteinID;
    protected char[] sequence;


    public FastaProtein(String proteinID, String sequence) {
        this.proteinID = proteinID;
        this.sequence = sequence.toCharArray();
    }

    public FastaProtein(String proteinID) {
        this.proteinID = proteinID;
        this.sequence = null;
    }

    public String getProteinID() {
        return proteinID;
    }

    public String getSequence() {
        return new String(sequence);
    }

    public char[] getCharArray() { return sequence; }

    public void setSequence(String sequence) {
        this.sequence = sequence.toCharArray();
    }

    public void setSequence(char[] sequence) {
        this.sequence = sequence;
    }


    public String getWTSequence(int start, int length, int offset) {

        if (start < 0 || start >= sequence.length) return "";

        start = (start-offset < 0)?0:start-offset;
        length = (start+length+2*offset < sequence.length)?length+2*offset:sequence.length-start;

        return String.copyValueOf(sequence, start, length);
    }

}
