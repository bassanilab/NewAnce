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

import java.io.Serializable;

/**
 * @author Markus Müller
 */

public class UniProtProtein implements Serializable{
    final protected String uniProtAC;
    final protected String uniProtName;
    final protected String description;
    final protected char[] sequence;
    final protected String geneName;
    final protected String dbFlag;


    public UniProtProtein(String uniProtAC, String uniProtName, String description, String geneName, String dbFlag, String sequence) {
        this.uniProtAC = uniProtAC;
        this.uniProtName = uniProtName;
        this.description = description;
        this.geneName = geneName;
        this.dbFlag = dbFlag;
        this.sequence = sequence.toCharArray();
    }

    public String getUniProtAC() {
        return uniProtAC;
    }

    public String getUniProtName() {
        return uniProtName;
    }

    public String getDescription() {
        return description;
    }

    public String getSequence() {
        return new String(sequence);
    }

    public char[] getCharArray() { return sequence; }

    public String getGeneName() {
        return geneName;
    }

    public String getDbFlag() {
        return dbFlag;
    }

    public String getFastaUniProtName() {
        if (dbFlag.equals("sp") || dbFlag.equals("tr"))
            return dbFlag+"|"+uniProtAC+"|"+uniProtName;
        else
            return uniProtAC;
    }
}
