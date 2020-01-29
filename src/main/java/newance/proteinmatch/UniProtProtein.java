/**
 * Copyright (C) 2019, SIB/LICR. All rights reserved
 *
 * SIB, Swiss Institute of Bioinformatics
 * Ludwig Institute for Cancer Research (LICR)
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer. Redistributions in binary
 * form must reproduce the above copyright notice, this list of conditions and
 * the following disclaimer in the documentation and/or other materials provided
 * with the distribution. Neither the name of the SIB/LICR nor the names of
 * its contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL SIB/LICR BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
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
