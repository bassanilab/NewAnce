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

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

/**
 * @author Markus Müller
 */

public class UniProtDB extends FastaDB {

    protected final ConcurrentHashMap<String, Map<String, List<PeptideUniProtSequenceMatch>>> searchedPeptides;


    public UniProtDB() {

        super(6);
        searchedPeptides = new ConcurrentHashMap<>();
    }

    public UniProtDB(String[] fastaLines) {

        super( fastaLines, 6);
        searchedPeptides = new ConcurrentHashMap<>();
    }

    public UniProtDB(String fastaFileName) {

        super(fastaFileName, 6);
        searchedPeptides = new ConcurrentHashMap<>();
    }

    public UniProtDB(List<String> fastaFileNames) {

        super(fastaFileNames, 6);
        searchedPeptides = new ConcurrentHashMap<>();
    }


    protected FastaProtein parseHeader(String header) {

        if (header.startsWith(">sp|") || header.startsWith(">tr|")) {
            return parseUniProtHeader(header);
        } else {
            return parseOtherHeader(header);
        }
    }

    //>ENSP00000334393.3 \VariantSimple=(141|A)
    //>LINE_L1_L1MC4_23_154353601_154353761.1-117
    protected static FastaProtein parseOtherHeader(String header) {

        String[] fields = header.split("\\s");
        String id = fields[0].substring(1);
        Map<String,String> headerFieldMap = new HashMap<>();

        UniProtProtein protein = new UniProtProtein(
                id,
                "",
                (fields.length>1)?fields[1]:"",
                "",
                "other");

        return protein;
    }


    //>tr|A0A024R3B9|A0A024R3B9_HUMAN Alpha-crystallin B chain OS=Homo sapiens GN=CRYAB PE=3 SV=1
    //>sp|A0A183|LCE6A_HUMAN Late cornified envelope protein 6A OS=Homo sapiens GN=LCE6A PE=2 SV=1
    protected static FastaProtein parseUniProtHeader(String header) {

        String[] fields = header.split("\\|");

        int idx = fields[2].indexOf(" ");
        int geneNameStart = fields[2].indexOf("GN=")+3;
        int geneNameEnd = fields[2].indexOf(" ",geneNameStart);
        if (geneNameEnd<0) geneNameEnd = fields[2].length();

        UniProtProtein protein = new UniProtProtein(
                fields[1],
                fields[2].substring(0,idx),
                fields[2].substring(idx+1),
                fields[2].substring(geneNameStart,geneNameEnd),
                fields[0].substring(1,3));

        return protein;
    }

    public void clear() {
        this.fastaFileNames.clear();
        this.proteins.clear();
        this.acIndex.clear();
        this.sequenceIndex.clear();
        this.searchedPeptides.clear();
    }


}
