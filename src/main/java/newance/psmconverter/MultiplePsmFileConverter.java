/*
Copyright (C) SIB - Swiss Institute of Bioinformatics, Lausanne, Switzerland
Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
*/

package newance.psmconverter;

import newance.proteinmatch.UniProtDB;

import java.io.IOException;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.regex.Pattern;

/**
 * @author Markus Müller
 */

public abstract class MultiplePsmFileConverter {

    protected final String psmRootDirName;
    protected final Pattern regex;
    protected final ConcurrentHashMap<String,List<PeptideSpectrumMatch>> psms;

    public MultiplePsmFileConverter(String psmRootDirName, Pattern regex) {

        this.psmRootDirName = psmRootDirName;
        this.regex = regex;
        this.psms = new ConcurrentHashMap<>();
    }


    public abstract void run() throws IOException;

    public void addDBProteins(UniProtDB uniProtDB) {

        if (uniProtDB==null) return;

        psms.forEach(100000,new AddUniProtIDs2Psm(uniProtDB));
    }

    public void reportDBProteins() {

        psms.forEach((id,psm) -> System.out.println(id+", "+psm.get(0).getPeptide().toString()+", "+psm.get(0).getProteinIDs().toString()));
    }

    public ConcurrentHashMap<String, List<PeptideSpectrumMatch>> getPsms() {
        return psms;
    }
}
