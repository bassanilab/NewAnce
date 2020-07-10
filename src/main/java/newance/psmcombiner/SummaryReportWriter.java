/*
Copyright (C) SIB - Swiss Institute of Bioinformatics, Lausanne, Switzerland
Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
*/

package newance.psmcombiner;

import newance.psmconverter.PeptideSpectrumMatch;
import newance.util.NewAnceParams;
import newance.util.ProcessPsmUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;

/**
 * @author Markus Müller
 */

public class SummaryReportWriter {
    private final String fileName;
    private BufferedWriter reportWriter;
    private boolean includeMaxQuant;


    public SummaryReportWriter(String fileName, boolean includeMaxQuant) {

        this.fileName = fileName;
        this.includeMaxQuant = includeMaxQuant;

        try {
            reportWriter = new BufferedWriter(new FileWriter(new File(fileName)));
            writeHeader();

        } catch (IOException e) {
            System.out.println("Unable to write report to file "+ fileName);
            reportWriter = null;
        }
    }

    private void writeHeader() throws IOException {
        if (includeMaxQuant) reportWriter.write("Group\tNrCometPSMs\tNrMaxQuantPSMs\tNrCombinedPSMs\tNrCometPepts\tNrMaxQuantPepts\tNrCombinedPepts\n");
        else reportWriter.write("Group\tNrCometPSMs\tNrCometPepts\n");
    }

    public void write(String group,
                      ConcurrentHashMap<String, List<PeptideSpectrumMatch>> combinedPsms,
                      ConcurrentHashMap<String, List<PeptideSpectrumMatch>>  cometPsms,
                      ConcurrentHashMap<String, List<PeptideSpectrumMatch>>  maxQuantPsms) {

        if (reportWriter==null) return;

        ConcurrentHashMap<String, List<PeptideSpectrumMatch>>  noDecoycombinedPsms = ProcessPsmUtils.removeDecoys(combinedPsms);
        ConcurrentHashMap<String, List<PeptideSpectrumMatch>>  noDecoyCometPsms = ProcessPsmUtils.removeDecoys(cometPsms);
        ConcurrentHashMap<String, List<PeptideSpectrumMatch>>  noDecoyMaxQuantPsms = ProcessPsmUtils.removeDecoys(maxQuantPsms);

        try {
            reportWriter.write(group+"\t"+ProcessPsmUtils.countPsms(noDecoyCometPsms)+"\t"+ ProcessPsmUtils.countPsms(noDecoyMaxQuantPsms)+"\t"+
                    ProcessPsmUtils.countPsms(noDecoycombinedPsms)+"\t"+ProcessPsmUtils.countUniquePeptides(noDecoyCometPsms)+"\t"+
                    ProcessPsmUtils.countUniquePeptides(noDecoyMaxQuantPsms)+"\t"+ ProcessPsmUtils.countUniquePeptides(noDecoycombinedPsms)+"\n");

        } catch (IOException e) {
        }
    }


    public void write(String group, ConcurrentHashMap<String, List<PeptideSpectrumMatch>> cometPsms) {

        ConcurrentHashMap<String, List<PeptideSpectrumMatch>>  noDecoyCometPsms = ProcessPsmUtils.removeDecoys(cometPsms);

        try {
            reportWriter.write(group+"\t"+ProcessPsmUtils.countPsms(noDecoyCometPsms)+"\t"+ProcessPsmUtils.countUniquePeptides(noDecoyCometPsms)+"\n");

        } catch (IOException e) {
        }
    }

    public void close() {
        try {
            if (reportWriter!=null) reportWriter.close();
            reportWriter = null;
        } catch (IOException e) {
            System.out.println("Unable to close file "+ fileName);
        }
    }

}
