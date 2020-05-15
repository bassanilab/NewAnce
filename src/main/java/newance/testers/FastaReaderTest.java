/*
Copyright (C) SIB - Swiss Institute of Bioinformatics, Lausanne, Switzerland
Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
*/

package newance.testers;

import newance.proteinmatch.FastaProtein;
import newance.proteinmatch.UniProtDB;
import newance.proteinmatch.UniProtProtein;
import newance.util.ExecutableOptions;
import newance.util.NewAnceParams;
import newance.util.RunTime2String;
import org.apache.commons.cli.*;

import java.io.IOException;

/**
 * @author Markus Müller
 */

public class FastaReaderTest extends ExecutableOptions {

    protected String fastaFile;
    protected NewAnceParams params;
    protected int maxNrDisplayedPSMs;


    public FastaReaderTest() {

        createOptions();
    }

    public static void main(String[] args) {

        FastaReaderTest fastaReaderTest =  new FastaReaderTest();
        try {
            fastaReaderTest.init(args).parseOptions(args).run();
        } catch (MissingOptionException e) {
            fastaReaderTest.checkHelpOption(args, "-h");
            fastaReaderTest.checkVersionOption(args, NewAnceParams.getInstance().getVersion(), "-v");
            fastaReaderTest.printOptions(args,e.getMessage());
        } catch (ParseException e) {
            fastaReaderTest.printOptions(args, e.getMessage());
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public FastaReaderTest init() throws IOException {

        return this;
    }

    public int run() throws IOException {

        long start = System.currentTimeMillis();
        UniProtDB uniProtDB = new UniProtDB(fastaFile);
        System.out.println("Reading fasta in " + RunTime2String.getTimeDiffString(System.currentTimeMillis() - start));

        int cnt = 0;

        System.out.println("NewAnceID\tProteinID\tProtName\tGeneName\tDescription\tDB\tSequence\n");

        for (FastaProtein prot : uniProtDB.getProteins()) {
            if(maxNrDisplayedPSMs>=0 && cnt>=maxNrDisplayedPSMs) break;

            UniProtProtein protein = (UniProtProtein) prot;
            System.out.print(protein.toString()+"\t"+prot.getProteinID()+"\t"+protein.getUniProtName()+"\t"+
                    protein.getGeneName()+"\t"+protein.getDescription()+"\t");
            System.out.println(protein.getDbFlag()+"\t"+prot.getSequence());
            cnt++;
        }

        return 0;
    }


    protected void createOptions() {

        this.cmdLineOpts = new Options();

        cmdLineOpts.addOption(Option.builder("fa").required(false).hasArg().longOpt("fastaFile").desc("Fasta file to be tested").build());
        cmdLineOpts.addOption(Option.builder("maxP").required(false).hasArg().longOpt("maxDisplayedPsms").desc("Maximal number of psms written to standard output").build());
        cmdLineOpts.addOption(Option.builder("h").required(false).hasArg(false).longOpt("help").desc("Help option for command line help").build());
        cmdLineOpts.addOption(Option.builder("v").required(false).hasArg(false).longOpt("version").desc("Version of NewAnce software").build());
    }

    @Override
    protected void check(CommandLine line) throws ParseException {

        this.params = NewAnceParams.getInstance();

        String maxPStr = getOptionString(line,"maxP");
        maxNrDisplayedPSMs = (maxPStr.isEmpty())?-1:Integer.parseInt(maxPStr);

        fastaFile = getOptionString(line,"fa");

    }

}
