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

package newance.testers;

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

        for (UniProtProtein prot : uniProtDB.getProteins()) {
            if(maxNrDisplayedPSMs>=0 && cnt>=maxNrDisplayedPSMs) break;

            System.out.print(prot.getFastaUniProtName()+"\t"+prot.getUniProtAC()+"\t"+prot.getUniProtName()+"\t"+prot.getGeneName()+"\t"+prot.getDescription()+"\t");
            System.out.println(prot.getDbFlag()+"\t"+prot.getSequence());
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
