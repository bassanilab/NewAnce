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
import newance.psmcombiner.ModificationPSMGrouper;
import newance.psmcombiner.Psm2StringFunction;
import newance.psmcombiner.RegExpProteinGrouper;
import newance.psmconverter.MaxQuantPsmConverter;
import newance.psmconverter.PeptideMatchData;
import newance.util.ExecutableOptions;
import newance.util.NewAnceParams;
import newance.util.PsmGrouper;
import newance.util.RunTime2String;
import org.apache.commons.cli.*;

import java.io.IOException;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.regex.Pattern;

/**
 * @author Markus Müller
 */

public class MaxQuantPsmReaderTest extends ExecutableOptions {

    protected String maxquantPsmDir;
    protected Pattern maxquantPsmRegExp;
    protected String uniprotFastaFile;
    protected PsmGrouper psmGrouper;
    protected NewAnceParams params;
    protected int maxNrDisplayedPSMs;
    protected final Psm2StringFunction stringFunction;


    public MaxQuantPsmReaderTest() {

        stringFunction = new Psm2StringFunction(Psm2StringFunction.TabStringMode.MAXQUANT);
        createOptions();
    }

    public static void main(String[] args) {

        MaxQuantPsmReaderTest maxQuantPsmReaderTest =  new MaxQuantPsmReaderTest();
        try {
            maxQuantPsmReaderTest.init(args).parseOptions(args).run();
        } catch (MissingOptionException e) {
            maxQuantPsmReaderTest.checkHelpOption(args, "-h");
            maxQuantPsmReaderTest.checkVersionOption(args, NewAnceParams.getInstance().getVersion(), "-v");
            maxQuantPsmReaderTest.printOptions(args,e.getMessage());
        } catch (ParseException e) {
            maxQuantPsmReaderTest.printOptions(args, e.getMessage());
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public MaxQuantPsmReaderTest init() throws IOException {

        return this;
    }

    public int run() throws IOException {

        System.out.println("Parsing MaxQuant msms.txt and peptides.txt files ...");
        MaxQuantPsmConverter maxQuantPsmConverter = new MaxQuantPsmConverter(maxquantPsmDir, maxquantPsmRegExp);
        maxQuantPsmConverter.run();

        System.out.println("Matching peptides to fasta file: "+uniprotFastaFile+" ...");
        System.out.println("Loading protein data from fasta file: "+uniprotFastaFile+" ...");

        UniProtDB uniProtDB = new UniProtDB(uniprotFastaFile);
        long start = System.currentTimeMillis();
        maxQuantPsmConverter.addDBProteins(uniProtDB);
        System.out.println("Comet DB matching ran in " + RunTime2String.getTimeDiffString(System.currentTimeMillis() - start));

        System.out.println("RunTime after Psm parsing: " + RunTime2String.getRuntimeString(Runtime.getRuntime()));

        if (params.getCodingProtRegExp()!=null) {
            if (params.getForcedNoncanonicalProts().isEmpty())
                psmGrouper = new RegExpProteinGrouper(params.getCodingProtRegExp(), params.getProtCodingGroup(), params.getNoncanonicalGroup());
            else
                psmGrouper = new RegExpProteinGrouper(params.getCodingProtRegExp(), params.getForcedNoncanonicalProts(), params.getProtCodingGroup(), params.getNoncanonicalGroup());
        } else {
            psmGrouper = new ModificationPSMGrouper();
        }

        ConcurrentHashMap<String, List<PeptideMatchData>>  psms = maxQuantPsmConverter.getPsms();

        int cnt = 0;

        System.out.println(stringFunction.getHeader()+"\tGroup");
        for (String specID : psms.keySet()) {
            if(cnt>=maxNrDisplayedPSMs) break;

            for (PeptideMatchData psm : psms.get(specID)) {

                if(cnt>=maxNrDisplayedPSMs) break;

                System.out.print(stringFunction.getTabString(specID,psm));
                System.out.println("\t"+psmGrouper.apply(specID,psm));

                cnt++;
            }
        }

        return 0;
    }


    protected void createOptions() {

        this.cmdLineOpts = new Options();

        cmdLineOpts.addOption(Option.builder("mqD").required(false).hasArg().longOpt("maxquantPsmDir").desc("MaxQuant psm root directory. If not provided only Comet is used.").hasArg().build());
        cmdLineOpts.addOption(Option.builder("mqRE").required(false).hasArg().longOpt("maxquantPsmRegex").desc("Regular expression of MaxQuant psm files (e.g. msms\\.txt$). If not provided only Comet is used.").build());
        cmdLineOpts.addOption(Option.builder("protG").required(false).hasArg().longOpt("proteinGroup").desc("Name of group with protein coding or canonical sequences").build());
        cmdLineOpts.addOption(Option.builder("crypG").required(false).hasArg().longOpt("crypticGroup").desc("Name of group with non-canonical or cryptic sequences").build());
        cmdLineOpts.addOption(Option.builder("stdRE").required(false).hasArg().longOpt("stdProtRegExp").desc("Regular expression to match fasta name of coding proteins (e.g. sp\\||tr\\| ").build());
        cmdLineOpts.addOption(Option.builder("spRE").required(false).hasArg().longOpt("spectrumFilter").desc("If this option is set, only spectrum files that match the regexp are used").build());
        cmdLineOpts.addOption(Option.builder("exclP").required(false).hasArg().longOpt("excludeCodingProts").desc("Comma separated list of protein names to be excluded from coding prots (e.g. PGBD5_HUMAN,POGZ_HUMAN,PGBD1_HUMAN)").build());
        cmdLineOpts.addOption(Option.builder("upFa").required(false).hasArg().longOpt("uniProtFastaFile").desc("Fasta file with coding or canonical proteins (e.g. UniProt fasta file)").build());
        cmdLineOpts.addOption(Option.builder("maxR").required(false).hasArg().longOpt("maxRank").desc("Maximal rank of peptide in list of spectrum matches (rank 1 = best) (default value: 1)").build());
        cmdLineOpts.addOption(Option.builder("minZ").required(false).hasArg().longOpt("minCharge").desc("Minimal charge of PSM (default value: 1)").build());
        cmdLineOpts.addOption(Option.builder("maxZ").required(false).hasArg().longOpt("maxCharge").desc("Maximal charge of PSM (default value: 5)").build());
        cmdLineOpts.addOption(Option.builder("minL").required(false).hasArg().longOpt("minLength").desc("Minimal length of peptide (default value: 8)").build());
        cmdLineOpts.addOption(Option.builder("maxL").required(false).hasArg().longOpt("maxLength").desc("Maximal length of peptide (default value: 25)").build());
        cmdLineOpts.addOption(Option.builder("maxP").required(false).hasArg().longOpt("maxDisplayedPsms").desc("Maximal number of psms written to standard output").build());
        cmdLineOpts.addOption(Option.builder("h").required(false).hasArg(false).longOpt("help").desc("Help option for command line help").build());
        cmdLineOpts.addOption(Option.builder("v").required(false).hasArg(false).longOpt("version").desc("Version of NewAnce software").build());
    }

    @Override
    protected void check(CommandLine line) throws ParseException {

        this.params = NewAnceParams.getInstance();

        params.add("cometPsmDir", getOptionString(line,"coD"));
        params.add("cometPsmRegExp", getOptionString(line,"coRE"));

        String mqDir = getOptionString(line,"mqD");
        params.add("includeMaxQuant", mqDir.isEmpty()?"false":"true");
        params.add("maxquantPsmDir", mqDir);

        params.add("reportHistos", getOptionString(line,"repH"));
        params.add("outputDir", getOptionString(line,"outD"));
        params.add("searchFastaFile", getOptionString(line,"seFa"));
        params.add("uniprotFastaFile", getOptionString(line,"upFa"));
        params.add("doPeptideProteinGrouping", getOptionString(line,"ppG"));
        params.add("writeParamsFile", getOptionString(line,"wP"));
        params.add("readParamsFile", getOptionString(line,"rP"));
        params.add("maxRank", getOptionString(line,"maxR"));
        params.add("minCharge", getOptionString(line,"minZ"));
        params.add("maxCharge", getOptionString(line,"maxZ"));
        params.add("minPeptideLength", getOptionString(line,"minL"));
        params.add("maxPeptideLength", getOptionString(line,"maxL"));
        params.add("fdrCometThreshold", getOptionString(line,"coFDR"));
        params.add("protCodingGroup", getOptionString(line,"protG"));
        params.add("noncanonicalGroup", getOptionString(line,"noncG"));
        params.add("excludedProtPattern", getOptionString(line,"exclP"));
        params.add("forcedNoncanonicalProts", getOptionString(line,"noncP"));
        params.add("spectrumRegExp", getOptionString(line,"spRE"));
        params.add("codingProtRegExp", getOptionString(line,"protRE"));
        params.add("outputPrefix", getOptionString(line,"outP"));
        params.add("modifications", getOptionString(line,"mod"));
        params.add("minNrPsmsPerHisto", getOptionString(line,"minPH"));
        params.add("minXCorr", getOptionString(line,"minXC"));
        params.add("maxXCorr", getOptionString(line,"maxXC"));
        params.add("nrXCorrBins", getOptionString(line,"nrXCB"));
        params.add("minDeltaCn", getOptionString(line,"minDC"));
        params.add("maxDeltaCn", getOptionString(line,"maxDC"));
        params.add("nrDeltaCnBins", getOptionString(line,"nrDCB"));
        params.add("minSpScore", getOptionString(line,"minSP"));
        params.add("maxSpScore", getOptionString(line,"maxSP"));
        params.add("nrSpScoreBins", getOptionString(line,"nrSPB"));
        params.add("nrThreads", getOptionString(line,"nrTh"));
        params.add("fdrControlMethod", getOptionString(line,"fdrM"));

        params.finalize();
    }

}
