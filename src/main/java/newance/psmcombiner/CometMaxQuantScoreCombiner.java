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

package newance.psmcombiner;

import newance.util.*;
import org.apache.commons.cli.*;
import newance.proteinmatch.OccamRazorSpectrumCounter;
import newance.proteinmatch.UniProtDB;
import newance.psmconverter.MaxQuantPsmConverter;
import newance.psmconverter.PeptideMatchData;
import newance.psmconverter.CometPsmConverter;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

/**
 * @author Markus Müller
 */

public class CometMaxQuantScoreCombiner extends ExecutableOptions {

    protected PsmGrouper psmGrouper;
    protected NewAnceParams params;
    protected SpectrumAccumulator spectrumAccumulator;
    protected String readParamsFile;


    public CometMaxQuantScoreCombiner() {

        createOptions();
        spectrumAccumulator = new SpectrumAccumulator();
    }

    public static void main(String[] args) {

        CometMaxQuantScoreCombiner cometScoreCombiner =  new CometMaxQuantScoreCombiner();
        try {
            cometScoreCombiner.init(args).parseOptions(args).run();
        } catch (MissingOptionException e) {
        } catch (ParseException e) {
            cometScoreCombiner.printOptions(args, e.getMessage());
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public int run() throws IOException {

        System.out.println("Parsing Comet pep.xml files ...");
        CometPsmConverter cometPsmConverter = new CometPsmConverter(params.getCometPsmDir(), params.getCometPsmRegExp());
        cometPsmConverter.run();

        MaxQuantPsmConverter maxQuantPsmConverter = null;
        if (params.isIncludeMaxQuant()) {
            System.out.println("Parsing MaxQuant msms.txt and peptides.txt files ...");
            maxQuantPsmConverter = new MaxQuantPsmConverter(params.getMaxquantPsmDir(), params.getMaxquantPsmRegExp());
            maxQuantPsmConverter.run();
        }

        System.out.println("Matching peptides to fasta file: "+params.getUniprotFastaFile()+" ...");
        System.out.println("Loading protein data from fasta file: "+params.getUniprotFastaFile()+" ...");

        UniProtDB uniProtDB = new UniProtDB(params.getUniprotFastaFile());
        long start = System.currentTimeMillis();
        cometPsmConverter.addDBProteins(uniProtDB);
        System.out.println("Comet DB matching ran in " + RunTime2String.getTimeDiffString(System.currentTimeMillis() - start));

        if (params.isIncludeMaxQuant()) {
            start = System.currentTimeMillis();
            maxQuantPsmConverter.addDBProteins(uniProtDB);
            System.out.println("MaxQuant DB matching ran in " + RunTime2String.getTimeDiffString(System.currentTimeMillis() - start));
        }

        GroupedFDRCalculator groupedFDRCalculator = buildGroupedFDRCalculator(cometPsmConverter.getPsms());

        if (params.getFdrControlMethod().equals("global")) {
            controlFDRGlobally(groupedFDRCalculator, cometPsmConverter, maxQuantPsmConverter);
        } else {
            controlFDRGroupwise(groupedFDRCalculator, cometPsmConverter, maxQuantPsmConverter);
        }

        writePeptideProteinGroupReport(uniProtDB);

//        groupedFDRCalculator.calcPValuesFDR();

        System.out.println("RunTime after Psm parsing: " + RunTime2String.getRuntimeString(Runtime.getRuntime()));


        return 0;
    }

    protected GroupedFDRCalculator buildGroupedFDRCalculator(ConcurrentHashMap<String, List<PeptideMatchData>> allPsms ) {

        if (params.getCodingProtRegExp()!=null) {
            if (params.getForcedNoncanonicalProts().isEmpty())
                psmGrouper = new RegExpProteinGrouper(params.getCodingProtRegExp(), params.getProtCodingGroup(), params.getNoncanonicalGroup());
            else
                psmGrouper = new RegExpProteinGrouper(params.getCodingProtRegExp(), params.getForcedNoncanonicalProts(), params.getProtCodingGroup(), params.getNoncanonicalGroup());
        } else {
            psmGrouper = new ModificationPSMGrouper();
        }

        GroupedFDRCalculator groupedFDRCalculator = new GroupedFDRCalculator(psmGrouper);
        groupedFDRCalculator.addAll(allPsms);
        groupedFDRCalculator.setCanCalculateFDR(params.getMinNrPsmsPerHisto());
        groupedFDRCalculator.calcClassProbs();
        if (params.isReportHistos()) groupedFDRCalculator.writeHistograms(params.getOutputDir()+File.separator+"histos", params.getOutputPrefix());
        groupedFDRCalculator.smoothHistogram(params.getSmoothDegree());
        groupedFDRCalculator.calcLocalFDR();
        if (params.isReportHistos()) groupedFDRCalculator.writeHistograms(params.getOutputDir()+File.separator+"histos", params.getOutputPrefix()+"_smoothed");

        return groupedFDRCalculator;
    }

    protected void controlFDRGlobally(GroupedFDRCalculator groupedFDRCalculator, CometPsmConverter cometPsmConverter, MaxQuantPsmConverter maxQuantPsmConverter) {

        float lFDRThreshold = groupedFDRCalculator.calcLocalFDRThreshold((float)params.getFdrCometThreshold());

        System.out.print(groupedFDRCalculator.printTree(lFDRThreshold));
        test(cometPsmConverter.getPsms(),groupedFDRCalculator,lFDRThreshold);

        SummaryReportWriter summaryReportWriter = new SummaryReportWriter(params.getOutputDir() +File.separator+ "SummaryReport.txt", params.isIncludeMaxQuant());

        for (String group : groupedFDRCalculator.getGroups()) {

            long start = System.currentTimeMillis();
            ConcurrentHashMap<String, List<PeptideMatchData>> filteredCometPsms = groupedFDRCalculator.filterPsms(cometPsmConverter.getPsms(), lFDRThreshold, group);
            System.out.println("Comet FDR filtering ran in " + RunTime2String.getTimeDiffString(System.currentTimeMillis() - start));


            test(filteredCometPsms,groupedFDRCalculator,lFDRThreshold,group);

            ConcurrentHashMap<String, List<PeptideMatchData>> combined;
            ConcurrentHashMap<String, List<PeptideMatchData>> maxQuantPsms = null;
            if (params.isIncludeMaxQuant()) {
                maxQuantPsms = ProcessPsmUtils.extractGroupPsms(psmGrouper,maxQuantPsmConverter.getPsms(),group);

                start = System.currentTimeMillis();
                combined = combine(filteredCometPsms,maxQuantPsms);
                System.out.println("Comet-MaxQuant combiner ran in " + RunTime2String.getTimeDiffString(System.currentTimeMillis() - start));

            } else {
                combined = filteredCometPsms;
            }

            combined.forEach(10000,spectrumAccumulator);

            writeToCombTabFile(combined, group+"_CometMaxQuantComb.txt");
            System.out.println(combined.size()+" spectra combined for group " + group);

            if (params.isIncludeMaxQuant()) {
                summaryReportWriter.write(group, combined, filteredCometPsms, maxQuantPsms);
            } else {
                summaryReportWriter.write(group, combined);
            }

            System.out.println("Write data to summary report file.");
        }

        summaryReportWriter.close();
    }

    protected void controlFDRGroupwise(GroupedFDRCalculator groupedFDRCalculator, CometPsmConverter cometPsmConverter, MaxQuantPsmConverter maxQuantPsmConverter) {

        Map<String, Float> grpThresholdMap = groupedFDRCalculator.calcGroupLocalFDRThreshold((float)params.getFdrCometThreshold());

        System.out.print(groupedFDRCalculator.printTree(grpThresholdMap));

        SummaryReportWriter summaryReportWriter = new SummaryReportWriter(params.getOutputDir() +File.separator+ "SummaryReport.txt", params.isIncludeMaxQuant());

        for (String group : groupedFDRCalculator.getGroups()) {

            long start = System.currentTimeMillis();
            ConcurrentHashMap<String, List<PeptideMatchData>> filteredCometPsms = groupedFDRCalculator.filterPsms(cometPsmConverter.getPsms(), grpThresholdMap.get(group), group);
            System.out.println("Comet FDR filtering ran in " + RunTime2String.getTimeDiffString(System.currentTimeMillis() - start));

            ConcurrentHashMap<String, List<PeptideMatchData>> combined;
            ConcurrentHashMap<String, List<PeptideMatchData>> maxQuantPsms = null;
            if (params.isIncludeMaxQuant()) {

                maxQuantPsms = ProcessPsmUtils.extractGroupPsms(psmGrouper,maxQuantPsmConverter.getPsms(),group);

                start = System.currentTimeMillis();
                combined = combine(filteredCometPsms,maxQuantPsms);
                System.out.println("Comet-MaxQuant combiner ran in " + RunTime2String.getTimeDiffString(System.currentTimeMillis() - start));

            } else {
                combined = filteredCometPsms;
            }

            combined.forEach(10000,spectrumAccumulator);

            writeToCombTabFile(combined, group+"_CometMaxQuantComb.txt");
            System.out.println(combined.size()+" spectra combined for group " + group);

            if (params.isIncludeMaxQuant()) {
                summaryReportWriter.write(group, combined, filteredCometPsms, maxQuantPsms);
            } else {
                summaryReportWriter.write(group, combined);
            }

            System.out.println("Write data to summary report file.");
        }

        summaryReportWriter.close();
    }

    protected void test(ConcurrentHashMap<String, List<PeptideMatchData>>  cometPsms, GroupedFDRCalculator groupedFDRCalculator, float lFDRThreshold) {

        ConcurrentHashMap<String, List<PeptideMatchData>>  filtered = groupedFDRCalculator.filterPsms(cometPsms, lFDRThreshold);

        int tCnt = 0, dCnt = 0;
        for (String specID : filtered.keySet()) {
            for (PeptideMatchData psm : filtered.get(specID)) {
                if (psm.isDecoy()) dCnt++;
                else tCnt++;
            }
        }

        float[] counts = groupedFDRCalculator.getTargetDecoyCounts(lFDRThreshold);

        System.out.println("root: test Psm count. tCnt= "+tCnt+"/"+counts[1]+", dCnt= "+dCnt+"/"+counts[0]);
    }

    protected void test(ConcurrentHashMap<String, List<PeptideMatchData>>  filtered, GroupedFDRCalculator groupedFDRCalculator, float lFDRThreshold, String group) {

        int tCnt = 0, dCnt = 0;
        for (String specID : filtered.keySet()) {
            for (PeptideMatchData psm : filtered.get(specID)) {
                if (psm.isDecoy()) dCnt++;
                else tCnt++;
            }
        }

        float[] counts = groupedFDRCalculator.getTargetDecoyCounts(lFDRThreshold,group);

        System.out.println(group+": test Psm count. tCnt= "+tCnt+"/"+counts[1]+", dCnt= "+dCnt+"/"+counts[0]);
    }

    protected ConcurrentHashMap<String, List<PeptideMatchData>> combine(ConcurrentHashMap<String, List<PeptideMatchData>>  cometPsms, ConcurrentHashMap<String, List<PeptideMatchData>> maxQuantPsms) {

        ConcurrentHashMap<String, List<PeptideMatchData>> combined = new ConcurrentHashMap<>();

        CometMaxQuantPsmCombiner combiner = new CometMaxQuantPsmCombiner(maxQuantPsms,combined);
        cometPsms.forEach(combiner);

        return combined;
    }

    protected void writeToCombTabFile(ConcurrentHashMap<String, List<PeptideMatchData>>  psms, String filename)  {

        final Psm2StringFunction stringFunction = new Psm2StringFunction("tab");

        StringFileWriter writer = new StringFileWriter(params.getOutputDir() + File.separator +filename, stringFunction);

        psms.forEach(10000, stringFunction, writer);

        writer.close();
    }

    protected void writeToCometTabFile(ConcurrentHashMap<String, List<PeptideMatchData>>  psms, String filename) throws IOException {

        final Psm2StringFunction stringFunction = new Psm2StringFunction("tab", Psm2StringFunction.TabStringMode.COMET);

        StringFileWriter writer = new StringFileWriter(params.getOutputDir() + File.separator +filename, stringFunction);

        psms.forEach(10000, stringFunction, writer);

        writer.close();
    }

    protected void writeToMaxQuantTabFile(ConcurrentHashMap<String, List<PeptideMatchData>>  psms, String filename) throws IOException {

        final Psm2StringFunction stringFunction = new Psm2StringFunction("tab", Psm2StringFunction.TabStringMode.MAXQUANT);

        StringFileWriter writer = new StringFileWriter(params.getOutputDir() + File.separator +filename, stringFunction);

        psms.forEach(10000, stringFunction, writer);

        writer.close();
    }


    protected void writeToFastaFile(ConcurrentHashMap<String, List<PeptideMatchData>>  psms, String filename) throws IOException {

        final Psm2StringFunction stringFunction = new Psm2StringFunction("fasta");

        StringFileWriter writer = new StringFileWriter(params.getOutputDir() + File.separator +filename, stringFunction, true);

        psms.forEach(10000, stringFunction, writer);

        writer.close();
    }


    protected void writeToTextFile(ConcurrentHashMap<String, List<PeptideMatchData>>  psms, String filename) throws IOException {

        final Psm2StringFunction stringFunction = new Psm2StringFunction("txt");

        StringFileWriter writer = new StringFileWriter(params.getOutputDir() + File.separator +filename, stringFunction, true);

        psms.forEach(10000, stringFunction, writer);

        writer.close();
    }

    protected void writePeptideProteinGroupReport(UniProtDB uniProtDB) {

        if (!params.isDoPeptideProteinGrouping()) return;

        if (params.getSearchFastaFile() !=null) uniProtDB.addFastaFile(params.getSearchFastaFile());
        OccamRazorSpectrumCounter spectrumCounter = new OccamRazorSpectrumCounter(spectrumAccumulator, uniProtDB);

        String reportFileName = params.getOutputDir() + File.separator+params.getOutputPrefix()+"PeptideProteinGroupingReport.txt";
        try {
            spectrumCounter.write(new File(reportFileName));
        } catch (IOException e) {
            System.out.println("Cannot writeHistograms protein group report to file "+reportFileName+".");
        }
    }


    protected void createOptions() {

        this.cmdLineOpts = new Options();

        cmdLineOpts.addOption(Option.builder("coD").required().hasArg().longOpt("cometPsmDir").desc("Comet psm root directory (required)").build());
        cmdLineOpts.addOption(Option.builder("mqD").required(false).hasArg().longOpt("maxquantPsmDir").desc("MaxQuant psm root directory. If not provided only Comet is used.").hasArg().build());
        cmdLineOpts.addOption(Option.builder("coRE").required().hasArg().longOpt("cometPsmRegex").desc("Regular expression of Comet psm files (e.g. \\.xml$)").build());
        cmdLineOpts.addOption(Option.builder("coFDR").required().hasArg().longOpt("cometFDR").desc("FDR for filtering Comet PSMs before combination (required)").build());
        cmdLineOpts.addOption(Option.builder("outD").required().hasArg().longOpt("outputDir").desc("Output directory for results (required)").build());
        cmdLineOpts.addOption(Option.builder("repH").required(false).hasArg(false).longOpt("reportHistogram").desc("Report histograms to text files").build());
        cmdLineOpts.addOption(Option.builder("protG").required(false).hasArg().longOpt("proteinGroup").desc("Name of group with protein coding or canonical sequences").build());
        cmdLineOpts.addOption(Option.builder("noncG").required(false).hasArg().longOpt("noncanonicalGroup").desc("Name of group with non-canonical or cryptic sequences").build());
        cmdLineOpts.addOption(Option.builder("protRE").required(false).hasArg().longOpt("protRegExp").desc("Regular expression to match fasta name of coding proteins (e.g. sp\\||tr\\| ").build());
        cmdLineOpts.addOption(Option.builder("spRE").required(false).hasArg().longOpt("spectrumFilter").desc("If this option is set, only spectrum files that match the regexp are used").build());
        cmdLineOpts.addOption(Option.builder("exclP").required(false).hasArg().longOpt("excludeProts").desc("Regular expression of proteins excluded from analysis").build());
        cmdLineOpts.addOption(Option.builder("noncP").required(false).hasArg().longOpt("noncanonicalProts").desc("Comma separated list of protein names to be included in noncanonical group (e.g. PGBD5_HUMAN,POGZ_HUMAN,PGBD1_HUMAN)").build());
        cmdLineOpts.addOption(Option.builder("mod").required(false).hasArg().longOpt("modifications").desc("Comma separated list of additional peptide modifications used in search (e.g. Cysteinyl:C3H5NO2S,Oxidation:O )").build());
        cmdLineOpts.addOption(Option.builder("seFa").required(false).hasArg().longOpt("searchFastaFile").desc("Fasta file used for search (for protein grouping export)").build());
        cmdLineOpts.addOption(Option.builder("upFa").required(false).hasArg().longOpt("uniProtFastaFile").desc("Fasta file with coding or canonical proteins (e.g. UniProt fasta file)").build());
        cmdLineOpts.addOption(Option.builder("ppG").required(false).hasArg(false).longOpt("peptideProteinGrouping").desc("Perform peptide protein grouping export").build());
        cmdLineOpts.addOption(Option.builder("maxR").required(false).hasArg().longOpt("maxRank").desc("Maximal rank of peptide in list of spectrum matches (rank 1 = best) (default value: 1)").build());
        cmdLineOpts.addOption(Option.builder("minZ").required(false).hasArg().longOpt("minCharge").desc("Minimal charge of PSM (default value: 1)").build());
        cmdLineOpts.addOption(Option.builder("maxZ").required(false).hasArg().longOpt("maxCharge").desc("Maximal charge of PSM (default value: 5)").build());
        cmdLineOpts.addOption(Option.builder("minL").required(false).hasArg().longOpt("minLength").desc("Minimal length of peptide (default value: 8)").build());
        cmdLineOpts.addOption(Option.builder("maxL").required(false).hasArg().longOpt("maxLength").desc("Maximal length of peptide (default value: 25)").build());
        cmdLineOpts.addOption(Option.builder("nrTh").required(false).hasArg().longOpt("nrThreads").desc("Number of threads used by NewAnce (default value: nr of available processors - 2)").build());
        cmdLineOpts.addOption(Option.builder("smD").required(false).hasArg().longOpt("smoothDegree").desc("Number of threads used by NewAnce (default value: nr of available processors - 2)").build());
        cmdLineOpts.addOption(Option.builder("outP").required(false).hasArg().longOpt("outputPrefix").desc("Number of smoothing steps. Higher value mean more smoothing (default value: 0).").build());
        cmdLineOpts.addOption(Option.builder("minPH").required(false).hasArg().longOpt("minPsm4Histo").desc("Minimal number of psms to calculate local FDR in histogram (default value: 100000).").build());
        cmdLineOpts.addOption(Option.builder("fdrM").required(false).hasArg().longOpt("fdrControlMethod").desc("Method to control pFDR: global or groupwise (default global).").build());
        cmdLineOpts.addOption(Option.builder("minXC").required(false).hasArg().longOpt("minXCorr").desc("Minimal Comet XCorr in histogram").build());
        cmdLineOpts.addOption(Option.builder("maxXC").required(false).hasArg().longOpt("maxXCorr").desc("Maximal Comet XCorr in histogram").build());
        cmdLineOpts.addOption(Option.builder("nrXCB").required(false).hasArg().longOpt("nrXCorrBins").desc("Number of Comet XCorr bins in histogram").build());
        cmdLineOpts.addOption(Option.builder("minSP").required(false).hasArg().longOpt("minSpScore").desc("Minimal Comet SpScore in histogram").build());
        cmdLineOpts.addOption(Option.builder("maxSP").required(false).hasArg().longOpt("maxSpScore").desc("Maximal Comet SpScore in histogram").build());
        cmdLineOpts.addOption(Option.builder("nrSPB").required(false).hasArg().longOpt("nrSpScoreBins").desc("Number of Comet SpScore bins in histogram").build());
        cmdLineOpts.addOption(Option.builder("minDC").required(false).hasArg().longOpt("minDeltaCn").desc("Minimal Comet DeltaCn in histogram").build());
        cmdLineOpts.addOption(Option.builder("maxDC").required(false).hasArg().longOpt("maxDeltaCn").desc("Maximal Comet DeltaCn in histogram").build());
        cmdLineOpts.addOption(Option.builder("nrDCB").required(false).hasArg().longOpt("nrDeltaCnBins").desc("Number of Comet DeltaCn bins in histogram").build());
        cmdLineOpts.addOption(Option.builder("wP").required(false).hasArg().longOpt("write2ParamFile").desc("Filename where parameters should to written.").build());
        cmdLineOpts.addOption(Option.builder("rP").required(false).hasArg().longOpt("readParamFile").desc("Name of file from which parameters should to read.").build());
        cmdLineOpts.addOption(Option.builder("h").required(false).hasArg(false).longOpt("help").desc("Help option for command line help").build());
        cmdLineOpts.addOption(Option.builder("v").required(false).hasArg(false).longOpt("version").desc("Version of NewAnce software").build());
    }

    @Override
    protected void check(CommandLine line) throws ParseException {

        if (optionsSet) return;

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
        params.add("smoothDegree", getOptionString(line,"smD"));
        params.add("fdrControlMethod", getOptionString(line,"fdrM"));

        params.finalize();
    }
}
