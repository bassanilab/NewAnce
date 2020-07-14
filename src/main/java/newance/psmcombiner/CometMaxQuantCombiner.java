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

import newance.proteinmatch.VariantProtDB;
import newance.psmconverter.PeptideSpectrumMatch;
import newance.util.*;
import org.apache.commons.cli.*;
import newance.proteinmatch.OccamRazorSpectrumCounter;
import newance.proteinmatch.UniProtDB;
import newance.psmconverter.MaxQuantMultipleMSMSFileConverter;
import newance.psmconverter.CometMultiplePepXMLFileConverter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;

/**
 * @author Markus Müller
 */

public class CometMaxQuantCombiner extends ExecutableOptions {

    protected PsmGrouper psmGrouper;
    protected NewAnceParams params;
    protected SpectrumAccumulator spectrumAccumulator;
    protected ConcurrentHashMap<String, List<PeptideSpectrumMatch>> cometPSMs;
    protected GroupedFDRCalculator cometGroupedFDRCalculator;
    protected ConcurrentHashMap<String, List<PeptideSpectrumMatch>> maxQuantPSMs;
    protected GroupedFDRCalculator maxquantGroupedFDRCalculator;
    protected ConcurrentHashMap<String, List<PeptideSpectrumMatch>> combinedPSMs;
    protected Map<String, Float> cometGrplFDRThresholdMap;
    protected Map<String, Float> maxquantGrplFDRThresholdMap;

    public CometMaxQuantCombiner() {

        params = NewAnceParams.getInstance();
        version = params.getVersion();

        createOptions();
        spectrumAccumulator = new SpectrumAccumulator();

        cometPSMs = null;
        maxQuantPSMs = null;
        cometGroupedFDRCalculator = null;
        maxquantGroupedFDRCalculator = null;
        cometGrplFDRThresholdMap = null;
        maxquantGrplFDRThresholdMap = null;
    }

    public static void main(String[] args) {

        CometMaxQuantCombiner cometScoreCombiner =  new CometMaxQuantCombiner();
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

        long start = System.currentTimeMillis();

        writeNewAnceInfo();
        if (params.isWriteParamsFile()) writeParams();

        UniProtDB uniProtDB = null;
        if (!params.getUniprotFastaFile().isEmpty()) {
            System.out.println("Load fasta sequences from "+params.getUniprotFastaFile()+" ...");
            uniProtDB = new UniProtDB(params.getUniprotFastaFile());
        }

        VariantProtDB variantProtDB = null;
        if (!params.getSearchFastaFile().isEmpty()) {
            System.out.println("Load fasta sequences from "+params.getSearchFastaFile()+" ...");
            variantProtDB = new VariantProtDB(params.getSearchFastaFile());
        }

        processPSMs(uniProtDB, variantProtDB, NewAnceParams.SearchTool.COMET);
        processPSMs(uniProtDB, variantProtDB, NewAnceParams.SearchTool.MAXQUANT);

        if (params.getFdrControlMethod().equals("combined")) {
            controlFDRCombined();
        } else {
            controlFDRSeparate();
        }

        writePeptideProteinGroupReport(uniProtDB);

        System.out.println("NewAnce ran in " +
                RunTime2String.getTimeDiffString(System.currentTimeMillis() - start));
        return 0;
    }


    protected void controlFDRCombined() {

        float cometlFDRThreshold = cometGroupedFDRCalculator.calcLocalFDRThreshold((float)params.getFdrCometThreshold());
        float maxquantlFDRThreshold = maxquantGroupedFDRCalculator.calcLocalFDRThreshold((float)params.getFdrMaxQuantThreshold());

        cometGrplFDRThresholdMap = new HashMap<>();
        maxquantGrplFDRThresholdMap = new HashMap<>();

        for (String grp : cometGroupedFDRCalculator.getGroups()) {
            cometGrplFDRThresholdMap.put(grp,cometlFDRThreshold);
            maxquantGrplFDRThresholdMap.put(grp,maxquantlFDRThreshold);
        }

        reportPSMs();

        String sumFileName = params.getOutputDir() +File.separator+params.getOutputTag()+"_SummaryReport.txt";
        SummaryReportWriter summaryReportWriter = new SummaryReportWriter(sumFileName, true);

        for (String group : cometGroupedFDRCalculator.getGroups()) {

            ConcurrentHashMap<String, List<PeptideSpectrumMatch>> filteredCometPsms =
                    cometGroupedFDRCalculator.filterPsms(cometPSMs, cometlFDRThreshold, group);

            ConcurrentHashMap<String, List<PeptideSpectrumMatch>> filteredMaxQuantPsms =
                    maxquantGroupedFDRCalculator.filterPsms(maxQuantPSMs, maxquantlFDRThreshold, group);


            combinedPSMs = combine(filteredCometPsms,filteredMaxQuantPsms);

            combinedPSMs.forEach(10000,spectrumAccumulator);

            writeToCombTabFile(group+"_"+params.getOutputTag()+"_NewAncePSMs.txt");
            System.out.println(combinedPSMs.size()+" spectra combined for group " + group);

            summaryReportWriter.write(group, combinedPSMs, filteredCometPsms, filteredMaxQuantPsms);

            System.out.println("Write data to summary report file.");
        }

        summaryReportWriter.close();
    }

    protected void controlFDRSeparate() {

        cometGrplFDRThresholdMap =
                cometGroupedFDRCalculator.calcGroupLocalFDRThreshold((float)params.getFdrCometThreshold());
        maxquantGrplFDRThresholdMap =
                maxquantGroupedFDRCalculator.calcGroupLocalFDRThreshold((float)params.getFdrMaxQuantThreshold());

        reportPSMs();

        String sumFileName = params.getOutputDir() +File.separator+params.getOutputTag()+"_SummaryReport.txt";
        SummaryReportWriter summaryReportWriter = new SummaryReportWriter(sumFileName, true);

        for (String group : cometGroupedFDRCalculator.getGroups()) {

            ConcurrentHashMap<String, List<PeptideSpectrumMatch>> filteredCometPsms =
                    cometGroupedFDRCalculator.filterPsms(cometPSMs, cometGrplFDRThresholdMap.get(group), group);

            ConcurrentHashMap<String, List<PeptideSpectrumMatch>> filteredMaxQuantPsms =
                    maxquantGroupedFDRCalculator.filterPsms(maxQuantPSMs, maxquantGrplFDRThresholdMap.get(group), group);


            combinedPSMs = combine(filteredCometPsms,filteredMaxQuantPsms);

            combinedPSMs.forEach(10000,spectrumAccumulator);

            writeToCombTabFile(group+"_"+params.getOutputTag()+"_NewAncePSMs.txt");
            System.out.println(combinedPSMs.size()+" spectra combined for group " + group);

            summaryReportWriter.write(group, combinedPSMs, filteredCometPsms, filteredMaxQuantPsms);

            System.out.println("Write data to summary report file.");
        }

        summaryReportWriter.close();

    }

    protected void processPSMs(UniProtDB uniProtDB, VariantProtDB variantProtDB,
                               NewAnceParams.SearchTool searchTool) throws IOException {

        GroupedFDRCalculator groupedFDRCalculator =
                GroupedFDRCalculator.buildGroupedFDRCalculator(uniProtDB, variantProtDB, searchTool);

        if (searchTool == NewAnceParams.SearchTool.COMET) {

            cometGroupedFDRCalculator = groupedFDRCalculator;
            CometMultiplePepXMLFileConverter cometMultiplePepXMLFileConverter =
                    new CometMultiplePepXMLFileConverter(params.getCometPsmDir(), params.getCometPsmRegExp(),
                            groupedFDRCalculator, false);

            cometMultiplePepXMLFileConverter.run();

            cometPSMs = cometMultiplePepXMLFileConverter.getPsms();
        } else {

            maxquantGroupedFDRCalculator = groupedFDRCalculator;
            MaxQuantMultipleMSMSFileConverter maxQuantMultipleMSMSConverter =
                    new MaxQuantMultipleMSMSFileConverter(params.getMaxquantPsmDir(), params.getMaxquantPsmRegExp(),
                            groupedFDRCalculator, false);

            maxQuantMultipleMSMSConverter.run();

            maxQuantPSMs = maxQuantMultipleMSMSConverter.getPsms();
        }

        groupedFDRCalculator.process(params.getMinNrPsmsPerHisto(), params.getSmoothDegree()); // calculate local fdr here

        String dir = (searchTool == NewAnceParams.SearchTool.COMET)?"cometHistos":"maxquantHistos";
        if (params.isReportHistos())
            groupedFDRCalculator.writeHistograms(params.getOutputDir()+File.separator+dir, params.getOutputTag());

    }

    protected void reportPSMs() {

        if (params.reportAllPSM()) {
            writePSMTabFile(cometPSMs, cometGroupedFDRCalculator, cometGrplFDRThresholdMap,
                    params.getOutputTag() + "_CometPSMs.txt", NewAnceParams.SearchTool.COMET);

            writePSMTabFile(maxQuantPSMs, maxquantGroupedFDRCalculator, maxquantGrplFDRThresholdMap,
                    params.getOutputTag() + "_MaxQuantPSMs.txt", NewAnceParams.SearchTool.MAXQUANT);
        }

        writeGroupHistoTree(cometGroupedFDRCalculator.printTree(cometGrplFDRThresholdMap), false);
        writeGroupHistoTree(cometGroupedFDRCalculator.printTree(maxquantGrplFDRThresholdMap), true);
    }

    protected void writeGroupHistoTree(String grpStatisticsInfo, boolean append) {

        String paramsFileName = params.getOutputDir() + File.separator + params.getOutputTag()+"_NewAnceStatistics.txt";

        try {
            BufferedWriter paramsWriter = new BufferedWriter(new FileWriter(new File(paramsFileName), append));
            paramsWriter.write(grpStatisticsInfo);
            paramsWriter.close();
        } catch (IOException e) {
            System.out.println("Cannot write NewAnce parameters to file "+paramsFileName+".");
        }
    }

    protected ConcurrentHashMap<String, List<PeptideSpectrumMatch>> combine(
            ConcurrentHashMap<String, List<PeptideSpectrumMatch>>  cometPsms,
            ConcurrentHashMap<String, List<PeptideSpectrumMatch>> maxQuantPsms) {

        if (maxQuantPsms == null) return cometPsms;

        ConcurrentHashMap<String, List<PeptideSpectrumMatch>> combined = new ConcurrentHashMap<>();

        CometMaxQuantPsmMerger combiner = new CometMaxQuantPsmMerger(maxQuantPsms,combined);
        cometPsms.forEach(combiner);

        return combined;
    }

    protected void writeToCombTabFile(String filename)  {

        final CombinedPsm2StringFunction stringFunction =
                new CombinedPsm2StringFunction(cometGroupedFDRCalculator, maxquantGroupedFDRCalculator,
                        cometGrplFDRThresholdMap, maxquantGrplFDRThresholdMap);

        StringFileWriter writer =
                new StringFileWriter(params.getOutputDir() + File.separator +filename, stringFunction);

        combinedPSMs.forEach(10000, stringFunction, writer);

        writer.close();
    }


    protected void writePSMTabFile(ConcurrentHashMap<String, List<PeptideSpectrumMatch>>  psms,
                                   GroupedFDRCalculator groupedFDRCalculator, Map<String, Float> grpThresholdMap,
                                   String filename, NewAnceParams.SearchTool searchTool)  {

        Psm2StringFunction stringFunction = null;

        if (searchTool == NewAnceParams.SearchTool.COMET)
            stringFunction = new CometPsm2StringFunction(groupedFDRCalculator,grpThresholdMap);
        else
            stringFunction = new MaxQuantPsm2StringFunction(groupedFDRCalculator,grpThresholdMap);

        StringFileWriter writer =
                new StringFileWriter(params.getOutputDir() + File.separator +filename, stringFunction);

        psms.forEach(10000, stringFunction, writer);

        writer.close();
    }


    protected void writePeptideProteinGroupReport(UniProtDB uniProtDB) {

        if (!params.isDoPeptideProteinGrouping()) return;

        if (params.getSearchFastaFile() !=null) uniProtDB.addFastaFile(params.getSearchFastaFile());
        OccamRazorSpectrumCounter spectrumCounter = new OccamRazorSpectrumCounter(spectrumAccumulator, uniProtDB);

        String reportFileName = params.getOutputDir() + File.separator + params.getOutputTag()+"_PeptideProteinGroupingReport.txt";
        try {
            spectrumCounter.write(new File(reportFileName));
        } catch (IOException e) {
            System.out.println("Cannot writeHistograms protein group report to file "+reportFileName+".");
        }
    }

    protected void writeParams() {

        String paramsFileName = params.getOutputDir() + File.separator + params.getOutputTag()+"_NewAnceParameters.txt";

        try {
            BufferedWriter paramsWriter = new BufferedWriter(new FileWriter(new File(paramsFileName)));
            paramsWriter.write(params.toString());
            paramsWriter.close();
        } catch (IOException e) {
            System.out.println("Cannot write NewAnce parameters to file "+paramsFileName+".");
        }
    }

    protected void writeNewAnceInfo() {

        String sepLine = "*************************************************************************************************************************";
        String newanceInfo = "Running NewAnce Version "+params.getVersion();
        int gap = (sepLine.length()-4-newanceInfo.length())/2;
        String gapStr = "**";
        for (int i=0;i<gap;i++) gapStr += " ";
        newanceInfo = gapStr + newanceInfo;
        for (int i=newanceInfo.length();i<sepLine.length()-2;i++) newanceInfo += " ";
        newanceInfo += "**";
        System.out.println("");
        System.out.println(sepLine);
        System.out.println(newanceInfo);
        System.out.println(sepLine);
        System.out.println("");
        System.out.println("Parameters:");
        System.out.println(params.toString());
        System.out.println("");
    }

    protected void createOptions() {

        this.cmdLineOpts = new Options();

        cmdLineOpts.addOption(Option.builder("coD").required().hasArg().longOpt("cometPsmDir").desc("Comet psm root directory (required)").build());
        cmdLineOpts.addOption(Option.builder("mqD").required().hasArg().longOpt("maxquantPsmDir").desc("MaxQuant psm root directory.").hasArg().build());
        cmdLineOpts.addOption(Option.builder("coRE").required().hasArg().longOpt("cometPsmRegex").desc("Regular expression of Comet psm files (e.g. \\.xml$) (required)").build());
        cmdLineOpts.addOption(Option.builder("coFDR").required(false).hasArg().longOpt("cometFDR").desc("FDR for filtering Comet PSMs before combination (default value 0.03)").build());
        cmdLineOpts.addOption(Option.builder("mqFDR").required(false).hasArg().longOpt("maxquantFDR").desc("FDR for filtering MaxQuant PSMs before combination (default value 0.03)").build());
        cmdLineOpts.addOption(Option.builder("outD").required().hasArg().longOpt("outputDir").desc("Output directory for results (required)").build());
        cmdLineOpts.addOption(Option.builder("repH").required(false).hasArg(false).longOpt("reportHistogram").desc("Report histograms to text files").build());
        cmdLineOpts.addOption(Option.builder("rCoH").required(false).hasArg().longOpt("readCometHistograms").desc("Directory where Comet histograms files are placed.").build());
        cmdLineOpts.addOption(Option.builder("rMqH").required(false).hasArg().longOpt("readMaxQuantHistograms").desc("Directory where MaxQuant histograms files are placed.").build());
        cmdLineOpts.addOption(Option.builder("fH").required(false).hasArg(false).longOpt("forceHistograms").desc("Histograms are imported even if enough PSMs are available.").build());
        cmdLineOpts.addOption(Option.builder("groupM").required(false).hasArg().longOpt("groupingMethod").desc("Method for PSM grouping: fasta or modif or none (default none).").build());
        cmdLineOpts.addOption(Option.builder("groupN").required(false).hasArg().longOpt("groupNames").desc("Comma separated list of names of sequence groups in fasta file (e.g. prot,lncRNA,TE ). Will be used as prefixes for output files.").build());
        cmdLineOpts.addOption(Option.builder("groupRE").required(false).hasArg().longOpt("groupRegEx").desc("Comma separated list of regular expression defining sequence groups of fasta headers (e.g. \"sp\\||tr\\|ENSP00\",\"ENST00\",\"SINE_|LINE_|LTR_|DNA_|Retroposon_\" ). Will be used as prefixes for output files.").build());
        cmdLineOpts.addOption(Option.builder("wAll").required(false).hasArg(false).longOpt("writeFullPSMExport").desc("If flag is set, all Comet and MaxQuant PSMs are written to a tab file.").build());
        cmdLineOpts.addOption(Option.builder("spRE").required(false).hasArg().longOpt("spectrumFilter").desc("If this option is set, only spectrum ids that match this regexp are used.  If not set no filtering is performed.").build());
        cmdLineOpts.addOption(Option.builder("exclP").required(false).hasArg().longOpt("excludeProts").desc("Regular expression of proteins excluded from analysis. If not set no proteins are excluded.").build());
        cmdLineOpts.addOption(Option.builder("groupF").required(false).hasArg().longOpt("groupProteinFile").desc("Tab file with protein group assignments which will override assignment by groupRE").build());
        cmdLineOpts.addOption(Option.builder("mod").required(false).hasArg().longOpt("modifications").desc("Comma separated list of peptide modifications used in search (e.g. Cysteinyl:C3H5NO2S,Oxidation:O)").build());
        cmdLineOpts.addOption(Option.builder("seFa").required(false).hasArg().longOpt("searchFastaFile").desc("Fasta file that was used for the search (required for protein grouping export and annotation of variants)").build());
        cmdLineOpts.addOption(Option.builder("upFa").required(false).hasArg().longOpt("uniProtFastaFile").desc("Fasta file with coding or canonical proteins (e.g. UniProt fasta file)").build());
        cmdLineOpts.addOption(Option.builder("ppG").required(false).hasArg(false).longOpt("peptideProteinGrouping").desc("Perform peptide protein grouping export.").build());
        cmdLineOpts.addOption(Option.builder("maxR").required(false).hasArg().longOpt("maxRank").desc("Maximal rank of peptide in list of spectrum matches (rank 1 = best) (default value: 1)").build());
        cmdLineOpts.addOption(Option.builder("minZ").required(false).hasArg().longOpt("minCharge").desc("Minimal charge of PSM (default value: 1)").build());
        cmdLineOpts.addOption(Option.builder("maxZ").required(false).hasArg().longOpt("maxCharge").desc("Maximal charge of PSM (default value: 5)").build());
        cmdLineOpts.addOption(Option.builder("minL").required(false).hasArg().longOpt("minLength").desc("Minimal length of peptide (default value: 8)").build());
        cmdLineOpts.addOption(Option.builder("maxL").required(false).hasArg().longOpt("maxLength").desc("Maximal length of peptide (default value: 25)").build());
        cmdLineOpts.addOption(Option.builder("nrTh").required(false).hasArg().longOpt("nrThreads").desc("Number of threads used by NewAnce (default value: nr of available processors - 2)").build());
        cmdLineOpts.addOption(Option.builder("smD").required(false).hasArg().longOpt("smoothDegree").desc("Degree of smoothing (0: no smoothing, n: n x smoothing) (default value 1)").build());
        cmdLineOpts.addOption(Option.builder("outT").required(false).hasArg().longOpt("outputTag").desc("Tag inserted into output file names after prefix.").build());
        cmdLineOpts.addOption(Option.builder("minPH").required(false).hasArg().longOpt("minPsm4Histo").desc("Minimal number of psms to calculate local FDR in histogram (default value: 100000).").build());
        cmdLineOpts.addOption(Option.builder("fdrM").required(false).hasArg().longOpt("fdrControlMethod").desc("Method to control pFDR: combined or separate (default combined).").build());
        cmdLineOpts.addOption(Option.builder("minXC").required(false).hasArg().longOpt("minXCorr").desc("Minimal Comet XCorr in histogram (default value 0)").build());
        cmdLineOpts.addOption(Option.builder("maxXC").required(false).hasArg().longOpt("maxXCorr").desc("Maximal Comet XCorr in histogram (default value 5)").build());
        cmdLineOpts.addOption(Option.builder("nrXCB").required(false).hasArg().longOpt("nrXCorrBins").desc("Number of Comet XCorr bins in histogram (default value 40)" ).build());
        cmdLineOpts.addOption(Option.builder("minSP").required(false).hasArg().longOpt("minSpScore").desc("Minimal Comet SpScore in histogram (default value 0)").build());
        cmdLineOpts.addOption(Option.builder("maxSP").required(false).hasArg().longOpt("maxSpScore").desc("Maximal Comet SpScore in histogram (default value 1)").build());
        cmdLineOpts.addOption(Option.builder("nrSPB").required(false).hasArg().longOpt("nrSpScoreBins").desc("Number of Comet SpScore bins in histogram (default value 40)").build());
        cmdLineOpts.addOption(Option.builder("minDC").required(false).hasArg().longOpt("minDeltaCn").desc("Minimal Comet DeltaCn in histogram (default value 0)").build());
        cmdLineOpts.addOption(Option.builder("maxDC").required(false).hasArg().longOpt("maxDeltaCn").desc("Maximal Comet DeltaCn in histogram (default value 2500)").build());
        cmdLineOpts.addOption(Option.builder("nrDCB").required(false).hasArg().longOpt("nrDeltaCnBins").desc("Number of Comet DeltaCn bins in histogram (default value 40)").build());
        cmdLineOpts.addOption(Option.builder("minScore").required(false).hasArg().longOpt("minScore").desc("Minimal MaxQuant Score in histogram (default value 0)").build());
        cmdLineOpts.addOption(Option.builder("maxScore").required(false).hasArg().longOpt("maxScore").desc("Maximal MaxQuant Score in histogram (default value 5)").build());
        cmdLineOpts.addOption(Option.builder("nrScoreB").required(false).hasArg().longOpt("nrScoreBins").desc("Number of MaxQuant Score bins in histogram (default value 40)" ).build());
        cmdLineOpts.addOption(Option.builder("minDS").required(false).hasArg().longOpt("minDeltaScore").desc("Minimal MaxQuant DeltaScore in histogram (default value 0)").build());
        cmdLineOpts.addOption(Option.builder("maxDS").required(false).hasArg().longOpt("maxDeltaScore").desc("Maximal MaxQuant DeltaScore in histogram (default value 1)").build());
        cmdLineOpts.addOption(Option.builder("nrDSB").required(false).hasArg().longOpt("nrDeltaScoreBins").desc("Number of MaxQuant DeltaScore bins in histogram (default value 40)").build());
        cmdLineOpts.addOption(Option.builder("minPEP").required(false).hasArg().longOpt("minPEP").desc("Minimal MaxQuant PEP in histogram (default value 0)").build());
        cmdLineOpts.addOption(Option.builder("maxPEP").required(false).hasArg().longOpt("maxPEP").desc("Maximal MaxQuant PEP in histogram (default value 2500)").build());
        cmdLineOpts.addOption(Option.builder("nrPEPB").required(false).hasArg().longOpt("nrPEPBins").desc("Number of MaxQuant PEP bins in histogram (default value 40)").build());
        cmdLineOpts.addOption(Option.builder("wP").required(false).hasArg(false).longOpt("write2ParamFile").desc("This option is set if parameters should be written to file.").build());
        cmdLineOpts.addOption(Option.builder("rP").required(false).hasArg().longOpt("readParamFile").desc("Name of file from which parameters should to read.").build());
        cmdLineOpts.addOption(Option.builder("d").required(false).hasArg(false).longOpt("debug").desc("Debug option").build());
        cmdLineOpts.addOption(Option.builder("h").required(false).hasArg(false).longOpt("help").desc("Help option for command line help").build());
        cmdLineOpts.addOption(Option.builder("v").required(false).hasArg(false).longOpt("version").desc("Version of NewAnce software").build());
    }

    @Override
    protected void check(CommandLine line) throws ParseException {

        if (optionsSet) return;

        this.params = NewAnceParams.getInstance();

        params.add("cometPsmDir", getOptionString(line, "coD"));
        params.add("cometPsmRegExp", getOptionString(line, "coRE"));

        String mqDir = getOptionString(line, "mqD");
        params.add("maxquantPsmDir", getOptionString(line, "mqD"));
        params.add("debug", getOptionString(line, "d"));
        params.add("forceHistos", getOptionString(line, "fH"));
        params.add("reportHistos", getOptionString(line, "repH"));
        params.add("readCometHistos", getOptionString(line, "rCoH"));
        params.add("readMaxQuantHistos", getOptionString(line, "rMqH"));
        params.add("outputDir", getOptionString(line, "outD"));
        params.add("searchFastaFile", getOptionString(line, "seFa"));
        params.add("uniprotFastaFile", getOptionString(line, "upFa"));
        params.add("doPeptideProteinGrouping", getOptionString(line, "ppG"));
        params.add("writeParamsFile", getOptionString(line, "wP"));
        params.add("readParamsFile", getOptionString(line, "rP"));
        params.add("maxRank", getOptionString(line, "maxR"));
        params.add("minCharge", getOptionString(line, "minZ"));
        params.add("maxCharge", getOptionString(line, "maxZ"));
        params.add("minPeptideLength", getOptionString(line, "minL"));
        params.add("maxPeptideLength", getOptionString(line, "maxL"));
        params.add("fdrCometThreshold", getOptionString(line, "coFDR"));
        params.add("excludedProtPattern", getOptionString(line, "exclP"));
        params.add("proteinGroupMapFile", getOptionString(line, "groupF"));
        params.add("spectrumRegExp", getOptionString(line, "spRE"));
        params.add("outputTag", getOptionString(line, "outT"));
        params.add("modifications", getOptionString(line, "mod"));
        params.add("minNrPsmsPerHisto", getOptionString(line, "minPH"));
        params.add("minXCorr", getOptionString(line, "minXC"));
        params.add("maxXCorr", getOptionString(line, "maxXC"));
        params.add("nrXCorrBins", getOptionString(line, "nrXCB"));
        params.add("minDeltaCn", getOptionString(line, "minDC"));
        params.add("maxDeltaCn", getOptionString(line, "maxDC"));
        params.add("nrDeltaCnBins", getOptionString(line, "nrDCB"));
        params.add("minSpScore", getOptionString(line, "minSP"));
        params.add("maxSpScore", getOptionString(line, "maxSP"));
        params.add("nrSpScoreBins", getOptionString(line, "nrSPB"));
        params.add("minScore", getOptionString(line, "minScore"));
        params.add("maxScore", getOptionString(line, "maxScore"));
        params.add("nrScoreBins", getOptionString(line, "nrScoreB"));
        params.add("minDeltaScore", getOptionString(line, "minDS"));
        params.add("maxDeltaScore", getOptionString(line, "maxDS"));
        params.add("nrDeltaScoreBins", getOptionString(line, "nrDSB"));
        params.add("minPEP", getOptionString(line, "minPEP"));
        params.add("maxPEP", getOptionString(line, "maxPEP"));
        params.add("nrPEPBins", getOptionString(line, "nrPEPB"));
        params.add("nrThreads", getOptionString(line, "nrTh"));
        params.add("smoothDegree", getOptionString(line, "smD"));
        params.add("fdrControlMethod", getOptionString(line, "fdrM"));
        params.add("groupingMethod", getOptionString(line, "groupM"));
        params.add("groupNames", getOptionString(line, "groupN"));
        params.add("groupRegExs", getOptionString(line, "groupRE"));
        params.add("reportAllPSM", getOptionString(line, "wAll"));

        params.finalize();

    }
}
