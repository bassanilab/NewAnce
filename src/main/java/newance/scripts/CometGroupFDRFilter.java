/*
Copyright (C) SIB - Swiss Institute of Bioinformatics, Lausanne, Switzerland
Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
*/

package newance.scripts;

import newance.proteinmatch.OccamRazorSpectrumCounter;
import newance.proteinmatch.UniProtDB;
import newance.proteinmatch.VariantProtDB;
import newance.psmcombiner.*;
import newance.psmconverter.CometMultiplePepXMLFileConverter;
import newance.psmconverter.PeptideSpectrumMatch;
import newance.util.*;
import org.apache.commons.cli.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * @author Markus Müller
 */

public class CometGroupFDRFilter extends ExecutableOptions {

    protected NewAnceParams params;
    protected SpectrumAccumulator spectrumAccumulator;
    protected ConcurrentHashMap<String, List<PeptideSpectrumMatch>> cometPSMs;
    protected GroupedFDRCalculator groupedFDRCalculator;
    protected Map<String, Float> grplFDRThresholdMap;

    public CometGroupFDRFilter() {

        params = NewAnceParams.getInstance();
        version = params.getVersion();

        createOptions();
        spectrumAccumulator = new SpectrumAccumulator();

        cometPSMs = null;
        groupedFDRCalculator = null;
        grplFDRThresholdMap = null;
    }

    public static void main(String[] args) {

        CometGroupFDRFilter cometGroupFDRFilter =  new CometGroupFDRFilter();
        try {
            cometGroupFDRFilter.init(args).parseOptions(args).run();
        } catch (MissingOptionException e) {
        } catch (ParseException e) {
            cometGroupFDRFilter.printOptions(args, e.getMessage());
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

        processPSMs(uniProtDB, variantProtDB);

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

        grplFDRThresholdMap = new HashMap<>();
        float lFDRThreshold = groupedFDRCalculator.calcLocalFDRThreshold((float)params.getFdrCometThreshold());
        for (String grp : groupedFDRCalculator.getGroups()) grplFDRThresholdMap.put(grp,lFDRThreshold);

        reportPSMs();

        String sumFileName = params.getOutputDir() +File.separator+params.getOutputTag()+"_SummaryReport.txt";
        SummaryReportWriter summaryReportWriter =
                new SummaryReportWriter(sumFileName, false);

        for (String group : groupedFDRCalculator.getGroups()) {

            ConcurrentHashMap<String, List<PeptideSpectrumMatch>> filteredCometPsms =
                    groupedFDRCalculator.filterPsms(cometPSMs, lFDRThreshold, group);

            filteredCometPsms.forEach(10000,spectrumAccumulator);

            String filename = group+"_"+params.getOutputTag()+"_NewAncePSMs.txt";
            writePSMTabFile(filteredCometPsms, groupedFDRCalculator, grplFDRThresholdMap, filename);
            System.out.println(filteredCometPsms.size()+" spectra for group " + group);

            summaryReportWriter.write(group, filteredCometPsms);
        }

        summaryReportWriter.close();
    }

    protected void controlFDRSeparate() {

        grplFDRThresholdMap = groupedFDRCalculator.calcGroupLocalFDRThreshold((float)params.getFdrCometThreshold());

        reportPSMs();

        String sumFileName = params.getOutputDir() +File.separator+params.getOutputTag()+"_SummaryReport.txt";
        SummaryReportWriter summaryReportWriter =
                new SummaryReportWriter(sumFileName, false);


        for (String group : groupedFDRCalculator.getGroups()) {

            ConcurrentHashMap<String, List<PeptideSpectrumMatch>> filteredCometPsms =
                    groupedFDRCalculator.filterPsms(cometPSMs, grplFDRThresholdMap.get(group), group);

            filteredCometPsms.forEach(10000,spectrumAccumulator);

            String filename = group+"_"+params.getOutputTag()+"_NewAncePSMs.txt";
            writePSMTabFile(filteredCometPsms, groupedFDRCalculator, grplFDRThresholdMap, filename);
            System.out.println(filteredCometPsms.size()+" spectra combined for group " + group);

            summaryReportWriter.write(group, filteredCometPsms);
        }

        summaryReportWriter.close();
    }

    protected void processPSMs(UniProtDB uniProtDB, VariantProtDB variantProtDB) throws IOException {

        System.out.println("Process PSMs ...");
        groupedFDRCalculator =
                GroupedFDRCalculator.buildGroupedFDRCalculator(uniProtDB, variantProtDB, NewAnceParams.SearchTool.COMET);

        CometMultiplePepXMLFileConverter cometMultiplePepXMLFileConverter =
                new CometMultiplePepXMLFileConverter(params.getCometPsmDir(), params.getCometPsmRegExp(),
                        groupedFDRCalculator, false);
        cometMultiplePepXMLFileConverter.run();

        cometPSMs = cometMultiplePepXMLFileConverter.getPsms();

        groupedFDRCalculator.process(params.getMinNrPsmsPerHisto(), params.getSmoothDegree()); // calculate local fdr here

        if (params.isReportHistos())
            groupedFDRCalculator.writeHistograms(params.getOutputDir()+File.separator+"cometHistos", params.getOutputTag());

    }

    protected void reportPSMs() {

        if (params.reportAllPSM()) {
            writePSMTabFile(cometPSMs, groupedFDRCalculator, grplFDRThresholdMap,
                    params.getOutputTag() + "_CometPSMs.txt");
        }

        writeGroupHistoTree(groupedFDRCalculator.printTree(grplFDRThresholdMap), false);
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

    protected void writePSMTabFile(ConcurrentHashMap<String, List<PeptideSpectrumMatch>>  psms,
                                   GroupedFDRCalculator groupedFDRCalculator,
                                   Map<String, Float> grpThresholdMap, String filename)  {

        final Psm2StringFunction stringFunction;

        if (params.isOutputPsql()){
            stringFunction = new CometPsm2PsqlStringFunction(groupedFDRCalculator,grpThresholdMap);
        } else{
            stringFunction = new CometPsm2StringFunction(groupedFDRCalculator,grpThresholdMap);
        }

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
        String newanceInfo = "Running NewAnce (Comet only) Version "+params.getVersion();
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
        cmdLineOpts.addOption(Option.builder("coRE").required().hasArg().longOpt("cometPsmRegex").desc("Regular expression of Comet psm files (e.g. \\.xml$) (required)").build());
        cmdLineOpts.addOption(Option.builder("coFDR").required(false).hasArg().longOpt("cometFDR").desc("FDR for filtering Comet PSMs before combination (default value 0.03)").build());
        cmdLineOpts.addOption(Option.builder("outD").required().hasArg().longOpt("outputDir").desc("Output directory for results (required)").build());
        cmdLineOpts.addOption(Option.builder("outPsql").required(false).hasArg(false).longOpt("outputPsql").desc("Output format (table) that supports direct import into PostgreSQL.").build());
        cmdLineOpts.addOption(Option.builder("repH").required(false).hasArg(false).longOpt("reportHistogram").desc("Report histograms to text files").build());
        cmdLineOpts.addOption(Option.builder("rCoH").required(false).hasArg().longOpt("readCometHistograms").desc("Directory where Comet histograms files are placed.").build());
        cmdLineOpts.addOption(Option.builder("alpha").required(false).hasArg().longOpt("combineHistoWeight").desc("Histograms are alpha*data_histo + (1-alpha)*prior_histo. 0 <0 alpha <= 1. Only used if rCoH option is set.").build());
        cmdLineOpts.addOption(Option.builder("fH").required(false).hasArg(false).longOpt("forceHistograms").desc("Histograms are imported even if enough PSMs are available.").build());
        cmdLineOpts.addOption(Option.builder("groupM").required(false).hasArg().longOpt("groupingMethod").desc("Method for PSM grouping: fasta, modif, famo (fasta&modif) or none (default none).").build());
        cmdLineOpts.addOption(Option.builder("groupN").required(false).hasArg().longOpt("groupNames").desc("Comma separated list of names of sequence groups in fasta file (e.g. prot,lncRNA,TE ). Will be used as prefixes for output files.").build());
        cmdLineOpts.addOption(Option.builder("groupRE").required(false).hasArg().longOpt("groupRegEx").desc("Comma separated list of regular expression defining sequence groups of fasta headers (e.g. \"sp\\||tr\\|ENSP00\",\"ENST00\",\"SINE_|LINE_|LTR_|DNA_|Retroposon_\" ). Will be used as prefixes for output files.").build());
        cmdLineOpts.addOption(Option.builder("wAll").required(false).hasArg(false).longOpt("writeFullCometExport").desc("If flag is set, all Comet PSMs are written to a tab file.").build());
        cmdLineOpts.addOption(Option.builder("spRE").required(false).hasArg().longOpt("spectrumFilter").desc("If this option is set, only spectrum ids that match this regexp are used.  If not set no filtering is performed.").build());
        cmdLineOpts.addOption(Option.builder("exclP").required(false).hasArg().longOpt("excludeProts").desc("Regular expression of proteins excluded from analysis. If not set no proteins are excluded.").build());
        cmdLineOpts.addOption(Option.builder("groupF").required(false).hasArg().longOpt("groupProteinFile").desc("Tab file with protein group assignments which will override assignment by groupRE").build());
        cmdLineOpts.addOption(Option.builder("mod").required(false).hasArg().longOpt("modifications").desc("Comma separated list of peptide modifications used in search (e.g. Cysteinyl:C3H5NO2S,Oxidation:O)").build());
        cmdLineOpts.addOption(Option.builder("seFa").required(false).hasArg().longOpt("searchFastaFile").desc("Fasta file that was used for the search (required for protein grouping export and annotation of variants in the Comet results)").build());
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
        cmdLineOpts.addOption(Option.builder("minXCPSM").required(false).hasArg().longOpt("minXCorrPSM").desc("Minimal Comet xcorr for PSM (default value 0.8)").build());
        cmdLineOpts.addOption(Option.builder("minSPPSM").required(false).hasArg().longOpt("minSpScorePSM").desc("Minimal Comet spscore for PSM (default value 50)").build());
        cmdLineOpts.addOption(Option.builder("minDCPSM").required(false).hasArg().longOpt("minDeltaCnPSM").desc("Minimal Comet deltacn for PSM (default value 0)").build());
        cmdLineOpts.addOption(Option.builder("minScorePSM").required(false).hasArg().longOpt("minScorePSM").desc("Minimal MaxQuant score for PSM (default value 10)").build());
        cmdLineOpts.addOption(Option.builder("minDSPSM").required(false).hasArg().longOpt("minDeltaScorePSM").desc("Minimal MaxQuant DeltaScore for PSM (default value 0)").build());
        cmdLineOpts.addOption(Option.builder("minPEPPSM").required(false).hasArg().longOpt("minPEPPSM").desc("Minimal MaxQuant PEP for PSM (default value 0)").build());
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

        params.add("debug", getOptionString(line, "d"));
        params.add("forceHistos", getOptionString(line, "fH"));
        params.add("reportHistos", getOptionString(line, "repH"));
        params.add("readCometHistos", getOptionString(line, "rCoH"));
        params.add("alpha", getOptionString(line, "alpha"));
        params.add("outputDir", getOptionString(line, "outD"));
        params.add("outputPsql", getOptionString(line, "outPsql"));
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
        params.add("minXCorrPSM", getOptionString(line, "minXCPSM"));
        params.add("minSpScorePSM", getOptionString(line, "minSPPSM"));
        params.add("minDeltaCnPSM", getOptionString(line, "minDCPSM"));
        params.add("minScorePSM",getOptionString(line, "minScorePSM"));
        params.add("minDeltaScorePSM", getOptionString(line, "minDSPSM"));
        params.add("minPEPPSM", getOptionString(line, "minPEPPSM"));
        params.add("nrSpScoreBins", getOptionString(line, "nrSPB"));
        params.add("smoothDegree", getOptionString(line, "smD"));
        params.add("fdrControlMethod", getOptionString(line, "fdrM"));
        params.add("groupingMethod", getOptionString(line, "groupM"));
        params.add("groupNames", getOptionString(line, "groupN"));
        params.add("groupRegExs", getOptionString(line, "groupRE"));
        params.add("reportAllPSM", getOptionString(line, "wAll"));


        params.finalize();

    }
}
