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

import newance.proteinmatch.OccamRazorSpectrumCounter;
import newance.proteinmatch.UniProtDB;
import newance.psmcombiner.*;
import newance.psmconverter.CometMultiplePepXMLFileConverter;
import newance.psmconverter.MaxQuantMultipleMSMSFileConverter;
import newance.psmconverter.PeptideSpectrumMatch;
import newance.util.*;
import org.apache.commons.cli.*;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

/**
 * @author Markus Müller
 */

public class CometMaxQuantScoreCombiner_PeptideTest extends ExecutableOptions {

    protected PsmGrouper psmGrouper;
    protected NewAnceParams params;
    protected SpectrumAccumulator spectrumAccumulator;

    protected Set<String> peptides;

    public CometMaxQuantScoreCombiner_PeptideTest() {

        createOptions();
        spectrumAccumulator = new SpectrumAccumulator();
    }

    public static void main(String[] args) {

        CometMaxQuantScoreCombiner_PeptideTest cometScoreCombiner =  new CometMaxQuantScoreCombiner_PeptideTest();
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
        CometMultiplePepXMLFileConverter cometMultiplePepXMLConverter = new CometMultiplePepXMLFileConverter(params.getCometPsmDir(), params.getCometPsmRegExp());
        cometMultiplePepXMLConverter.run();

        MaxQuantMultipleMSMSFileConverter maxQuantMultipleMSMSConverter = null;
        if (params.isIncludeMaxQuant()) {
            System.out.println("Parsing MaxQuant msms.txt and peptides.txt files ...");
            maxQuantMultipleMSMSConverter = new MaxQuantMultipleMSMSFileConverter(params.getMaxquantPsmDir(), params.getMaxquantPsmRegExp());
            maxQuantMultipleMSMSConverter.run();
        }

        System.out.println("Matching peptides to fasta file: "+params.getUniprotFastaFile()+" ...");
        System.out.println("Loading protein data from fasta file: "+params.getUniprotFastaFile()+" ...");

        UniProtDB uniProtDB = new UniProtDB(params.getUniprotFastaFile());
        long start = System.currentTimeMillis();
        cometMultiplePepXMLConverter.addDBProteins(uniProtDB);
        System.out.println("Comet DB matching ran in " + RunTime2String.getTimeDiffString(System.currentTimeMillis() - start));

        if (params.isIncludeMaxQuant()) {
            start = System.currentTimeMillis();
            maxQuantMultipleMSMSConverter.addDBProteins(uniProtDB);
            System.out.println("MaxQuant DB matching ran in " + RunTime2String.getTimeDiffString(System.currentTimeMillis() - start));
        }

        GroupedFDRCalculator groupedFDRCalculator = buildGroupedFDRCalculator(cometMultiplePepXMLConverter.getPsms());

        if (params.getFdrControlMethod().equals("combined")) {
            controlFDRCombined(groupedFDRCalculator, cometMultiplePepXMLConverter, maxQuantMultipleMSMSConverter);
        } else {
            controlFDRSeparate(groupedFDRCalculator, cometMultiplePepXMLConverter, maxQuantMultipleMSMSConverter);
        }

        writePeptideProteinGroupReport(uniProtDB);

        System.out.println("RunTime after Psm parsing: " + RunTime2String.getRuntimeString(Runtime.getRuntime()));


        return 0;
    }

    protected GroupedFDRCalculator buildGroupedFDRCalculator(ConcurrentHashMap<String, List<PeptideSpectrumMatch>> allPsms ) {

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
        if (params.isReportHistos()) groupedFDRCalculator.writeHistograms(params.getOutputDir()+File.separator+"histos", params.getOutputTag());
        groupedFDRCalculator.smoothHistogram(params.getSmoothDegree());
        groupedFDRCalculator.calcLocalFDR();
        if (params.isReportHistos()) groupedFDRCalculator.writeHistograms(params.getOutputDir()+File.separator+"histos", params.getOutputTag()+"_smoothed");

        return groupedFDRCalculator;
    }

    protected void controlFDRCombined(GroupedFDRCalculator groupedFDRCalculator, CometMultiplePepXMLFileConverter cometMultiplePepXMLConverter, MaxQuantMultipleMSMSFileConverter maxQuantMultipleMSMSConverter) {

        float lFDRThreshold = groupedFDRCalculator.calcLocalFDRThreshold((float)params.getFdrCometThreshold());


        for (String group : groupedFDRCalculator.getGroups()) {

            ConcurrentHashMap<String, List<PeptideSpectrumMatch>> filteredCometPsms = groupedFDRCalculator.filterPsms(cometMultiplePepXMLConverter.getPsms(), 2f, group);

            reportComet(filteredCometPsms,groupedFDRCalculator,lFDRThreshold,group);

            ConcurrentHashMap<String, List<PeptideSpectrumMatch>> combined;
            ConcurrentHashMap<String, List<PeptideSpectrumMatch>> maxQuantPsms = null;
            if (params.isIncludeMaxQuant()) {
                maxQuantPsms = ProcessPsmUtils.extractGroupPsms(psmGrouper, maxQuantMultipleMSMSConverter.getPsms(),group);
                reportMaxQuant(maxQuantPsms,group);

                combined = combine(filteredCometPsms,maxQuantPsms);

                reportCombined(combined, groupedFDRCalculator, lFDRThreshold, group);

            }
        }
    }

    protected void controlFDRSeparate(GroupedFDRCalculator groupedFDRCalculator, CometMultiplePepXMLFileConverter cometMultiplePepXMLConverter, MaxQuantMultipleMSMSFileConverter maxQuantMultipleMSMSConverter) {

        Map<String, Float> grpThresholdMap = groupedFDRCalculator.calcGroupLocalFDRThreshold((float)params.getFdrCometThreshold());

        for (String group : groupedFDRCalculator.getGroups()) {

            ConcurrentHashMap<String, List<PeptideSpectrumMatch>> filteredCometPsms = groupedFDRCalculator.filterPsms(cometMultiplePepXMLConverter.getPsms(), 2f, group);

            reportComet(filteredCometPsms,groupedFDRCalculator,grpThresholdMap.get(group),group);

            ConcurrentHashMap<String, List<PeptideSpectrumMatch>> combined;
            ConcurrentHashMap<String, List<PeptideSpectrumMatch>> maxQuantPsms = null;
            if (params.isIncludeMaxQuant()) {

                maxQuantPsms = ProcessPsmUtils.extractGroupPsms(psmGrouper, maxQuantMultipleMSMSConverter.getPsms(),group);
                reportMaxQuant(maxQuantPsms,group);

                combined = combine(filteredCometPsms,maxQuantPsms);
                reportCombined(combined, groupedFDRCalculator, grpThresholdMap.get(group), group);

            }
        }

    }

    protected void reportComet(ConcurrentHashMap<String, List<PeptideSpectrumMatch>>  psms, GroupedFDRCalculator groupedFDRCalculator, float lFDRThreshold, String group) {

        System.out.println("==============================");
        System.out.println("Comet "+group+" psms: ");
        System.out.println("==============================");
        System.out.println("");
        Psm2StringFunction psm2StringFunction = new Psm2StringFunction(Psm2StringFunction.TabStringMode.COMET, groupedFDRCalculator);
        System.out.println("Group\t"+psm2StringFunction.getHeader()+"\t"+"lFDR.Threshold");
        List<PeptideSpectrumMatch> tmp = new ArrayList<>();
        for (String specID : psms.keySet()) {
            for (PeptideSpectrumMatch psm : psms.get(specID)) {
                if (peptides.contains(psm.getPeptide().toSymbolString())) {
                    tmp.clear();
                    tmp.add(psm);
                    System.out.print(group+"\t"+psm2StringFunction.apply(specID,tmp));
                    System.out.println("\t"+String.format("%.4f ",lFDRThreshold));
                }
            }
        }
        System.out.println("");
    }


    protected void reportMaxQuant(ConcurrentHashMap<String, List<PeptideSpectrumMatch>>  psms, String group) {

        System.out.println("==============================");
        System.out.println("MaxQuant "+group+" psms: ");
        System.out.println("==============================");
        System.out.println("");
        Psm2StringFunction psm2StringFunction = new Psm2StringFunction(Psm2StringFunction.TabStringMode.MAXQUANT);
        System.out.println("Group\t"+psm2StringFunction.getHeader());
        List<PeptideSpectrumMatch> tmp = new ArrayList<>();
        for (String specID : psms.keySet()) {
            for (PeptideSpectrumMatch psm : psms.get(specID)) {
                if (peptides.contains(psm.getPeptide().toSymbolString())) {
                    tmp.clear();
                    tmp.add(psm);
                    System.out.println(group+"\t"+psm2StringFunction.apply(specID,tmp)) ;
                }
            }
        }
        System.out.println("");
    }

    protected void reportCombined(ConcurrentHashMap<String, List<PeptideSpectrumMatch>>  psms, GroupedFDRCalculator groupedFDRCalculator, float lFDRThreshold, String group) {

        System.out.println("==============================");
        System.out.println("Combined "+group+" psms: ");
        System.out.println("==============================");
        System.out.println("");
        Psm2StringFunction psm2StringFunction = new Psm2StringFunction(Psm2StringFunction.TabStringMode.COMBINED, groupedFDRCalculator);
        System.out.println("Group\t"+psm2StringFunction.getHeader()+"\t"+"Comet.lFDR.Threshold");
        List<PeptideSpectrumMatch> tmp = new ArrayList<>();
        for (String specID : psms.keySet()) {
            for (PeptideSpectrumMatch psm : psms.get(specID)) {
                if (peptides.contains(psm.getPeptide().toSymbolString())) {
                    tmp.clear();
                    tmp.add(psm);
                    System.out.print(group+"\t"+psm2StringFunction.apply(specID,tmp));
                    System.out.println("\t"+String.format("%.4f ",lFDRThreshold));
                }
            }
        }
        System.out.println("");
    }

    protected ConcurrentHashMap<String, List<PeptideSpectrumMatch>> combine(ConcurrentHashMap<String, List<PeptideSpectrumMatch>>  cometPsms, ConcurrentHashMap<String, List<PeptideSpectrumMatch>> maxQuantPsms) {

        ConcurrentHashMap<String, List<PeptideSpectrumMatch>> combined = new ConcurrentHashMap<>();

        CometMaxQuantPsmMerger combiner = new CometMaxQuantPsmMerger(maxQuantPsms,combined);
        cometPsms.forEach(combiner);

        return combined;
    }


    protected void writePeptideProteinGroupReport(UniProtDB uniProtDB) {

        if (!params.isDoPeptideProteinGrouping()) return;

        if (params.getSearchFastaFile() !=null) uniProtDB.addFastaFile(params.getSearchFastaFile());
        OccamRazorSpectrumCounter spectrumCounter = new OccamRazorSpectrumCounter(spectrumAccumulator, uniProtDB);

        String reportFileName = params.getOutputDir() + File.separator+params.getOutputTag()+"PeptideProteinGroupingReport.txt";
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
        cmdLineOpts.addOption(Option.builder("fdrM").required(false).hasArg().longOpt("fdrControlMethod").desc("Method to control pFDR: combined or separate (default combined).").build());
        cmdLineOpts.addOption(Option.builder("minXC").required(false).hasArg().longOpt("minXCorr").desc("Minimal Comet XCorr in histogram").build());
        cmdLineOpts.addOption(Option.builder("maxXC").required(false).hasArg().longOpt("maxXCorr").desc("Maximal Comet XCorr in histogram").build());
        cmdLineOpts.addOption(Option.builder("nrXCB").required(false).hasArg().longOpt("nrXCorrBins").desc("Number of Comet XCorr bins in histogram").build());
        cmdLineOpts.addOption(Option.builder("minSP").required(false).hasArg().longOpt("minSpScore").desc("Minimal Comet SpScore in histogram").build());
        cmdLineOpts.addOption(Option.builder("maxSP").required(false).hasArg().longOpt("maxSpScore").desc("Maximal Comet SpScore in histogram").build());
        cmdLineOpts.addOption(Option.builder("nrSPB").required(false).hasArg().longOpt("nrSpScoreBins").desc("Number of Comet SpScore bins in histogram").build());
        cmdLineOpts.addOption(Option.builder("minDC").required(false).hasArg().longOpt("minDeltaCn").desc("Minimal Comet DeltaCn in histogram").build());
        cmdLineOpts.addOption(Option.builder("maxDC").required(false).hasArg().longOpt("maxDeltaCn").desc("Maximal Comet DeltaCn in histogram").build());
        cmdLineOpts.addOption(Option.builder("nrDCB").required(false).hasArg().longOpt("nrDeltaCnBins").desc("Number of Comet DeltaCn bins in histogram").build());
        cmdLineOpts.addOption(Option.builder("rP").required(false).hasArg().longOpt("readParamFile").desc("Name of file from which parameters should to read.").build());
        cmdLineOpts.addOption(Option.builder("h").required(false).hasArg(false).longOpt("help").desc("Help option for command line help").build());
        cmdLineOpts.addOption(Option.builder("v").required(false).hasArg(false).longOpt("version").desc("Version of NewAnce software").build());
        cmdLineOpts.addOption(Option.builder("pepts").required(false).hasArg().longOpt("peptides").desc("Comma separated list of peptides to be printed.").build());
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

        this.peptides = getSetValue(getOptionString(line, "pepts"));

        params.finalize();
    }

    private Set<String> getSetValue(String value) {

        if (value.equals("[]")) {
            return new HashSet<>();
        } else if (value.startsWith("[") && value.endsWith("]") ) {
            String s = value.substring(1,value.length()-1);
            return new HashSet<>(Arrays.asList(s.split(",")));
        } else {
            throw new RuntimeException("Invalid value "+value+" for option pepts. Not a comma separated list value.");
        }
    }


}
