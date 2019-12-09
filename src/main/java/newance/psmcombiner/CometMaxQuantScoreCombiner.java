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
import java.util.regex.Pattern;

/**
 * Created by markusmueller on 12.03.18.
 */
public class CometMaxQuantScoreCombiner extends ExecutableOptions {

    protected String cometPsmDir;
    protected String maxquantPsmDir;
    protected Pattern cometPsmRegExp;
    protected Pattern maxquantPsmRegExp;
    protected boolean includeMaxQuant;
    protected String outputDir;
    protected boolean reportHistos;
    protected String searchFastaFile;
    protected String uniprotFastaFile;
    protected PsmGrouper psmGrouper;
    protected NewAnceParams params;
    protected SpectrumAccumulator spectrumAccumulator;
    protected boolean doPeptideProteinGrouping;


    public CometMaxQuantScoreCombiner() {

        createOptions();
        spectrumAccumulator = new SpectrumAccumulator();
    }

    public static void main(String[] args) {

        CometMaxQuantScoreCombiner cometScoreCombiner =  new CometMaxQuantScoreCombiner();
        try {
            cometScoreCombiner.parseOptions(args).init().run();
        } catch (MissingOptionException e) {
            cometScoreCombiner.checkHelpOption(args, "-h");
            cometScoreCombiner.checkVersionOption(args, NewAnceParams.getInstance().getVersion(), "-v");
            cometScoreCombiner.printOptions(args,e.getMessage());
        } catch (ParseException e) {
            cometScoreCombiner.printOptions(args, e.getMessage());
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public CometMaxQuantScoreCombiner init() throws IOException {

        return this;
    }

    public int run() throws IOException {

        System.out.println("Parsing Comet pep.xml files ...");
        CometPsmConverter cometPsmConverter = new CometPsmConverter(cometPsmDir, cometPsmRegExp);
        cometPsmConverter.run();

        MaxQuantPsmConverter maxQuantPsmConverter = null;
        if (includeMaxQuant) {
            System.out.println("Parsing MaxQuant msms.txt and peptides.txt files ...");
            maxQuantPsmConverter = new MaxQuantPsmConverter(maxquantPsmDir, maxquantPsmRegExp);
            maxQuantPsmConverter.run();
        }

        System.out.println("Matching peptides to fasta file: "+uniprotFastaFile+" ...");
        System.out.println("Loading protein data from fasta file: "+uniprotFastaFile+" ...");

        UniProtDB uniProtDB = new UniProtDB(uniprotFastaFile);
        long start = System.currentTimeMillis();
        cometPsmConverter.addDBProteins(uniProtDB);
        System.out.println("Comet DB matching ran in " + RunTime2String.getTimeDiffString(System.currentTimeMillis() - start));

        if (includeMaxQuant) {
            start = System.currentTimeMillis();
            maxQuantPsmConverter.addDBProteins(uniProtDB);
            System.out.println("MaxQuant DB matching ran in " + RunTime2String.getTimeDiffString(System.currentTimeMillis() - start));
        }

        GroupedFDRCalculator groupedFDRCalculator = buildGroupedFDRCalculator(cometPsmConverter.getPsms());

        if (params.getFdrControlMethod().equals("global")) {
            controlFDRGlobally(groupedFDRCalculator, cometPsmConverter, maxQuantPsmConverter);
        } else {

        }

        writePeptideProteinGroupReport(uniProtDB);

//        groupedFDRCalculator.calcPValuesFDR();

        System.out.println("RunTime after Psm parsing: " + RunTime2String.getRuntimeString(Runtime.getRuntime()));


        return 0;
    }

    protected GroupedFDRCalculator buildGroupedFDRCalculator(ConcurrentHashMap<String, List<PeptideMatchData>> allPsms ) {

        if (params.getStdProtRegExp()!=null) {
            if (params.getExcludedProtPattern().isEmpty())
                psmGrouper = new RegExpProteinGrouper(params.getStdProtRegExp(), params.getStandartGroup(), params.getCrypticGroup());
            else
                psmGrouper = new RegExpProteinGrouper(params.getStdProtRegExp(), params.getForcedCrypticProts(), params.getStandartGroup(), params.getCrypticGroup());
        } else {
            psmGrouper = new ModificationPSMGrouper();
        }

        GroupedFDRCalculator groupedFDRCalculator = new GroupedFDRCalculator(psmGrouper);
        groupedFDRCalculator.addAll(allPsms);
        groupedFDRCalculator.setCanCalculateFDR(params.getMinNrPsmsPerHisto());
        groupedFDRCalculator.calcClassProbs();
        if (reportHistos) groupedFDRCalculator.writeHistograms(outputDir+File.separator+"histos", params.getOutputPrefix());
        groupedFDRCalculator.smoothHistogram(3);
        groupedFDRCalculator.calcLocalFDR();
        if (reportHistos) groupedFDRCalculator.writeHistograms(outputDir+File.separator+"histos", params.getOutputPrefix()+"_smoothed");

//        float dfdr = 1f/100;
//        for (int i=0;i<100;i++) {
//            float fdr = i*dfdr;
//            System.out.println("lFDR = "+fdr+", pFDR = "+groupedFDRCalculator.calcGlobalFDR(fdr));
//        }

        return groupedFDRCalculator;
    }

    protected void controlFDRGlobally(GroupedFDRCalculator groupedFDRCalculator, CometPsmConverter cometPsmConverter, MaxQuantPsmConverter maxQuantPsmConverter) {

        float lFDRThreshold = groupedFDRCalculator.calcLocalFDRThreshold((float)params.getFdrCometThreshold());

        System.out.print(groupedFDRCalculator.printTree(lFDRThreshold));
        test(cometPsmConverter.getPsms(),groupedFDRCalculator,lFDRThreshold);

        SummaryReportWriter summaryReportWriter = new SummaryReportWriter(outputDir +File.separator+ "SummaryReport.txt", includeMaxQuant);

        for (String group : groupedFDRCalculator.getGroups()) {

            long start = System.currentTimeMillis();
            ConcurrentHashMap<String, List<PeptideMatchData>> filteredCometPsms = groupedFDRCalculator.filterPsms(cometPsmConverter.getPsms(), lFDRThreshold, group);
            System.out.println("Comet FDR filtering ran in " + RunTime2String.getTimeDiffString(System.currentTimeMillis() - start));


            test(filteredCometPsms,groupedFDRCalculator,lFDRThreshold,group);

            ConcurrentHashMap<String, List<PeptideMatchData>> combined;
            ConcurrentHashMap<String, List<PeptideMatchData>> maxQuantPsms = null;
            if (includeMaxQuant) {

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

            if (includeMaxQuant) {
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

        SummaryReportWriter summaryReportWriter = new SummaryReportWriter(outputDir +File.separator+ "SummaryReport.txt", includeMaxQuant);

        for (String group : groupedFDRCalculator.getGroups()) {

            long start = System.currentTimeMillis();
            ConcurrentHashMap<String, List<PeptideMatchData>> filteredCometPsms = groupedFDRCalculator.filterPsms(cometPsmConverter.getPsms(), grpThresholdMap.get(group), group);
            System.out.println("Comet FDR filtering ran in " + RunTime2String.getTimeDiffString(System.currentTimeMillis() - start));

            ConcurrentHashMap<String, List<PeptideMatchData>> combined;
            ConcurrentHashMap<String, List<PeptideMatchData>> maxQuantPsms = null;
            if (includeMaxQuant) {

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

            if (includeMaxQuant) {
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

        StringFileWriter writer = new StringFileWriter(outputDir +"/"+filename, stringFunction);

        psms.forEach(10000, stringFunction, writer);

        writer.close();
    }

    protected void writeToCometTabFile(ConcurrentHashMap<String, List<PeptideMatchData>>  psms, String filename) throws IOException {

        final Psm2StringFunction stringFunction = new Psm2StringFunction("tab", Psm2StringFunction.TabStringMode.COMET);

        StringFileWriter writer = new StringFileWriter(outputDir +"/"+filename, stringFunction);

        psms.forEach(10000, stringFunction, writer);

        writer.close();
    }

    protected void writeToMaxQuantTabFile(ConcurrentHashMap<String, List<PeptideMatchData>>  psms, String filename) throws IOException {

        final Psm2StringFunction stringFunction = new Psm2StringFunction("tab", Psm2StringFunction.TabStringMode.MAXQUANT);

        StringFileWriter writer = new StringFileWriter(outputDir +"/"+filename, stringFunction);

        psms.forEach(10000, stringFunction, writer);

        writer.close();
    }


    protected void writeToFastaFile(ConcurrentHashMap<String, List<PeptideMatchData>>  psms, String filename) throws IOException {

        final Psm2StringFunction stringFunction = new Psm2StringFunction("fasta");

        StringFileWriter writer = new StringFileWriter(outputDir +"/"+filename, stringFunction, true);

        psms.forEach(10000, stringFunction, writer);

        writer.close();
    }


    protected void writeToTextFile(ConcurrentHashMap<String, List<PeptideMatchData>>  psms, String filename) throws IOException {

        final Psm2StringFunction stringFunction = new Psm2StringFunction("txt");

        StringFileWriter writer = new StringFileWriter(outputDir +"/"+filename, stringFunction, true);

        psms.forEach(10000, stringFunction, writer);

        writer.close();
    }

    protected void writePeptideProteinGroupReport(UniProtDB uniProtDB) {

        if (searchFastaFile !=null) uniProtDB.addFastaFile(searchFastaFile);
        OccamRazorSpectrumCounter spectrumCounter = new OccamRazorSpectrumCounter(spectrumAccumulator, uniProtDB);

        String reportFileName = outputDir + File.separator +"PeptideProteinGroupingReport.txt";
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
        cmdLineOpts.addOption(Option.builder("mqRE").required(false).hasArg().longOpt("maxquantPsmRegex").desc("Regular expression of MaxQuant psm files (e.g. msms\\.txt$). If not provided only Comet is used.").build());
        cmdLineOpts.addOption(Option.builder("coFDR").required().hasArg().longOpt("cometFDR").desc("FDR for filtering Comet PSMs before combination (required)").build());
        cmdLineOpts.addOption(Option.builder("outD").required().hasArg().longOpt("outputDir").desc("Output directory for results (required)").build());
        cmdLineOpts.addOption(Option.builder("repH").required(false).hasArg(false).longOpt("reportHistogram").desc("Report histograms to text files").build());
        cmdLineOpts.addOption(Option.builder("protG").required(false).hasArg().longOpt("proteinGroup").desc("Name of group with protein coding or canonical sequences").build());
        cmdLineOpts.addOption(Option.builder("crypG").required(false).hasArg().longOpt("crypticGroup").desc("Name of group with non-canonical or cryptic sequences").build());
        cmdLineOpts.addOption(Option.builder("stdRE").required(false).hasArg().longOpt("stdProtRegExp").desc("Regular expression to match fasta name of coding proteins (e.g. sp\\||tr\\| ").build());
        cmdLineOpts.addOption(Option.builder("spRE").required(false).hasArg().longOpt("spectrumFilter").desc("If this option is set, only spectrum files that match the regexp are used").build());
        cmdLineOpts.addOption(Option.builder("exclP").required(false).hasArg().longOpt("excludeCodingProts").desc("Comma separated list of protein names to be excluded from coding prots (e.g. PGBD5_HUMAN,POGZ_HUMAN,PGBD1_HUMAN)").build());
        cmdLineOpts.addOption(Option.builder("seFa").required(false).hasArg().longOpt("searchFastaFile").desc("Fasta file used for search (for protein grouping export)").build());
        cmdLineOpts.addOption(Option.builder("upFa").required(false).hasArg().longOpt("uniProtFastaFile").desc("Fasta file with coding or canonical proteins (e.g. UniProt fasta file)").build());
        cmdLineOpts.addOption(Option.builder("ppG").required(false).hasArg(false).longOpt("peptideProteinGrouping").desc("Perform peptide protein grouping export").build());
        cmdLineOpts.addOption(Option.builder("maxR").required(false).hasArg().longOpt("maxRank").desc("Maximal rank of peptide in list of spectrum matches (rank 1 = best) (default value: 1)").build());
        cmdLineOpts.addOption(Option.builder("minZ").required(false).hasArg().longOpt("minCharge").desc("Minimal charge of PSM (default value: 1)").build());
        cmdLineOpts.addOption(Option.builder("maxZ").required(false).hasArg().longOpt("maxCharge").desc("Maximal charge of PSM (default value: 5)").build());
        cmdLineOpts.addOption(Option.builder("minL").required(false).hasArg().longOpt("minLength").desc("Minimal length of peptide (default value: 8)").build());
        cmdLineOpts.addOption(Option.builder("maxL").required(false).hasArg().longOpt("maxLength").desc("Maximal length of peptide (default value: 25)").build());
        cmdLineOpts.addOption(Option.builder("nrTh").required(false).hasArg().longOpt("nrThreads").desc("Number of threads used by NewAnce (default value: nr of available processors - 2)").build());
        cmdLineOpts.addOption(Option.builder("outP").required(false).hasArg().longOpt("outputPrefix").desc("Prefix for output files. If this option is not set no prefix is used.").build());
        cmdLineOpts.addOption(Option.builder("minPH").required(false).hasArg().longOpt("minPsm4Histo").desc("Minimal number of psms to calculate local FDR in histogram (default 100000).").build());
        cmdLineOpts.addOption(Option.builder("fdrM").required(false).hasArg().longOpt("fdrControlMethod").desc("Method to control pFDR: global or groupwise (default global).").build());
        cmdLineOpts.addOption(Option.builder("h").required(false).hasArg(false).longOpt("help").desc("Help option for command line help").build());
        cmdLineOpts.addOption(Option.builder("v").required(false).hasArg(false).longOpt("version").desc("Version of NewAnce software").build());
    }

    @Override
    protected void check(CommandLine line) throws ParseException {

        this.cometPsmDir = checkExistFileOption(line,"coD");
        this.cometPsmRegExp = checkRegExOption(line, "coRE");

        this.maxquantPsmDir = checkExistFileOption(line,"mqD");
        this.maxquantPsmRegExp = checkRegExOption(line, "mqRE");

        if (maxquantPsmDir.isEmpty() || maxquantPsmRegExp==null) this.includeMaxQuant = false;
        else this.includeMaxQuant = true;

        this.reportHistos = checkBooleanOption(line, "repH");
        this.outputDir = checkOutputDirOption(line, "outD");
        this.searchFastaFile = checkExistFileOption(line, "seFa");
        this.uniprotFastaFile = checkExistFileOption(line, "upFa");
        this.doPeptideProteinGrouping = checkBooleanOption(line, "ppG");

        this.params = NewAnceParams.getInstance();

        this.params.setMaxRank(checkIntOption(line,"maxR",1,100,1));
        this.params.setMinCharge(checkIntOption(line,"minZ",1,100,1));
        this.params.setMaxCharge(checkIntOption(line,"maxZ",1,100,5));
        this.params.setMinPeptideLength(checkIntOption(line,"minL",1,1000,8));
        this.params.setMaxPeptideLength(checkIntOption(line,"maxL",1,1000,25));
        this.params.setFdrCometThreshold(checkDoubleOption(line,"coFDR",0.0,100000000.0,0.03));
        this.params.setStandartGroup(checkStringOption(line,"protG"));
        this.params.setCrypticGroup(checkStringOption(line,"crypG"));
        this.params.setForcedCrypticProts(checkStringListOption(line,"exclP"));
        this.params.setSpectrumRegExp(checkRegExOption(line,"spRE"));
        this.params.setStdProtRegExp(checkRegExOption(line,"stdRE"));
        this.params.setOutputPrefix(checkStringOption(line,"outP"));
        this.params.setMinNrPsmsPerHisto(checkIntOption(line,"minPH",1,1000000000,50000));

        int nrProcs = Runtime.getRuntime().availableProcessors();
        this.params.setNrThreads(checkIntOption(line,"nrTh",1,nrProcs,nrProcs-2));

        Set<String> allowedStrValues = new HashSet<>();
        allowedStrValues.add("global");
        allowedStrValues.add("groupwise");
        this.params.setFdrControlMethod(checkDefinedStringOption(line,"fdrM", allowedStrValues));
    }

}
