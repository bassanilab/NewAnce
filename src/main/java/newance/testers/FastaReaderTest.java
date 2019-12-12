package newance.testers;

import newance.proteinmatch.UniProtDB;
import newance.proteinmatch.UniProtProtein;
import newance.psmcombiner.ModificationPSMGrouper;
import newance.psmcombiner.Psm2StringFunction;
import newance.psmcombiner.RegExpProteinGrouper;
import newance.psmconverter.CometPsmConverter;
import newance.psmconverter.PeptideMatchData;
import newance.util.ExecutableOptions;
import newance.util.NewAnceParams;
import newance.util.PsmGrouper;
import newance.util.RunTime2String;
import org.apache.commons.cli.*;

import java.io.IOException;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Copyright (C) 2019
 * @author Markus Müller
 * @Institutions: SIB, Swiss Institute of Bioinformatics; Ludwig Institute for Cancer Research
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

        System.out.println("ID\tProtName\tGeneName\tDescription\tDB\tSequence\n");

        for (UniProtProtein prot : uniProtDB.getProteins()) {
            if(cnt>=maxNrDisplayedPSMs) break;

            System.out.print(prot.getUniProtAC()+"\t"+prot.getUniProtName()+"\t"+prot.getGeneName()+"\t"+prot.getDescription()+"\t");
            System.out.println(prot.getDbFlag()+"\t"+prot.getSequence());
            cnt++;
        }

        return 0;
    }


    protected void createOptions() {

        this.cmdLineOpts = new Options();

        cmdLineOpts.addOption(Option.builder("fa").required(false).hasArg().longOpt("uniProtFastaFile").desc("Fasta file with coding or canonical proteins (e.g. UniProt fasta file)").build());
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
