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
import newance.psmcombiner.*;
import newance.psmconverter.CometMultiplePepXMLFileConverter;
import newance.psmconverter.MaxQuantMultipleMSMSFileConverter;
import newance.psmconverter.PeptideSpectrumMatch;
import newance.util.*;
import org.apache.commons.cli.*;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * @author Markus Müller
 */

public class CometHistogramCalculator extends ExecutableOptions {

    protected PsmGrouper psmGrouper;
    protected NewAnceParams params;
    protected SpectrumAccumulator spectrumAccumulator;

    public CometHistogramCalculator() {

        version = NewAnceParams.getInstance().getVersion();
        createOptions();
        spectrumAccumulator = new SpectrumAccumulator();
    }

    public static void main(String[] args) {

        CometHistogramCalculator cometScoreCombiner =  new CometHistogramCalculator();
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

        System.out.println("Build groupedFDRCalculator ...");

        GroupedFDRCalculator groupedFDRCalculator = new GroupedFDRCalculator();

        CometMultiplePepXMLFileConverter cometMultiplePepXMLConverter = new CometMultiplePepXMLFileConverter(params.getCometPsmDir(), params.getCometPsmRegExp(), groupedFDRCalculator, true);
        cometMultiplePepXMLConverter.run();

        System.out.println("Write histograms ...");
        groupedFDRCalculator.writeHistograms(params.getOutputDir(), params.getOutputTag(),1);

        return 0;
    }

    protected void createOptions() {

        this.cmdLineOpts = new Options();

        cmdLineOpts.addOption(Option.builder("coD").required().hasArg().longOpt("cometPsmDir").desc("Comet psm root directory (required)").build());
        cmdLineOpts.addOption(Option.builder("coRE").required().hasArg().longOpt("cometPsmRegex").desc("Regular expression of Comet psm files (e.g. \\.xml$) (required)").build());
        cmdLineOpts.addOption(Option.builder("outD").required().hasArg().longOpt("outputDir").desc("Output directory for results (required)").build());
        cmdLineOpts.addOption(Option.builder("readH").required(false).hasArg().longOpt("readHistograms").desc("Directory where histograms files are placed.").build());
        cmdLineOpts.addOption(Option.builder("spRE").required(false).hasArg().longOpt("spectrumFilter").desc("If this option is set, only spectrum ids that match this regexp are used.  If not set no filtering is performed.").build());
        cmdLineOpts.addOption(Option.builder("mod").required(false).hasArg().longOpt("modifications").desc("Comma separated list of peptide modifications used in search (e.g. Cysteinyl:C3H5NO2S,Oxidation:O)").build());
        cmdLineOpts.addOption(Option.builder("maxR").required(false).hasArg().longOpt("maxRank").desc("Maximal rank of peptide in list of spectrum matches (rank 1 = best) (default value: 1)").build());
        cmdLineOpts.addOption(Option.builder("minZ").required(false).hasArg().longOpt("minCharge").desc("Minimal charge of PSM (default value: 1)").build());
        cmdLineOpts.addOption(Option.builder("maxZ").required(false).hasArg().longOpt("maxCharge").desc("Maximal charge of PSM (default value: 5)").build());
        cmdLineOpts.addOption(Option.builder("minL").required(false).hasArg().longOpt("minLength").desc("Minimal length of peptide (default value: 8)").build());
        cmdLineOpts.addOption(Option.builder("maxL").required(false).hasArg().longOpt("maxLength").desc("Maximal length of peptide (default value: 25)").build());
        cmdLineOpts.addOption(Option.builder("nrTh").required(false).hasArg().longOpt("nrThreads").desc("Number of threads used by NewAnce (default value: nr of available processors - 2)").build());
        cmdLineOpts.addOption(Option.builder("outT").required(false).hasArg().longOpt("outputTag").desc("Tag inserted into output file names after prefix.").build());
        cmdLineOpts.addOption(Option.builder("minXC").required(false).hasArg().longOpt("minXCorr").desc("Minimal Comet XCorr in histogram (default value 0)").build());
        cmdLineOpts.addOption(Option.builder("maxXC").required(false).hasArg().longOpt("maxXCorr").desc("Maximal Comet XCorr in histogram (default value 5)").build());
        cmdLineOpts.addOption(Option.builder("nrXCB").required(false).hasArg().longOpt("nrXCorrBins").desc("Number of Comet XCorr bins in histogram (default value 40)" ).build());
        cmdLineOpts.addOption(Option.builder("minSP").required(false).hasArg().longOpt("minSpScore").desc("Minimal Comet SpScore in histogram (default value 0)").build());
        cmdLineOpts.addOption(Option.builder("maxSP").required(false).hasArg().longOpt("maxSpScore").desc("Maximal Comet SpScore in histogram (default value 1)").build());
        cmdLineOpts.addOption(Option.builder("nrSPB").required(false).hasArg().longOpt("nrSpScoreBins").desc("Number of Comet SpScore bins in histogram (default value 40)").build());
        cmdLineOpts.addOption(Option.builder("minDC").required(false).hasArg().longOpt("minDeltaCn").desc("Minimal Comet DeltaCn in histogram (default value 0)").build());
        cmdLineOpts.addOption(Option.builder("maxDC").required(false).hasArg().longOpt("maxDeltaCn").desc("Maximal Comet DeltaCn in histogram (default value 2500)").build());
        cmdLineOpts.addOption(Option.builder("nrDCB").required(false).hasArg().longOpt("nrDeltaCnBins").desc("Number of Comet DeltaCn bins in histogram (default value 40)").build());
        cmdLineOpts.addOption(Option.builder("wP").required(false).hasArg().longOpt("write2ParamFile").desc("Filename where parameters should be written to.").build());
        cmdLineOpts.addOption(Option.builder("rP").required(false).hasArg().longOpt("readParamFile").desc("Name of file from which parameters should to read.").build());
        cmdLineOpts.addOption(Option.builder("h").required(false).hasArg(false).longOpt("help").desc("Help option for command line help").build());
        cmdLineOpts.addOption(Option.builder("v").required(false).hasArg(false).longOpt("version").desc("Version of NewAnce software").build());
    }

    @Override
    protected void check(CommandLine line) throws ParseException {

        if (optionsSet) return;

        this.params = NewAnceParams.getInstance();

        params.add("cometPsmDir", getOptionString(line, "coD"));
        params.add("cometPsmRegExp", getOptionString(line, "coRE"));
        params.add("readHistos", getOptionString(line, "readH"));
        params.add("outputDir", getOptionString(line, "outD"));
        params.add("writeParamsFile", getOptionString(line, "wP"));
        params.add("readParamsFile", getOptionString(line, "rP"));
        params.add("maxRank", getOptionString(line, "maxR"));
        params.add("minCharge", getOptionString(line, "minZ"));
        params.add("maxCharge", getOptionString(line, "maxZ"));
        params.add("minPeptideLength", getOptionString(line, "minL"));
        params.add("maxPeptideLength", getOptionString(line, "maxL"));
        params.add("spectrumRegExp", getOptionString(line, "spRE"));
        params.add("outputTag", getOptionString(line, "outT"));
        params.add("modifications", getOptionString(line, "mod"));
        params.add("minXCorr", getOptionString(line, "minXC"));
        params.add("maxXCorr", getOptionString(line, "maxXC"));
        params.add("nrXCorrBins", getOptionString(line, "nrXCB"));
        params.add("minDeltaCn", getOptionString(line, "minDC"));
        params.add("maxDeltaCn", getOptionString(line, "maxDC"));
        params.add("nrDeltaCnBins", getOptionString(line, "nrDCB"));
        params.add("minSpScore", getOptionString(line, "minSP"));
        params.add("maxSpScore", getOptionString(line, "maxSP"));
        params.add("nrSpScoreBins", getOptionString(line, "nrSPB"));
        params.add("nrThreads", getOptionString(line, "nrTh"));

        params.finalize();

    }
}
