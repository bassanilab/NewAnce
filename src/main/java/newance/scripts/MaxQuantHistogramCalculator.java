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

import newance.psmcombiner.GroupedFDRCalculator;
import newance.psmcombiner.SpectrumAccumulator;
import newance.psmconverter.MaxQuantMultipleMSMSFileConverter;
import newance.util.ExecutableOptions;
import newance.util.NewAnceParams;
import newance.util.PsmGrouper;
import org.apache.commons.cli.*;

import java.io.IOException;

/**
 * @author Markus Müller
 */

public class MaxQuantHistogramCalculator extends ExecutableOptions {

    protected NewAnceParams params;
    protected SpectrumAccumulator spectrumAccumulator;

    public MaxQuantHistogramCalculator() {

        version = NewAnceParams.getInstance().getVersion();
        createOptions();
        spectrumAccumulator = new SpectrumAccumulator();
    }

    public static void main(String[] args) {

        MaxQuantHistogramCalculator maxQuantScoreCombiner =  new MaxQuantHistogramCalculator();
        try {
            maxQuantScoreCombiner.init(args).parseOptions(args).run();
        } catch (MissingOptionException e) {
        } catch (ParseException e) {
            maxQuantScoreCombiner.printOptions(args, e.getMessage());
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public int run() throws IOException {

        System.out.println("Build groupedFDRCalculator ...");

        GroupedFDRCalculator groupedFDRCalculator = new GroupedFDRCalculator(null, NewAnceParams.SearchTool.MAXQUANT);

        MaxQuantMultipleMSMSFileConverter maxQuantMultiplePepXMLConverter =
                new MaxQuantMultipleMSMSFileConverter(params.getMaxquantPsmDir(), params.getMaxquantPsmRegExp(), groupedFDRCalculator, true);
        maxQuantMultiplePepXMLConverter.run();

        System.out.println("Write histograms ...");
        groupedFDRCalculator.writeHistograms(params.getOutputDir(), params.getOutputTag(),1);

        return 0;
    }

    protected void createOptions() {

        this.cmdLineOpts = new Options();

        cmdLineOpts.addOption(Option.builder("mqD").required().hasArg().longOpt("maxquantPsmDir").desc("MaxQuant psm root directory (required)").build());
        cmdLineOpts.addOption(Option.builder("outD").required().hasArg().longOpt("outputDir").desc("Output directory for results (required)").build());
        cmdLineOpts.addOption(Option.builder("rMqH").required(false).hasArg().longOpt("readMaxQuantHistograms").desc("Directory where histograms files are placed.").build());
        cmdLineOpts.addOption(Option.builder("spRE").required(false).hasArg().longOpt("spectrumFilter").desc("If this option is set, only spectrum ids that match this regexp are used.  If not set no filtering is performed.").build());
        cmdLineOpts.addOption(Option.builder("mod").required(false).hasArg().longOpt("modifications").desc("Comma separated list of peptide modifications used in search (e.g. Cysteinyl:C3H5NO2S,Oxidation:O)").build());
        cmdLineOpts.addOption(Option.builder("maxR").required(false).hasArg().longOpt("maxRank").desc("Maximal rank of peptide in list of spectrum matches (rank 1 = best) (default value: 1)").build());
        cmdLineOpts.addOption(Option.builder("minZ").required(false).hasArg().longOpt("minCharge").desc("Minimal charge of PSM (default value: 1)").build());
        cmdLineOpts.addOption(Option.builder("maxZ").required(false).hasArg().longOpt("maxCharge").desc("Maximal charge of PSM (default value: 5)").build());
        cmdLineOpts.addOption(Option.builder("minL").required(false).hasArg().longOpt("minLength").desc("Minimal length of peptide (default value: 8)").build());
        cmdLineOpts.addOption(Option.builder("maxL").required(false).hasArg().longOpt("maxLength").desc("Maximal length of peptide (default value: 25)").build());
        cmdLineOpts.addOption(Option.builder("nrTh").required(false).hasArg().longOpt("nrThreads").desc("Number of threads used by NewAnce (default value: nr of available processors - 2)").build());
        cmdLineOpts.addOption(Option.builder("outT").required(false).hasArg().longOpt("outputTag").desc("Tag inserted into output file names after prefix.").build());
        cmdLineOpts.addOption(Option.builder("minScore").required(false).hasArg().longOpt("minScore").desc("Minimal MaxQuant Score in histogram (default value 0)").build());
        cmdLineOpts.addOption(Option.builder("maxScore").required(false).hasArg().longOpt("maxScore").desc("Maximal MaxQuant Score in histogram (default value 5)").build());
        cmdLineOpts.addOption(Option.builder("nrScoreB").required(false).hasArg().longOpt("nrScoreBins").desc("Number of MaxQuant Score bins in histogram (default value 40)" ).build());
        cmdLineOpts.addOption(Option.builder("minDS").required(false).hasArg().longOpt("minDeltaScore").desc("Minimal MaxQuant DeltaScore in histogram (default value 0)").build());
        cmdLineOpts.addOption(Option.builder("maxDS").required(false).hasArg().longOpt("maxDeltaScore").desc("Maximal MaxQuant DeltaScore in histogram (default value 1)").build());
        cmdLineOpts.addOption(Option.builder("nrDSB").required(false).hasArg().longOpt("nrDeltaScoreBins").desc("MaxQuant of MaxQuant DeltaScore bins in histogram (default value 40)").build());
        cmdLineOpts.addOption(Option.builder("minPEP").required(false).hasArg().longOpt("minPEP").desc("Minimal MaxQuant PEP in histogram (default value 0)").build());
        cmdLineOpts.addOption(Option.builder("maxPEP").required(false).hasArg().longOpt("maxPEP").desc("Maximal MaxQuant PEP in histogram (default value 2500)").build());
        cmdLineOpts.addOption(Option.builder("nrPEPB").required(false).hasArg().longOpt("nrPEPBins").desc("Number of MaxQuant PEP bins in histogram (default value 40)").build());
        cmdLineOpts.addOption(Option.builder("wP").required(false).hasArg().longOpt("write2ParamFile").desc("Filename where parameters should be written to.").build());
        cmdLineOpts.addOption(Option.builder("rP").required(false).hasArg().longOpt("readParamFile").desc("Name of file from which parameters should to read.").build());
        cmdLineOpts.addOption(Option.builder("h").required(false).hasArg(false).longOpt("help").desc("Help option for command line help").build());
        cmdLineOpts.addOption(Option.builder("v").required(false).hasArg(false).longOpt("version").desc("Version of NewAnce software").build());
    }

    @Override
    protected void check(CommandLine line) throws ParseException {

        if (optionsSet) return;

        this.params = NewAnceParams.getInstance();

        params.add("maxquantPsmDir", getOptionString(line, "mqD"));
        params.add("readMaxQuantHistos", getOptionString(line, "rMqH"));
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

        params.finalize();
    }
}
