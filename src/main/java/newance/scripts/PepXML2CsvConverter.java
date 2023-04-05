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

import newance.psmcombiner.*;
import newance.psmconverter.CometMultiplePepXMLFileConverter;
import newance.psmconverter.PeptideSpectrumMatch;
import newance.util.ExecutableOptions;
import newance.util.NewAnceParams;
import newance.util.StringFileWriter;
import org.apache.commons.cli.*;

import java.io.File;
import java.io.IOException;
import java.nio.file.InvalidPathException;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * @author Markus Müller
 */

public class PepXML2CsvConverter extends ExecutableOptions {

    protected NewAnceParams params;
    protected SpectrumAccumulator spectrumAccumulator;
    protected File csvFile;

    public PepXML2CsvConverter() {

        version = NewAnceParams.getInstance().getVersion();
        createOptions();
        spectrumAccumulator = new SpectrumAccumulator();
    }

    public static void main(String[] args) {

        PepXML2CsvConverter cometScoreCombiner =  new PepXML2CsvConverter();
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

        CometMultiplePepXMLFileConverter cometMultiplePepXMLConverter =
                new CometMultiplePepXMLFileConverter(params.getCometPsmDir(), params.getCometPsmRegExp(), groupedFDRCalculator, true);
        cometMultiplePepXMLConverter.run();

        System.out.println("Write results to csv ...");
        writePSMTabFile(cometMultiplePepXMLConverter.getPsms(), this.csvFile);

        return 0;
    }

    protected void createOptions() {

        this.cmdLineOpts = new Options();

        cmdLineOpts.addOption(Option.builder("coD").required().hasArg().longOpt("cometPsmDir").desc("Comet psm root directory (required)").build());
        cmdLineOpts.addOption(Option.builder("coRE").required().hasArg().longOpt("cometPsmRegex").desc("Regular expression of Comet psm files (e.g. \\.xml$) (required)").build());
        cmdLineOpts.addOption(Option.builder("csvF").required().hasArg().longOpt("csvFile").desc("Output csv file (required)").build());
        cmdLineOpts.addOption(Option.builder("spRE").required(false).hasArg().longOpt("spectrumFilter").desc("If this option is set, only spectrum ids that match this regexp are used.  If not set no filtering is performed.").build());
        cmdLineOpts.addOption(Option.builder("mod").required(false).hasArg().longOpt("modifications").desc("Comma separated list of peptide modifications used in search (e.g. Cysteinyl:C3H5NO2S,Oxidation:O)").build());
        cmdLineOpts.addOption(Option.builder("maxR").required(false).hasArg().longOpt("maxRank").desc("Maximal rank of peptide in list of spectrum matches (rank 1 = best) (default value: 1)").build());
        cmdLineOpts.addOption(Option.builder("minZ").required(false).hasArg().longOpt("minCharge").desc("Minimal charge of PSM (default value: 1)").build());
        cmdLineOpts.addOption(Option.builder("maxZ").required(false).hasArg().longOpt("maxCharge").desc("Maximal charge of PSM (default value: 5)").build());
        cmdLineOpts.addOption(Option.builder("minL").required(false).hasArg().longOpt("minLength").desc("Minimal length of peptide (default value: 8)").build());
        cmdLineOpts.addOption(Option.builder("maxL").required(false).hasArg().longOpt("maxLength").desc("Maximal length of peptide (default value: 25)").build());
        cmdLineOpts.addOption(Option.builder("nrTh").required(false).hasArg().longOpt("nrThreads").desc("Number of threads used by NewAnce (default value: nr of available processors - 2)").build());
        cmdLineOpts.addOption(Option.builder("h").required(false).hasArg(false).longOpt("help").desc("Help option for command line help").build());
        cmdLineOpts.addOption(Option.builder("v").required(false).hasArg(false).longOpt("version").desc("Version of NewAnce software").build());
    }

    @Override
    protected void check(CommandLine line) throws ParseException {

        if (optionsSet) return;

        this.params = NewAnceParams.getInstance();

        params.add("cometPsmDir", getOptionString(line, "coD"));
        params.add("cometPsmRegExp", getOptionString(line, "coRE"));
        params.add("maxRank", getOptionString(line, "maxR"));
        params.add("minCharge", getOptionString(line, "minZ"));
        params.add("maxCharge", getOptionString(line, "maxZ"));
        params.add("minPeptideLength", getOptionString(line, "minL"));
        params.add("maxPeptideLength", getOptionString(line, "maxL"));
        params.add("spectrumRegExp", getOptionString(line, "spRE"));
        params.add("modifications", getOptionString(line, "mod"));
        params.add("nrThreads", getOptionString(line, "nrTh"));

        this.csvFile = new File(getOptionString(line, "csvF"));

        params.finalize();
    }

    protected void writePSMTabFile(ConcurrentHashMap<String, List<PeptideSpectrumMatch>> psms, File filename)  {

        final Psm2StringFunction stringFunction;

        stringFunction = new PepXMLPSm2StringFunction();

        StringFileWriter writer =
                new StringFileWriter(filename, stringFunction);

        psms.forEach(10000, stringFunction, writer);

        writer.close();
    }


}
