/*
Copyright (C) SIB - Swiss Institute of Bioinformatics, Lausanne, Switzerland
Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
*/

package newance.testers;

import newance.proteinmatch.UniProtDB;
import newance.psmcombiner.*;
import newance.psmconverter.CometMultiplePepXMLFileConverter;
import newance.psmconverter.PeptideSpectrumMatch;
import newance.util.*;
import org.apache.commons.cli.*;

import java.io.IOException;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.regex.Pattern;

/**
 * @author Markus Müller
 */

public class CometXMLReaderTest extends ExecutableOptions {

    protected NewAnceParams params;
    protected int maxNrDisplayedPSMs;
    protected final Psm2StringFunction stringFunction;


    public CometXMLReaderTest() {

        stringFunction = new Psm2StringFunction(Psm2StringFunction.TabStringMode.COMET);
        maxNrDisplayedPSMs = -1;
        createOptions();
    }

    public static void main(String[] args) {

        CometXMLReaderTest cometXMLReaderTest =  new CometXMLReaderTest();
        try {
            cometXMLReaderTest.init(args).parseOptions(args).run();
        } catch (MissingOptionException e) {
        } catch (ParseException e) {
            cometXMLReaderTest.printOptions(args, e.getMessage());
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public CometXMLReaderTest init() throws IOException {

        return this;
    }

    public int run() throws IOException {

        System.out.println("Parsing Comet pep.xml files ...");
        CometMultiplePepXMLFileConverter cometMultiplePepXMLConverter = new CometMultiplePepXMLFileConverter(params.getCometPsmDir(), params.getCometPsmRegExp());
        cometMultiplePepXMLConverter.run();


        System.out.println("RunTime after Psm parsing: " + RunTime2String.getRuntimeString(Runtime.getRuntime()));

        ConcurrentHashMap<String, List<PeptideSpectrumMatch>>  psms = cometMultiplePepXMLConverter.getPsms();

        int cnt = 0;

        System.out.println(stringFunction.getHeader());
        for (String specID : psms.keySet()) {
            if(maxNrDisplayedPSMs>=0 && cnt>=maxNrDisplayedPSMs) break;

            for (PeptideSpectrumMatch psm : psms.get(specID)) {

                if(maxNrDisplayedPSMs>=0 && cnt>=maxNrDisplayedPSMs) break;

                System.out.println(stringFunction.getTabString(specID,psm));

                cnt++;
            }
        }

        return 0;
    }


    protected void createOptions() {

        this.cmdLineOpts = new Options();

        cmdLineOpts.addOption(Option.builder("coD").required().hasArg().longOpt("cometPsmDir").desc("Comet psm root directory (required)").build());
        cmdLineOpts.addOption(Option.builder("coRE").required().hasArg().longOpt("cometPsmRegex").desc("Regular expression of Comet psm files (e.g. \\.xml$)").build());
        cmdLineOpts.addOption(Option.builder("maxR").required(false).hasArg().longOpt("maxRank").desc("Maximal rank of peptide in list of spectrum matches (rank 1 = best) (default value: 1)").build());
        cmdLineOpts.addOption(Option.builder("minZ").required(false).hasArg().longOpt("minCharge").desc("Minimal charge of PSM (default value: 1)").build());
        cmdLineOpts.addOption(Option.builder("maxZ").required(false).hasArg().longOpt("maxCharge").desc("Maximal charge of PSM (default value: 5)").build());
        cmdLineOpts.addOption(Option.builder("minL").required(false).hasArg().longOpt("minLength").desc("Minimal length of peptide (default value: 8)").build());
        cmdLineOpts.addOption(Option.builder("maxL").required(false).hasArg().longOpt("maxLength").desc("Maximal length of peptide (default value: 25)").build());
        cmdLineOpts.addOption(Option.builder("maxP").required(false).hasArg().longOpt("maxDisplayedPsms").desc("Maximal number of psms written to standard output").build());
        cmdLineOpts.addOption(Option.builder("mod").required(false).hasArg().longOpt("modifications").desc("Comma separated list of peptide modifications used in search (e.g. Cysteinyl:C3H5NO2S,Oxidation:O)").build());
        cmdLineOpts.addOption(Option.builder("spRE").required(false).hasArg().longOpt("spectrumFilter").desc("If this option is set, only spectrum ids that match this regexp are used.  If not set no filtering is performed.").build());
        cmdLineOpts.addOption(Option.builder("h").required(false).hasArg(false).longOpt("help").desc("Help option for command line help").build());
        cmdLineOpts.addOption(Option.builder("v").required(false).hasArg(false).longOpt("version").desc("Version of NewAnce software").build());
    }

    @Override
    protected void check(CommandLine line) throws ParseException {

        this.params = NewAnceParams.getInstance();

        params.add("cometPsmDir", getOptionString(line,"coD"));
        params.add("cometPsmRegExp", getOptionString(line,"coRE"));

        params.add("maxRank", getOptionString(line,"maxR"));
        params.add("minCharge", getOptionString(line,"minZ"));
        params.add("maxCharge", getOptionString(line,"maxZ"));
        params.add("minPeptideLength", getOptionString(line,"minL"));
        params.add("maxPeptideLength", getOptionString(line,"maxL"));
        params.add("spectrumRegExp", getOptionString(line,"spRE"));
        params.add("modifications", getOptionString(line,"mod"));

        String maxPStr = getOptionString(line,"maxP");
        maxNrDisplayedPSMs = (maxPStr.isEmpty())?-1:Integer.parseInt(maxPStr);

        params.finalize();


    }

}
