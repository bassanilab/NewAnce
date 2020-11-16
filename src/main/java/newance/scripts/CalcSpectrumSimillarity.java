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

import newance.mzjava.mol.*;
import newance.mzjava.mol.modification.*;
import newance.mzjava.ms.consensus.PeptideConsensusSpectrum;
import newance.mzjava.ms.io.mgf.MgfReader;
import newance.mzjava.ms.io.mgf.MgfWriter;
import newance.mzjava.ms.io.mgf.MsConvertTitleParser;
import newance.mzjava.ms.peaklist.PeakList;
import newance.mzjava.ms.peaklist.peakfilter.AbstractMergePeakFilter;
import newance.mzjava.ms.spectrasim.MatchedNdpSimFunc;
import newance.mzjava.ms.spectrasim.NdpSimFunc;
import newance.mzjava.ms.spectrum.MsnSpectrum;
import newance.mzjava.ms.spectrum.PepFragAnnotation;
import newance.util.ExecutableOptions;
import org.apache.commons.cli.*;
import org.apache.commons.io.FilenameUtils;

import java.io.*;
import java.net.URLEncoder;
import java.sql.*;
import java.util.*;

/**
 * @author Markus Müller
 */

public class CalcSpectrumSimillarity extends ExecutableOptions {

    private File mgfFile;
    private final Map<String,MsnSpectrum> consSpectrumMap;
    private final Map<String,MsnSpectrum> prositSpectrumMap;

    public CalcSpectrumSimillarity(String version) {

        this.version = version;
        this.consSpectrumMap = new HashMap<>();
        this.prositSpectrumMap = new HashMap<>();
        createOptions();
    }

    @Override
    protected void createOptions() {

        this.cmdLineOpts = new Options();

        cmdLineOpts.addOption(Option.builder("mgf").required().hasArg().longOpt("mgfDir").desc("Directory containing mgf files (required)").build());
        cmdLineOpts.addOption(Option.builder("h").required(false).hasArg(false).longOpt("help").desc("Help option for command line help").build());
        cmdLineOpts.addOption(Option.builder("v").required(false).hasArg(false).longOpt("version").desc("Version of NewAnce software").build());
    }

    @Override
    protected void check(CommandLine line) throws ParseException {

        if (optionsSet) return;

        mgfFile = new File(getOptionString(line,"mgf"));
    }

    @Override
    public int run() throws IOException {

        readSpectra();
        calcCorrelation();

        return 0;
    }


    private void readSpectra() {

        try {

            System.out.println("Parsing and matching mgf file : " + mgfFile.getAbsolutePath());

            MgfReader reader = new MgfReader(mgfFile, new MsConvertTitleParser());

            MsnSpectrum spectrum = null;
            while (reader.hasNext()) {
                try {
                    spectrum = reader.next();
                    String comment = spectrum.getComment();
                    if (comment.endsWith("cons")) {
                        comment = comment.replace(".cons","");
                        consSpectrumMap.putIfAbsent(comment,spectrum);
                    }
                    else
                        prositSpectrumMap.putIfAbsent(comment,spectrum);

                } catch (IOException e) {
                    System.out.println("WARNING : " + e.getMessage());
                }
            }
        }
        catch(IOException e) {
            System.out.println("Cannot parse mgf file : " + mgfFile.getAbsolutePath());
        }
    }

    private void calcCorrelation() throws IOException {

        String simFileName = mgfFile.getAbsolutePath();
        simFileName = simFileName.replace(".mgf","_similarity.txt");

        FileWriter fileWriter = new FileWriter(new File(simFileName));
        fileWriter.write("peptide\tcharge\tpep_mass\tlength\tmodification\tmatched_similarity\tall_similarity\n");

        for (String peptideKey: prositSpectrumMap.keySet()) {
            if (consSpectrumMap.containsKey(peptideKey)) {
                MsnSpectrum prositSpectrum = prositSpectrumMap.get(peptideKey);
                MsnSpectrum consSpectrum = consSpectrumMap.get(peptideKey);

                MatchedNdpSimFunc matchedNdpSimFunc = new MatchedNdpSimFunc(5, new AbsoluteTolerance(0.1));

                double matchedSim = matchedNdpSimFunc.calcSimilarity(prositSpectrum,consSpectrum);

                NdpSimFunc ndpSimFunc = new NdpSimFunc(5, new AbsoluteTolerance(0.1));

                double allSim = ndpSimFunc.calcSimilarity(prositSpectrum,consSpectrum);

                String[] fields = peptideKey.split("\\.");
                Peptide peptide = Peptide.parse(fields[0]);

                fileWriter.write(peptide.toString()+"\t"+prositSpectrum.getPrecursor().getCharge()+"\t"+
                        prositSpectrum.getPrecursor().getMass()+"\t"+peptide.size()+"\t"+
                        peptide.hasModifications()+"\t"+matchedSim+"\t"+allSim+"\n");
            }
        }

        fileWriter.close();
    }

    public ExecutableOptions init(String[] args) throws IOException {

        optionsSet = false;
        checkHelpOption(args, "-h");
        checkVersionOption(args, version, "-v");

        return this;
    }

    public static void main(String[] args) {

        CalcSpectrumSimillarity createSpectrumLibrary =  new CalcSpectrumSimillarity("Version 1.0.0");
        try {
            createSpectrumLibrary.init(args).parseOptions(args).run();
        } catch (MissingOptionException e) {
        } catch (ParseException e) {
            createSpectrumLibrary.printOptions(args, e.getMessage());
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}

