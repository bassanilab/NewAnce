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

import newance.mzjava.mol.Peptide;
import newance.mzjava.mol.modification.AbsoluteTolerance;
import newance.mzjava.mol.modification.PpmTolerance;
import newance.mzjava.ms.io.mgf.MgfReader;
import newance.mzjava.ms.io.mgf.MsConvertTitleParser;
import newance.mzjava.ms.peaklist.IntervalList;
import newance.mzjava.ms.peaklist.PeakAnnotation;
import newance.mzjava.ms.peaklist.peakfilter.MzRangeFilter;
import newance.mzjava.ms.peaklist.peakfilter.NPeaksPerSlidingWindowFilter;
import newance.mzjava.ms.peaklist.peaktransformer.HighestPeakPerBinNormalizer;
import newance.mzjava.ms.peaklist.peaktransformer.LogTransformer;
import newance.mzjava.ms.spectrasim.MatchedNdpSimFunc;
import newance.mzjava.ms.spectrasim.NdpSimFunc;
import newance.mzjava.ms.spectrum.MsnSpectrum;
import newance.util.ExecutableOptions;
import org.apache.commons.cli.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Markus Müller
 */

public class CalcSpectrumSimilarity2 extends ExecutableOptions {

    private File mgfFile;
    private File outputFile;
    private final List<MsnSpectrum> spectrumList;


    public CalcSpectrumSimilarity2(String version) {

        this.version = version;
        this.spectrumList = new ArrayList<>();
        createOptions();
    }

    @Override
    protected void createOptions() {

        this.cmdLineOpts = new Options();

        cmdLineOpts.addOption(Option.builder("mgf").required().hasArg().longOpt("mgfFile").desc("mgf file containing spectra (required)").build());
        cmdLineOpts.addOption(Option.builder("o").required().hasArg().longOpt("simFile").desc("output file name for similarities (required)").build());
        cmdLineOpts.addOption(Option.builder("h").required(false).hasArg(false).longOpt("help").desc("Help option for command line help").build());
        cmdLineOpts.addOption(Option.builder("v").required(false).hasArg(false).longOpt("version").desc("Version of NewAnce software").build());
    }

    @Override
    protected void check(CommandLine line) throws ParseException {

        if (optionsSet) return;

        mgfFile = new File(getOptionString(line,"mgf"));
        outputFile = new File(getOptionString(line,"o"));
    }

    @Override
    public int run() throws IOException {

        readSpectra();
        calcCorrelation();

        return 0;
    }

    private MsnSpectrum discardPrecB1(MsnSpectrum spectrum) {
        IntervalList intervals = new IntervalList();

        intervals.addInterval(100.0,150.0);
        intervals.addInterval(440.0,450.0);

        MzRangeFilter<PeakAnnotation> mzRangeFilter = new MzRangeFilter<PeakAnnotation>(intervals, false);
        spectrum.apply(mzRangeFilter);

        return spectrum;
    }


    private MsnSpectrum transform(MsnSpectrum spectrum) {
//        LogTransformer<PeakAnnotation> processor = new LogTransformer<PeakAnnotation>();
        NPeaksPerSlidingWindowFilter<PeakAnnotation> processor =
                new NPeaksPerSlidingWindowFilter<PeakAnnotation>(2, 5.0, 0.1);;
        spectrum.apply(processor);

        return spectrum;
    }


    private void readSpectra() {

        try {

            System.out.println("Parsing and matching mgf file : " + mgfFile.getAbsolutePath());

            MgfReader reader = new MgfReader(mgfFile, new MsConvertTitleParser());

            MsnSpectrum spectrum = null;
            while (reader.hasNext()) {
                try {
                    spectrum = reader.next();
                    spectrum = transform(spectrum);
                    this.spectrumList.add(spectrum);

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

        FileWriter fileWriter = new FileWriter(outputFile);
        for (MsnSpectrum spectrum : spectrumList) {
            fileWriter.write("\t");
            fileWriter.write(spectrum.getComment());
        }
        fileWriter.write("\n");

        MatchedNdpSimFunc matchedNdpSimFunc = new MatchedNdpSimFunc(5, new PpmTolerance(20.0));

        NdpSimFunc ndpSimFunc = new NdpSimFunc(5, new PpmTolerance(20.0));

        double[][] sims = new double[spectrumList.size()][spectrumList.size()];
        for (int i=0; i<spectrumList.size(); i++) {
            sims[i][i] = 1.0;
            for (int j=i+1; j<spectrumList.size(); j++) {
                sims[i][j] = matchedNdpSimFunc.calcSimilarity(spectrumList.get(i),spectrumList.get(j));
                sims[j][i] = sims[i][j];
            }
        }

        for (int i=0; i<spectrumList.size(); i++) {
            fileWriter.write(spectrumList.get(i).getComment()+"\t");
            for (int j=0; j<spectrumList.size(); j++) {
                String tab = (j==0)?"":"\t";
                fileWriter.write(String.format("%s%f", tab, sims[i][j]));
            }
            fileWriter.write("\n");
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

        CalcSpectrumSimilarity2 createSpectrumLibrary =  new CalcSpectrumSimilarity2("Version 1.0.0");
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

