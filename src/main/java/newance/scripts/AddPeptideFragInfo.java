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

import java.util.Map.Entry;
import com.google.common.collect.Lists;
import newance.mzjava.mol.*;
import newance.mzjava.mol.modification.AbsoluteTolerance;
import newance.mzjava.ms.consensus.PeptideFragmentCoverage;
import newance.mzjava.ms.io.mgf.MgfReader;
import newance.mzjava.ms.io.mgf.MgfWriter;
import newance.mzjava.ms.io.mgf.MsConvertTitleParser;
import newance.mzjava.ms.peaklist.DoublePeakList;
import newance.mzjava.ms.peaklist.PeakList;
import newance.mzjava.ms.peaklist.peakfilter.NPeaksPerSlidingWindowFilter;
import newance.mzjava.ms.spectrum.MsnSpectrum;
import newance.mzjava.ms.spectrum.PepFragAnnotation;
import newance.util.ExecutableOptions;
import newance.util.StringFileWriter;
import org.apache.commons.cli.*;
import org.apache.commons.io.FilenameUtils;

import java.io.*;
import java.net.URLEncoder;
import java.sql.*;
import java.util.*;

/**
 * @author Markus Müller
 */

public class AddPeptideFragInfo extends ExecutableOptions {

    private class StringPair implements Entry<String, String> {

        String key;
        String value;

        StringPair(String key, String value) {
            this.key = key;
            this.value = value;
        }

        @Override
        public String getKey() {
            return key;
        }

        @Override
        public String getValue() {
            return value;
        }

        @Override
        public String setValue(String value) {
            this.value = value;
            return value;
        }
    }

    private File mgfDir;
    private File newAnceResultFile;
    private Set<IonType> ionTypes;
    private boolean includeH2Oloss;
    private boolean includeH3Nloss;
    private int nrPeaksPerWindow;
    private final Set<String> peptides;
    private final Map<String,Double> bestScoringPSM;
    private final Map<String,Map<String,Peptide>> spectrumPSMMap;
    private final Map<String,StringPair> peptideSpectrumMap;
    private final Map<String,File> mgfFileMap;
    private PeptideFragmentAnnotator annotator;

    public AddPeptideFragInfo(String version) {

        this.version = version;
        this.spectrumPSMMap = new HashMap<>();
        this.mgfFileMap = new HashMap<>();
        this.bestScoringPSM = new HashMap<>();
        this.peptideSpectrumMap = new HashMap<>();
        this.peptides = new HashSet();
        this.nrPeaksPerWindow = -1;


        createOptions();
    }

    @Override
    protected void createOptions() {

        this.cmdLineOpts = new Options();

        cmdLineOpts.addOption(Option.builder("mgf").required().hasArg().longOpt("mgfDir").desc("Directory containing mgf files (required)").build());
        cmdLineOpts.addOption(Option.builder("naf").required().hasArg().longOpt("newAnceResultFile").desc("Result file from NewAnce analysis (required)").build());
        cmdLineOpts.addOption(Option.builder("wl").longOpt("waterloss").desc("Include water loss on S, T, D, E").build());
        cmdLineOpts.addOption(Option.builder("al").longOpt("ammoniumloss").desc("Include ammonium loss on Q, K, R, N").build());
        cmdLineOpts.addOption(Option.builder("it").hasArg().longOpt("iontypes").desc("Which fragment ions to include (default: aby)").build());
        cmdLineOpts.addOption(Option.builder("np").hasArg().longOpt("nrPeaksWindow").desc("Retain only np highest peaks in 20Da window").build());
        cmdLineOpts.addOption(Option.builder("p").hasArg().longOpt("peptides").desc("Annotations for all these peptides will be provided. " +
                "Otherwise all peptides in the newance file will be processed.").build());
        cmdLineOpts.addOption(Option.builder("h").required(false).hasArg(false).longOpt("help").desc("Help option for command line help").build());
        cmdLineOpts.addOption(Option.builder("v").required(false).hasArg(false).longOpt("version").desc("Version of NewAnce software").build());
    }

    @Override
    protected void check(CommandLine line) throws ParseException {

        if (optionsSet) return;

        newAnceResultFile = new File(getOptionString(line,"naf"));
        mgfDir = new File(getOptionString(line,"mgf"));
        includeH2Oloss = checkBooleanOption(line,"wl");
        includeH3Nloss = checkBooleanOption(line,"al");

        String npStr = getOptionString(line,"np");
        if (!npStr.isEmpty()) {
            nrPeaksPerWindow = Integer.parseInt(npStr);
        }


        String itStr = getOptionString(line,"it");
        if (itStr.isEmpty()) {
            ionTypes = EnumSet.noneOf(IonType.class);
            ionTypes.add(IonType.a);
            ionTypes.add(IonType.b);
            ionTypes.add(IonType.y);
        } else if (!(itStr.contains("a") || itStr.contains("b") || itStr.contains("y"))) {
            throw new ParseException("No valid ion types in -it option: "+ itStr);
        } else {

            ionTypes = EnumSet.noneOf(IonType.class);
            for (char it : itStr.toCharArray()) {
                if (it=='a') ionTypes.add(IonType.a);
                if (it=='b') ionTypes.add(IonType.b);
                if (it=='y') ionTypes.add(IonType.y);
            }
        }

        String pepStr = getOptionString(line,"p");
        for (String p : pepStr.split(",")) {
            peptides.add(p);
        }

    }

    @Override
    public int run() throws IOException {

        setFragmentAnnotator();
        findMgfFiles();
        readNewAnceFile();
        createSpectrumPSMMap();

        System.out.println("Spectrum\tPeptide\tBestSpScore\tSequenceCoverage\tSpectrumCoverage\t" +
                "PeptideAnnotShort\tPeptideAnnot\tPeptideAnnotLong");
        for (String mgfName : spectrumPSMMap.keySet()) {
            annotateMatchedPeptides(mgfName, spectrumPSMMap.get(mgfName));
        }

        return 0;
    }

    private void createSpectrumPSMMap() {

        for (String p : peptideSpectrumMap.keySet()) {
            Peptide peptide = Peptide.parse(p);

            Entry<String,String> spectrumInfo = peptideSpectrumMap.get(p);
            spectrumPSMMap.putIfAbsent(spectrumInfo.getKey(),new HashMap<>());
            spectrumPSMMap.get(spectrumInfo.getKey()).put(spectrumInfo.getValue(),peptide);
        }
    }

    private void setFragmentAnnotator() {

        List<PeptidePeakGenerator<PepFragAnnotation>> peakGenerators = Lists.newArrayList();

        peakGenerators.add(new BackbonePeakGenerator(ionTypes, 1));

        if (includeH2Oloss) {
            Mass waterLoss = Composition.parseComposition("H-2O-1");

            peakGenerators.add(new PeptideNeutralLossPeakGenerator(waterLoss,
                    EnumSet.of(AminoAcid.S, AminoAcid.T, AminoAcid.D, AminoAcid.E), ionTypes, 1));
        }

        if (includeH3Nloss) {
            Mass ammoniumLoss = Composition.parseComposition("H-3N-1");

            peakGenerators.add(new PeptideNeutralLossPeakGenerator(ammoniumLoss,
                    EnumSet.of(AminoAcid.Q, AminoAcid.K, AminoAcid.R, AminoAcid.N), ionTypes, 1));
        }

        PeptideFragmenter fragmenter = new PeptideFragmenter(peakGenerators, PeakList.Precision.DOUBLE);

        this.annotator = new PeptideFragmentAnnotator(fragmenter, new AbsoluteTolerance(0.1));
    }

    private void findMgfFiles() {

        FilenameFilter mgfFilter = new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {

                return name.toLowerCase().endsWith(".mgf");
            }
        };

        for(File mgfFile : mgfDir.listFiles(mgfFilter)) {
            mgfFileMap.put(mgfFile.getName().substring(0,mgfFile.getName().indexOf('.')),mgfFile);
        }
    }

    private void annotateMatchedPeptides(String rawFileName, Map<String,Peptide> psms) {

        File mgfFile = new File(mgfDir.getAbsolutePath()+File.separator+rawFileName+".mgf");

        try {

            MgfReader reader = new MgfReader(mgfFile, new MsConvertTitleParser());

            MsnSpectrum spectrum = null;
            while (reader.hasNext()) {
                try {
                    spectrum = reader.next();

                    if (psms.containsKey(spectrum.getComment())) {
                        annotateMatchedPeptide(spectrum,  psms.get(spectrum.getComment()));
                    }

                } catch (IOException e) {
                    System.out.println("WARNING : " + e.getMessage());
                }
            }
        }
        catch(IOException e) {
            System.out.println("Cannot parse mgf file : " + mgfFile.getAbsolutePath());
        }
    }

    private void annotateMatchedPeptide(MsnSpectrum spectrum,Peptide peptide) {

        PeakList<PepFragAnnotation> peakList = new DoublePeakList<>(1);
        peakList.addPeaksNoAnnotations(spectrum);
        peakList.setPrecursor(spectrum.getPrecursor());

        if (nrPeaksPerWindow>0)
            peakList.apply(new NPeaksPerSlidingWindowFilter<PepFragAnnotation>(nrPeaksPerWindow,  20.0, 1.2));

        annotator.annotate(peakList, peptide);

        String peptStr = peptide.toString();
        PeptideFragmentCoverage peptideFragmentCoverage  = new PeptideFragmentCoverage(peakList, peptide);

        String seqcov = String.format("%.4f",peptideFragmentCoverage.getSequenceCoverage());
        String speccov = String.format("%.4f",peptideFragmentCoverage.getSpectrumCoverage());
        System.out.println(spectrum.getComment()+"\t"+peptStr+"\t"+bestScoringPSM.get(peptStr)+"\t"+
                seqcov+"\t"+speccov+"\t"+ peptideFragmentCoverage.getPeptideFragInfoShort()+"\t"+
                peptideFragmentCoverage.getPeptideFragInfo()+"\t"+peptideFragmentCoverage.getPeptideFragInfoLong());
    }

    private void readNewAnceFile() {

        try {
            // load the driver into memory
            Class.forName("org.relique.jdbc.csv.CsvDriver");

            Properties props = new Properties();
            props.put("fileExtension", "." + FilenameUtils.getExtension(newAnceResultFile.getName()));
            // create a connection. The first command line parameter is assumed to
            //  be the directory in which the .csv files are held
            Connection conn = DriverManager.getConnection("jdbc:relique:csv:" + newAnceResultFile.getParent() + "?separator=" + URLEncoder.encode("\t", "UTF-8"), props);

            // create a Statement object to execute the query with
            Statement stmt = conn.createStatement();

            // Select the ID and NAME columns from sample.csv
            ResultSet results = stmt.executeQuery("SELECT * FROM " + FilenameUtils.removeExtension(newAnceResultFile.getName()));

            String annotFileName = newAnceResultFile.getAbsolutePath();
            annotFileName = annotFileName.replaceAll(".txt$","_FragAnnot.txt");

//            FileWriter fileWriter = new FileWriter(new File(annotFileName));
//            fileWriter.write("spectrum_title\tpeptide\tcharge\tpep_mass\tmz\tmodification\n");

            while (results.next()) {

                String spectrum = results.getString("Spectrum");
                String peptideStr = results.getString("Peptide");
                double spscore = results.getDouble("SpScore");
                String msmsFileName = spectrum.substring(0,spectrum.indexOf('.'));

                if (!peptides.isEmpty() && !peptides.contains(peptideStr)) continue;

                if (!bestScoringPSM.containsKey(peptideStr) || spscore > bestScoringPSM.get(peptideStr)) {

                    bestScoringPSM.put(peptideStr, spscore);
                    peptideSpectrumMap.put(peptideStr,new StringPair(msmsFileName,spectrum));
                }
             }

            results.close();
            stmt.close();
            conn.close();
        } catch (ClassNotFoundException | SQLException | UnsupportedEncodingException e) {
            throw new IllegalStateException(e);
        } catch (IOException e) {
            System.out.println("Cannot parse NewAnce file : " + mgfDir.getAbsolutePath());
        }
    }

    public ExecutableOptions init(String[] args) throws IOException {

        optionsSet = false;
        checkHelpOption(args, "-h");
        checkVersionOption(args, version, "-v");

        return this;
    }

    public static void main(String[] args) {

        AddPeptideFragInfo addPeptideFragInfo =  new AddPeptideFragInfo("Version 1.0.0");
        try {
            addPeptideFragInfo.init(args).parseOptions(args).run();
        } catch (MissingOptionException e) {
        } catch (ParseException e) {
            addPeptideFragInfo.printOptions(args, e.getMessage());
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}

