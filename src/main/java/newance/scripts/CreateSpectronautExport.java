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

import newance.mzjava.ms.io.mgf.MgfReader;
import newance.mzjava.ms.io.mgf.MsConvertTitleParser;
import newance.mzjava.ms.spectrum.MsnSpectrum;
import newance.util.ExecutableOptions;
import newance.util.NewAnceParams;
import org.apache.commons.cli.*;
import org.apache.commons.collections15.map.HashedMap;
import org.apache.commons.io.FilenameUtils;

import java.io.*;
import java.net.URLEncoder;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.*;
import java.util.*;
import java.util.regex.Pattern;

/**
 * Created by markusmueller on 17.02.20.
 */
public class CreateSpectronautExport extends ExecutableOptions {

    private File mgfDir;
    private File maxquantPsmDir;
    private File newAnceResultFile;
    private List<String> psmFeatures;
    private List<String> spectronautFeatures;
    private Set<String> newAnceSpectra;
    private Map<String,String> psmLines;
    private Map<String,File> mgfFileMap;
    private Map<String,Integer> spectrumScanEventNrMap;
    private Map<String,Set<Integer>> spectrumScanNrMap;
    private Map<String,Double> spectrumIntensityMap;

    public CreateSpectronautExport(String version) {

        this.version = version;
        createOptions();
    }

    @Override
    protected void createOptions() {

        this.cmdLineOpts = new Options();

        cmdLineOpts.addOption(Option.builder("mgf").required(false).hasArg().longOpt("mgfDir").desc("Directory containing mgf files (for peptide intensities)").build());
        cmdLineOpts.addOption(Option.builder("mqD").required(false).hasArg().longOpt("maxquantPsmDir").desc("MaxQuant psm root directory (for peptide intensities).").hasArg().build());
        cmdLineOpts.addOption(Option.builder("naf").required().hasArg().longOpt("newAnceResultFile").desc("Result file from NewAnce analysis (required)").build());
        cmdLineOpts.addOption(Option.builder("h").required(false).hasArg(false).longOpt("help").desc("Help option for command line help").build());
        cmdLineOpts.addOption(Option.builder("v").required(false).hasArg(false).longOpt("version").desc("Version of NewAnce software").build());
    }

    @Override
    protected void check(CommandLine line) throws ParseException {

        if (optionsSet) return;

        newAnceResultFile = new File(NewAnceParams.getFileValue("newAnceResultFile",getOptionString(line, "naf")));
        mgfDir = getOptionString(line, "mgf").isEmpty()?null:new File(NewAnceParams.getDirectoryValue("mgfDir",getOptionString(line, "mgf")));
        maxquantPsmDir = getOptionString(line, "mqD").isEmpty()?null:new File(NewAnceParams.getDirectoryValue("maxquantPsmDir",getOptionString(line, "mqD")));

        spectronautFeatures = new ArrayList<>();
        spectronautFeatures.add("Raw File");
        spectronautFeatures.add("Precursor Charge");
        spectronautFeatures.add("Stripped Sequence");
        spectronautFeatures.add("Protein Group Id");
        spectronautFeatures.add("Labeled Sequence");
        spectronautFeatures.add("Scan Event");
        spectronautFeatures.add("Scan Number");
        spectronautFeatures.add("Retention Time");
        spectronautFeatures.add("MS1 Intensity");

        if (maxquantPsmDir != null) {

            psmFeatures = new ArrayList<>();

            psmFeatures.add("Raw file");
            psmFeatures.add("Charge");
            psmFeatures.add("Sequence");
            psmFeatures.add("Proteins");
            psmFeatures.add("Modified sequence");
            psmFeatures.add("Scan event number");
            psmFeatures.add("Scan number");
            psmFeatures.add("Retention time");
            psmFeatures.add("Precursor Intensity");
        } else if (mgfDir != null) {

            psmFeatures = new ArrayList<>();

            psmFeatures.add("Spectrum");
            psmFeatures.add("Charge");
            psmFeatures.add("Sequence");
            psmFeatures.add("ProteinList");
            psmFeatures.add("Modified sequence");
            psmFeatures.add("Scan event number");
            psmFeatures.add("ScanNr");
            psmFeatures.add("RT");
            psmFeatures.add("Precursor Intensity");
        }
    }

    @Override
    public int run() throws IOException {

        readNewAnceFile();
        if (maxquantPsmDir != null) {
            readMaxQuantFiles();
        } else if (mgfDir != null) {
            findMgfFiles();

            spectrumScanEventNrMap = new HashMap<>();
            spectrumIntensityMap = new HashMap<>();

            for (String mgfFileName : mgfFileMap.keySet()) {
                parseMgfFile(mgfFileMap.get(mgfFileName), spectrumScanNrMap.get(mgfFileName));
            }
            rereadNewAnceFile();
        }

        writeSpectronautFile();

        return 0;
    }

    private void findMgfFiles() {

        mgfFileMap = new HashMap<>();

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


    private List<MsnSpectrum> parseMgfFile(File mgfFile, Set<Integer> scanNumbers) {

        List<MsnSpectrum> spectra = new ArrayList<>();

        System.out.println("Parsing and matching mgf file : " + mgfFile.getAbsolutePath());
        String fileName = mgfFile.getName();
        fileName = fileName.substring(0,fileName.lastIndexOf('.'));

        try {
            MgfReader reader = new MgfReader(mgfFile, new MsConvertTitleParser());

            MsnSpectrum spectrum = null;
            int prevScanNr = 0;
            int scanEventNr = 0;
            while (reader.hasNext()) {
                try {
                    spectrum = reader.next();

                    int scanNr = spectrum.getScanNumbers().getFirst().getValue();

                    if (scanNr == prevScanNr+1) scanEventNr++;
                    else scanEventNr = 1;

                    prevScanNr = scanNr;

                    if (scanNumbers.contains(scanNr)) {
                        spectrum.setComment(fileName+"."+scanNr+"."+scanNr+"."+spectrum.getPrecursor().getCharge());
                        spectra.add(spectrum);
                        spectrumScanEventNrMap.putIfAbsent(spectrum.getComment(),scanEventNr);
                        spectrumIntensityMap.putIfAbsent(spectrum.getComment(),spectrum.getPrecursor().getIntensity());
                    }
                } catch (IOException e) {
                    System.out.println("WARNING : " + e.getMessage());
                }
            }
        }
        catch(IOException e) {
            System.out.println("Cannot parse mgf file : " + mgfFile.getAbsolutePath());
        }

        return spectra;
    }


    private void readMaxQuantFiles() {

        try {
            final Pattern regex = NewAnceParams.getInstance().getMaxquantPsmRegExp();
            final List<File> psmFileList = new ArrayList<>();
            Files.walk(Paths.get(maxquantPsmDir.getAbsolutePath()))
                    .filter(Files::isRegularFile)
                    .filter(f -> regex.matcher(f.getFileName().toString()).find())
                    .forEach(f -> psmFileList.add(f.toFile()));

            for (File msmsFile : psmFileList) {
                readMaxQuantFile(msmsFile);
            }

        } catch (IOException e) {
        }
    }

    private void readMaxQuantFile(File file) {

        psmLines = new HashMap<>();

        try {
            // load the driver into memory
            Class.forName("org.relique.jdbc.csv.CsvDriver");

            Properties props = new Properties();
            props.put("fileExtension", "." + FilenameUtils.getExtension(file.getName()));
            // create a connection. The first command line parameter is assumed to
            //  be the directory in which the .csv files are held
            Connection conn = DriverManager.getConnection("jdbc:relique:csv:" + file.getParent() + "?separator=" + URLEncoder.encode("\t", "UTF-8"), props);

            // create a Statement object to execute the query with
            Statement stmt = conn.createStatement();

            // Select the ID and NAME columns from sample.csv
            ResultSet results = stmt.executeQuery("SELECT * FROM " + FilenameUtils.removeExtension(file.getName()));


            for (String feature : psmFeatures) {
                if (results.findColumn(feature)==0) {
                    System.out.println("ERROR: feature "+feature+" is not a valid msms.txt column name. This feature is ignored.");
                }
            }

            while (results.next()) {

                int charge = results.getInt("Charge");
                int scanNumber = results.getInt("Scan number");
                String rawFile = results.getString("Raw file");

                String specID = rawFile + "." + scanNumber + "." + scanNumber + "." + charge;

                if (newAnceSpectra.contains(specID)) {

                    String values = "";
                    for (String feature : psmFeatures) {
                        if (results.findColumn(feature) != 0)
                            values += values.isEmpty()?results.getString(feature):"\t"+results.getString(feature);
                        else
                            values += "\tNA";
                    }
                    psmLines.put(specID,values);
                }
            }

            // clean up
            results.close();
            stmt.close();
            conn.close();
        } catch (ClassNotFoundException | SQLException | UnsupportedEncodingException e) {
            throw new IllegalStateException(e);
        }
    }

    private void readNewAnceFile() {

        spectrumScanNrMap = new HashMap<>();
        newAnceSpectra = new HashSet<>();
        try {
            BufferedReader lineReader = new BufferedReader(new FileReader(newAnceResultFile ));

            String line;
            while ((line = lineReader.readLine()) != null) {

                if (line.isEmpty() || line.startsWith("Spectrum\tScanNr\tCharge")) continue;

                String[] fields = line.split("\t");

                newAnceSpectra.add(fields[0]);

                String msmsFileName = fields[0].substring(0,fields[0].indexOf('.'));

                spectrumScanNrMap.putIfAbsent(msmsFileName,new HashSet<>());
                spectrumScanNrMap.get(msmsFileName).add(Integer.parseInt(fields[1]));
            }
            lineReader.close();

        } catch (FileNotFoundException e) {
        } catch (IOException e) {
        }
    }

    private String change2SpectronautFormat(String peptide) {

        if (peptide.indexOf('(')<0) return "_"+peptide+"_";

        String modifiedPeptide = peptide.replaceAll("Oxidation","Oxidation (M)");
        modifiedPeptide = modifiedPeptide.replaceAll("Phopsho","Phospho (STY)");
        modifiedPeptide = modifiedPeptide.replaceAll("\\(Acetyl\\)_","(Acetyl (Protein N-term))");

        modifiedPeptide = "_"+modifiedPeptide+"_";
        return modifiedPeptide;
    }

    private void rereadNewAnceFile() {

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
            ResultSet results =
                    stmt.executeQuery("SELECT * FROM " + FilenameUtils.removeExtension(newAnceResultFile.getName()));


            psmLines = new HashMap<>();

            while (results.next()) {

                String spectrum = results.getString("Spectrum");

                String values = "";
                for (String feature : psmFeatures) {
                    if (!values.isEmpty()) values += "\t";

                    if (feature.startsWith("ProteinList")) {
                        String proteins = results.getString("Proteins");
                        proteins = proteins.replace("[","");
                        proteins = proteins.replace("]","");
                        proteins = proteins.replaceAll(" ","");
                        proteins = proteins.replaceAll(";",",");
                        values += proteins;
                    } else if (feature.startsWith("Modified sequence")) {
                        String modifiedPeptide = change2SpectronautFormat(results.getString("Peptide"));
                        values += modifiedPeptide;
                    } else if (feature.startsWith("Scan event number")) {
                        if (spectrumScanEventNrMap.containsKey(spectrum))
                            values += spectrumScanEventNrMap.get(spectrum);
                        else
                            values += "NA";

                    } else if (feature.startsWith("Precursor Intensity")) {
                        if (spectrumIntensityMap.containsKey(spectrum))
                            values += spectrumIntensityMap.get(spectrum);
                        else
                            values += "NA";

                    } else if (results.findColumn(feature) != 0) {
                        values += results.getString(feature);
                    }
                    else {
                        values += "NA";
                    }
                }

                psmLines.put(spectrum,values);
            }

            // clean up
            results.close();
            stmt.close();
            conn.close();
        } catch (ClassNotFoundException | SQLException | UnsupportedEncodingException e) {
            throw new IllegalStateException(e);
        } catch (IOException e) {
            System.out.println("Cannot parse NewAnce file : " + mgfDir.getAbsolutePath());
        }
    }

    private void writeSpectronautFile() {
        try {

            String newNewAnceFile = newAnceResultFile.getAbsolutePath().replace(".txt","_Spectronaut.txt");
            BufferedWriter writer = new BufferedWriter(new FileWriter(newNewAnceFile));

            String header = "";
            for (String feature : spectronautFeatures) {
                header += (header.isEmpty())?feature:"\t"+feature;
            }
            writer.write(header+"\n");

            for  (String key : newAnceSpectra) {

                writer.write(psmLines.get(key)+"\n");
            }
            writer.close();

        } catch (IOException e) {
        }
    }

    public ExecutableOptions init(String[] args) throws IOException {

        optionsSet = false;
        checkHelpOption(args, "-h");
        checkVersionOption(args, version, "-v");

        return this;
    }

    public static void main(String[] args) {

        CreateSpectronautExport createSpectronautExport =  new CreateSpectronautExport("Version 1.0.0");
        try {
            createSpectronautExport.init(args).parseOptions(args).run();
        } catch (MissingOptionException e) {
        } catch (ParseException e) {
            createSpectronautExport.printOptions(args, e.getMessage());
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
