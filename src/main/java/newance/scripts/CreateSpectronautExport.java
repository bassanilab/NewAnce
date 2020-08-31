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

import newance.util.ExecutableOptions;
import newance.util.NewAnceParams;
import org.apache.commons.cli.*;
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

    private File maxquantPsmDir;
    private File newAnceResultFile;
    private List<String> mqFeatures;
    private List<String> spectronautFeatures;
    private Set<String> newAnceSpectra;
    private Map<String,String> maxQuantLines;

    public CreateSpectronautExport(String version) {

        this.version = version;
        createOptions();
    }

    @Override
    protected void createOptions() {

        this.cmdLineOpts = new Options();

        cmdLineOpts.addOption(Option.builder("mqD").required(false).hasArg().longOpt("maxquantPsmDir").desc("MaxQuant psm root directory.").hasArg().build());
        cmdLineOpts.addOption(Option.builder("naf").required().hasArg().longOpt("newAnceResultFile").desc("Result file from NewAnce analysis (required)").build());
        cmdLineOpts.addOption(Option.builder("h").required(false).hasArg(false).longOpt("help").desc("Help option for command line help").build());
        cmdLineOpts.addOption(Option.builder("v").required(false).hasArg(false).longOpt("version").desc("Version of NewAnce software").build());
    }

    @Override
    protected void check(CommandLine line) throws ParseException {

        if (optionsSet) return;

        newAnceResultFile = new File(NewAnceParams.getFileValue("newAnceResultFile",getOptionString(line, "naf")));
        maxquantPsmDir = new File(NewAnceParams.getDirectoryValue("maxquantPsmDir",getOptionString(line, "mqD")));

        spectronautFeatures = new ArrayList<>();
        mqFeatures = new ArrayList<>();

        spectronautFeatures.add("Raw File");
        mqFeatures.add("Raw file");

        spectronautFeatures.add("Precursor Charge");
        mqFeatures.add("Charge");

        spectronautFeatures.add("Stripped Sequence");
        mqFeatures.add("Sequence");

        spectronautFeatures.add("Protein Group Id");
        mqFeatures.add("Proteins");

        spectronautFeatures.add("Labeled Sequence");
        mqFeatures.add("Modified sequence");

        spectronautFeatures.add("Scan Event");
        mqFeatures.add("Scan event number");

        spectronautFeatures.add("Scan Number");
        mqFeatures.add("Scan number");

        spectronautFeatures.add("Retention Time");
        mqFeatures.add("Retention time");

        spectronautFeatures.add("MS1 Intensity");
        mqFeatures.add("Precursor Intensity");
    }

    @Override
    public int run() throws IOException {

        readNewAnceFile();
        readMaxQuantFiles();

        writeSpectronautFile();

        return 0;
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

        maxQuantLines = new HashMap<>();

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


            for (String feature : mqFeatures) {
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
                    for (String feature : mqFeatures) {
                        if (results.findColumn(feature) != 0)
                            values += values.isEmpty()?results.getString(feature):"\t"+results.getString(feature);
                        else
                            values += "\tNA";
                    }
                    maxQuantLines.put(specID,values);
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

        newAnceSpectra = new HashSet<>();
        try {
            BufferedReader lineReader = new BufferedReader(new FileReader(newAnceResultFile ));

            String line;
            while ((line = lineReader.readLine()) != null) {

                if (line.isEmpty() || line.startsWith("Spectrum\tScanNr\tCharge")) continue;

                String[] fields = line.split("\t");

                newAnceSpectra.add(fields[0]);

            }
            lineReader.close();
        } catch (FileNotFoundException e) {
        } catch (IOException e) {
        }
    }

    private void writeSpectronautFile() {
        try {

            String newNewAnceFile = newAnceResultFile.getAbsolutePath().replace("PSMs.txt","_Spectronaut.txt");
            BufferedWriter writer = new BufferedWriter(new FileWriter(newNewAnceFile));

            String header = "";
            for (String feature : spectronautFeatures) {
                header += (header.isEmpty())?feature:"\t"+feature;
            }
            writer.write(header+"\n");

            for  (String key : newAnceSpectra) {

                writer.write(maxQuantLines.get(key)+"\n");
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

        CreateSpectronautExport addMaxQuantFeatures =  new CreateSpectronautExport("Version 1.0.0");
        try {
            addMaxQuantFeatures.init(args).parseOptions(args).run();
        } catch (MissingOptionException e) {
        } catch (ParseException e) {
            addMaxQuantFeatures.printOptions(args, e.getMessage());
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
