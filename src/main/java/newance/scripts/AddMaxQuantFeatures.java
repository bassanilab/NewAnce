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
public class AddMaxQuantFeatures extends ExecutableOptions {

    private File maxquantPsmDir;
    private File newAnceResultFile;
    private Set<String> maxQuantFeatures;
    private Map<String,String> newAnceLines;
    private Map<String,String> maxQuantLines;

    public AddMaxQuantFeatures(String version) {

        this.version = version;
        createOptions();
    }

    @Override
    protected void createOptions() {

        this.cmdLineOpts = new Options();

        cmdLineOpts.addOption(Option.builder("mqD").required(false).hasArg().longOpt("maxquantPsmDir").desc("MaxQuant psm root directory.").hasArg().build());
        cmdLineOpts.addOption(Option.builder("mqF").required(false).hasArg().longOpt("maxquantFeatures").desc("Comma separated list [feature1,feature2,feature3] of MaxQuant features to be added to NewAnce result file.").hasArg().build());
        cmdLineOpts.addOption(Option.builder("naf").required().hasArg().longOpt("newAnceResultFile").desc("Result file from NewAnce analysis (required)").build());
        cmdLineOpts.addOption(Option.builder("h").required(false).hasArg(false).longOpt("help").desc("Help option for command line help").build());
        cmdLineOpts.addOption(Option.builder("v").required(false).hasArg(false).longOpt("version").desc("Version of NewAnce software").build());
    }

    @Override
    protected void check(CommandLine line) throws ParseException {

        if (optionsSet) return;

        newAnceResultFile = new File(NewAnceParams.getFileValue("newAnceResultFile",getOptionString(line, "naf")));
        maxquantPsmDir = new File(NewAnceParams.getDirectoryValue("newAnceResultFile",getOptionString(line, "mqD")));
        maxQuantFeatures = NewAnceParams.getSetValue("newAnceResultFile",getOptionString(line, "mqF"));
    }

    @Override
    public int run() throws IOException {

        readNewAnceFile();
        readMaxQuantFiles();

        writeNewAnceFile();

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


            Set<String> invalidFeatures = new HashSet<>();
            for (String feature : maxQuantFeatures) {
                if (results.findColumn(feature)==0) {
                    invalidFeatures.add(feature);
                    System.out.println("WARNING: feature "+feature+" is not a valid msms.txt column name. This feature is ignored.");
                }
            }
            maxQuantFeatures.removeAll(invalidFeatures);

            while (results.next()) {

                int charge = results.getInt("Charge");
                int scanNumber = results.getInt("Scan number");
                String rawFile = results.getString("Raw file");

                String specID = rawFile + "." + scanNumber + "." + scanNumber + "." + charge;

                if (newAnceLines.containsKey(specID)) {

                    String values = "";
                    for (String feature : maxQuantFeatures) {
                        values += "\t"+results.getString(feature);
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

        newAnceLines = new HashMap<>();
        try {
            BufferedReader lineReader = new BufferedReader(new FileReader(newAnceResultFile ));

            String line = lineReader.readLine();
            newAnceLines.put("TITLE", line);

            while ((line = lineReader.readLine()) != null) {

                if (line.isEmpty()) continue;

                String[] fields = line.split("\t");

                newAnceLines.put(fields[0],line);

            }
            lineReader.close();
        } catch (FileNotFoundException e) {
        } catch (IOException e) {
        }
    }

    private void writeNewAnceFile() {
        try {

            String newNewAnceFile = newAnceResultFile.getAbsolutePath().replace("PSMs.txt","PSMs_extend.txt");
            BufferedWriter writer = new BufferedWriter(new FileWriter(newNewAnceFile));

            String header = newAnceLines.get("TITLE");
            for (String feature : maxQuantFeatures) {
                header += "\t"+feature;
            }
            writer.write(header+"\n");

            for  (String key : newAnceLines.keySet()) {

                if (!key.equals("TITLE")) {
                    String line = newAnceLines.get(key)+maxQuantLines.get(key);
                    writer.write(line+"\n");
                }
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

        AddMaxQuantFeatures addMaxQuantFeatures =  new AddMaxQuantFeatures("Version 1.0.0");
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
