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

import newance.mzjava.mol.MassCalculator;
import newance.mzjava.ms.io.mgf.MgfReader;
import newance.mzjava.ms.io.mgf.MgfWriter;
import newance.mzjava.ms.io.mgf.MsConvertTitleParser;
import newance.mzjava.ms.peaklist.PeakList;
import newance.mzjava.ms.spectrum.MsnSpectrum;
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

public class CreatePDVExport extends ExecutableOptions {


    private File mgfDir;
    private File newAnceResultFile;
    private final Map<String,Set<Integer>> spectrumFileMap;
    private final Map<String,File> mgfFileMap;

    public CreatePDVExport(String version) {

        this.version = version;
        this.spectrumFileMap = new HashMap<>();
        this.mgfFileMap = new HashMap<>();
        createOptions();
    }

    @Override
    protected void createOptions() {

        this.cmdLineOpts = new Options();

        cmdLineOpts.addOption(Option.builder("mgf").required().hasArg().longOpt("mgfDir").desc("Directory containing mgf files (required)").build());
        cmdLineOpts.addOption(Option.builder("naf").required().hasArg().longOpt("newAnceResultFile").desc("Result file from NewAnce analysis (required)").build());
        cmdLineOpts.addOption(Option.builder("h").required(false).hasArg(false).longOpt("help").desc("Help option for command line help").build());
        cmdLineOpts.addOption(Option.builder("v").required(false).hasArg(false).longOpt("version").desc("Version of NewAnce software").build());
    }

    @Override
    protected void check(CommandLine line) throws ParseException {

        if (optionsSet) return;

        newAnceResultFile = new File(getOptionString(line,"naf"));
        mgfDir = new File(getOptionString(line,"mgf"));
    }

    @Override
    public int run() throws IOException {

        readNewAnceFile();
        findMgfFiles();

        String mgfFileName = newAnceResultFile.getAbsolutePath();

        if (mgfFileName.contains("PSMs"))
            mgfFileName = mgfFileName.replace("PSMs.txt","PDV.mgf");
        else
            mgfFileName = mgfFileName.replace(".txt",".mgf");

        MgfWriter writer = new MgfWriter(new File(mgfFileName), PeakList.Precision.DOUBLE);
        for (String mgfName : spectrumFileMap.keySet()) {
            List<MsnSpectrum> spectra = selectMatchedSpectra(mgfName, spectrumFileMap.get(mgfName));
            writeSpectra(writer,spectra);
        }
        writer.close();

        return 0;
    }

    private void writeSpectra(MgfWriter writer , List<MsnSpectrum> spectra) {

        spectra.forEach(sp -> {
            try {
                writer.write(sp);
            } catch (IOException e) {
                e.printStackTrace();
            }
        });
    }

    private void writePDVTable(FileWriter fileWriter, ResultSet results) throws IOException, SQLException {

        fileWriter.write(results.getString("Spectrum")+"\t"+results.getString("Sequence")+"\t");

        int charge = Integer.parseInt(results.getString("Charge"));
        double pm = Double.parseDouble(results.getString("PeptideMass"));
        String mdStr = results.getString("Massdiff");
        if (mdStr.startsWith("+")) mdStr = mdStr.substring(1);
        double dm = Double.parseDouble(mdStr);

        fileWriter.write(charge+"\t"+results.getString("PeptideMass")+"\t");

        double mz = (pm+dm+charge*MassCalculator.PROTON_MASS)/charge;
        fileWriter.write(String.format("%.5f",mz)+"\t");

        String modifStr = "-";
        if (!results.getString("ModifName").equals("NA")) {
            String[] modifNames = results.getString("ModifName").split(",");
            String[] modifAAs = results.getString("ModifAA").split(",");
            String[] modifPos = results.getString("ModifPosition").split(",");
            String[] modifMasses = results.getString("ModifMass").split(",");

            modifStr = "";
            for (int i=0;i<modifNames.length;i++) {
                if (!modifStr.isEmpty()) modifStr += ";";
                modifStr += modifNames[i]+" of "+modifAAs[i]+"@"+modifPos[i]+"["+modifMasses[i]+"]";
            }
        }

        fileWriter.write(modifStr+"\n");

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

    private List<MsnSpectrum> selectMatchedSpectra(String rawFileName, Set<Integer> scanNumbers) {

        PeakList.Precision precision = PeakList.Precision.DOUBLE;
        List<MsnSpectrum> spectra = new ArrayList<>();

        File mgfFile = new File(mgfDir.getAbsolutePath()+File.separator+rawFileName+".mgf");

        System.out.println("Parsing and matching mgf file : " + mgfFile.getAbsolutePath());

        try {

            MgfReader reader = new MgfReader(mgfFile, new MsConvertTitleParser());

            MsnSpectrum spectrum = null;
            while (reader.hasNext()) {
                try {
                    spectrum = reader.next();

                    int scanNr = spectrum.getScanNumbers().getFirst().getValue();
                    if (scanNumbers.contains(scanNr)) {
                        spectrum.setComment(rawFileName+"."+scanNr+"."+scanNr+"."+spectrum.getPrecursor().getCharge());
                        spectra.add(spectrum);
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


            String pdvFileName = newAnceResultFile.getAbsolutePath();
            pdvFileName = pdvFileName.replaceAll(".txt$","_PDV.txt");

            FileWriter fileWriter = new FileWriter(new File(pdvFileName));
            fileWriter.write("spectrum_title\tpeptide\tcharge\tpep_mass\tmz\tmodification\n");

            while (results.next()) {

                String spectrum = results.getString("Spectrum");
                System.out.println(spectrum);
                String msmsFileName = spectrum.substring(0,spectrum.indexOf('.'));
                int scanNr = Integer.parseInt(results.getString("ScanNr"));

                spectrumFileMap.putIfAbsent(msmsFileName,new HashSet<>());
                spectrumFileMap.get(msmsFileName).add(scanNr);

                writePDVTable(fileWriter, results);
             }

            // clean up
            fileWriter.close();
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

        CreatePDVExport createPDVExport =  new CreatePDVExport("Version 1.0.0");
        try {
            createPDVExport.init(args).parseOptions(args).run();
        } catch (MissingOptionException e) {
        } catch (ParseException e) {
            createPDVExport.printOptions(args, e.getMessage());
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}

