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
import newance.mzjava.mol.modification.ModAttachment;
import newance.mzjava.mol.modification.Modification;
import newance.mzjava.mol.modification.ModificationList;
import newance.mzjava.ms.consensus.PeptideConsensusSpectrum;
import newance.mzjava.ms.io.mgf.MgfReader;
import newance.mzjava.ms.io.mgf.MgfWriter;
import newance.mzjava.ms.io.mgf.MsConvertTitleParser;
import newance.mzjava.ms.peaklist.PeakList;
import newance.mzjava.ms.peaklist.peakfilter.AbstractMergePeakFilter;
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

public class CreateSpectrumLibrary extends ExecutableOptions {

    private File mgfDir;
    private File newAnceResultFile;
    private final Map<String,Set<Integer>> spectrumFileMap;
    private final Map<String,MsnSpectrum> spectrumMap;
    private final Map<String,List<String>> peptideSpectrumMap;
    private final Map<String,Peptide> peptideMap;
    private final Map<String,File> mgfFileMap;

    public CreateSpectrumLibrary(String version) {

        this.version = version;
        this.spectrumFileMap = new HashMap<>();
        this.spectrumMap = new HashMap<>();
        this.peptideSpectrumMap = new HashMap<>();
        this.peptideMap = new HashMap<>();
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
        readSpectra();

        List<PeptideConsensusSpectrum> consensusSpectra = calculateConsensus();
        writeSpectra(consensusSpectra);
        writePDVTable(consensusSpectra);

        return 0;
    }

    private void writeSpectra(List<PeptideConsensusSpectrum> consensusSpectra) throws IOException  {

        String mgfFileName = newAnceResultFile.getAbsolutePath();

        if (mgfFileName.contains("PSMs"))
            mgfFileName = mgfFileName.replace("PSMs.txt","_SpecLib.mgf");
        else
            mgfFileName = mgfFileName.replace(".txt","_SpecLib.mgf");


        MgfWriter mgfWriter = new MgfWriter(new File(mgfFileName), PeakList.Precision.DOUBLE);
        for (PeptideConsensusSpectrum consensusSpectrum : consensusSpectra) {
            MsnSpectrum msnSpectrum = new MsnSpectrum();

            msnSpectrum.setComment(consensusSpectrum.getComment());
            msnSpectrum.setPrecursor(consensusSpectrum.getPrecursor());
            for (int i=0;i<consensusSpectrum.size();i++) {
                msnSpectrum.add(consensusSpectrum.getMz(i),consensusSpectrum.getIntensity(i));
            }
            mgfWriter.write(msnSpectrum);
        }

        mgfWriter.close();

    }

    private void writePDVTable(List<PeptideConsensusSpectrum> consensusSpectra) throws IOException {

        String pdvFileName = newAnceResultFile.getAbsolutePath();

        if (pdvFileName.contains("PSMs"))
            pdvFileName = pdvFileName.replace("PSMs.txt$","_cons_PDV.txt");
        else
            pdvFileName = pdvFileName.replace(".txt","_cons_PDV.txt");


        FileWriter fileWriter = new FileWriter(new File(pdvFileName));
        fileWriter.write("spectrum_title\tpeptide\tcharge\tpep_mass\tmz\tmodification\n");

        for (PeptideConsensusSpectrum consensusSpectrum: consensusSpectra) {

            Peptide peptide = consensusSpectrum.getPeptide();
            String sequence = peptide.toSymbolString();
            int charge = consensusSpectrum.getPrecursor().getCharge();
            String peptideKey = consensusSpectrum.getPeptide().toString() + "." + charge;
            double mass = consensusSpectrum.getPrecursor().getMass();
            double mz = consensusSpectrum.getPrecursor().getMz();

            fileWriter.write(peptideKey+".cons\t"+sequence+"\t"+charge+String.format("\t%.5f\t%.5f\t",mass,mz));

            String modifStr = "-";
            if (peptide.hasModifications()) {

                modifStr = "";
                int[] indexes = peptide.getModificationIndexes(ModAttachment.sideChainSet);
                for (int i : indexes) {
                    ModificationList modificationList = peptide.getModificationsAt(i, ModAttachment.sideChainSet);
                    for (int j=0;j<modificationList.size();j++) {
                        if (!modifStr.isEmpty()) modifStr += ";";
                        Modification modification = modificationList.get(j);
                        modifStr += modification.getLabel()+" of "+peptide.getSymbol(i)+"@"+(i+1)+
                                "["+modification.getMolecularMass()+"]";
                    }
                }

            }

            fileWriter.write(modifStr+"\n");
        }

        fileWriter.close();

    }

    private void readSpectra() {

        for (String mgfName : spectrumFileMap.keySet()) {
            File mgfFile = mgfFileMap.get(mgfName);

            try {

                System.out.println("Parsing and matching mgf file : " + mgfFile.getAbsolutePath());

                Set<Integer> scanNumbers = spectrumFileMap.get(mgfName);

                MgfReader reader = new MgfReader(mgfFile, new MsConvertTitleParser());

                MsnSpectrum spectrum = null;
                while (reader.hasNext()) {
                    try {
                        spectrum = reader.next();

                        int scanNr = spectrum.getScanNumbers().getFirst().getValue();
                        if (scanNumbers.contains(scanNr)) {
                            String specId = mgfName+"."+scanNr+"."+scanNr+"."+spectrum.getPrecursor().getCharge();
                            spectrum.setComment(specId);
                            spectrumMap.put(specId,spectrum);
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


            while (results.next()) {

                String spectrum = results.getString("Spectrum");
                System.out.println(spectrum);
                String msmsFileName = spectrum.substring(0,spectrum.indexOf('.'));
                int scanNr = Integer.parseInt(results.getString("ScanNr"));

                spectrumFileMap.putIfAbsent(msmsFileName,new HashSet<>());
                spectrumFileMap.get(msmsFileName).add(scanNr);

                String peptide = results.getString("Peptide");

                String peptideKey = peptide+"."+results.getInt("Charge");

                peptideSpectrumMap.putIfAbsent(peptideKey, new ArrayList<>());
                peptideSpectrumMap.get(peptideKey).add(spectrum);

                peptideMap.put(peptideKey,Peptide.parse(peptide));
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

    private List<PeptideConsensusSpectrum> calculateConsensus() {

        List<PeptideConsensusSpectrum> consensusSpectra = new ArrayList<>();

        for (String peptideKey : peptideSpectrumMap.keySet()) {
            List<MsnSpectrum> spectra = new ArrayList<>();
            for (String spectrumId: peptideSpectrumMap.get(peptideKey)) {
                spectra.add(spectrumMap.get(spectrumId));
            }

            List<PeptidePeakGenerator<PepFragAnnotation>> peakGeneratorList = new ArrayList<>();
            peakGeneratorList.add(new BackbonePeakGenerator(EnumSet.of(IonType.a, IonType.b, IonType.y), 10));
            Mass waterLoss = Composition.parseComposition("H-2O-1");
            Set<AminoAcid> aa = EnumSet.of(AminoAcid.S, AminoAcid.T, AminoAcid.D, AminoAcid.E);
            Set<IonType> ions = EnumSet.of(IonType.b, IonType.y);
            peakGeneratorList.add(new PeptideNeutralLossPeakGenerator(waterLoss, aa, ions, 5));

            PeptideConsensusSpectrum consensusSpectrum =
                    PeptideConsensusSpectrum.builder().
                            peptide(peptideMap.get(peptideKey)).
                            fragmenter(new PeptideFragmenter(peakGeneratorList, PeakList.Precision.DOUBLE)).
                            spectra(spectra).
                            setPeakFilterParams(0.2, 2).
                            fragMzTolerance(0.2).
                            intensityCombMethod(AbstractMergePeakFilter.IntensityMode.MEAN_ALL_INTENSITY).
                            build();

            consensusSpectrum.setComment(peptideKey+".cons");

            consensusSpectra.add(consensusSpectrum);
        }

        return consensusSpectra;
    }

    public ExecutableOptions init(String[] args) throws IOException {

        optionsSet = false;
        checkHelpOption(args, "-h");
        checkVersionOption(args, version, "-v");

        return this;
    }

    public static void main(String[] args) {

        CreateSpectrumLibrary createSpectrumLibrary =  new CreateSpectrumLibrary("Version 1.0.0");
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

