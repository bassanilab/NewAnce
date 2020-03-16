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

import newance.proteinmatch.PeptideUniProtSequenceMatch;
import newance.proteinmatch.UniProtDB;
import newance.util.ExecutableOptions;
import newance.util.NewAnceParams;
import org.apache.commons.cli.*;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

/**
 * Created by markusmueller on 25.03.17.
 */
public class AnnotatePeptides extends ExecutableOptions {

    private static final Logger LOGGER = Logger.getLogger(AnnotatePeptides.class.getName());

    private UniProtDB uniProtDB;
    private String columnName;
    private File fastaFile;
    private File peptideFile;
    private String header;
    private final List<String> peptides;
    private final Map<String,Map<String, List<PeptideUniProtSequenceMatch>>> matches;
    private final List<String> lines;


    public AnnotatePeptides(String fastaFile) {

        peptides = new ArrayList<>();
        matches = new HashMap<>();
        lines = new ArrayList<>();

        createOptions();
    }

    @Override
    protected void createOptions() {

        this.cmdLineOpts = new Options();

        cmdLineOpts.addOption(Option.builder("fa").required(true).hasArg().longOpt("fastaFile").desc("Fasta file with protein sequences.").hasArg().build());
        cmdLineOpts.addOption(Option.builder("p").required(true).hasArg().longOpt("peptideFile").desc("Peptide file in tab format. See -col argument.").hasArg().build());
        cmdLineOpts.addOption(Option.builder("col").required(false).hasArg().longOpt("peptideColumnName").desc("Name of coumn which contains peptide sequences to be searched. If not provided, first column without header is used.").hasArg().build());
        cmdLineOpts.addOption(Option.builder("h").required(false).hasArg(false).longOpt("help").desc("Help option for command line help").build());
        cmdLineOpts.addOption(Option.builder("v").required(false).hasArg(false).longOpt("version").desc("Version of NewAnce software").build());
    }

    @Override
    protected void check(CommandLine line) throws ParseException {

        if (optionsSet) return;

        fastaFile = new File(NewAnceParams.getFileValue("fastFileName",getOptionString(line, "fa")));
        peptideFile = new File(NewAnceParams.getFileValue("peptideFileName",getOptionString(line, "p")));
        columnName = getOptionString(line, "col");
    }

    @Override
    public int run() throws IOException {

        uniProtDB = new UniProtDB(fastaFile.getAbsolutePath());
        readPeptidesTxt(peptideFile);

        processPeptides();

        String resultFileName = peptideFile.getAbsolutePath();
        resultFileName = resultFileName.replaceFirst("\\..*$","");
        resultFileName += "_annot.txt";
        reportPeptides(new File(resultFileName));

        return 0;
    }

    private void processPeptides(){

        for (String peptide : peptides) {
            System.out.println("Processing peptide "+peptide+" ...");
            if (!matches.containsKey(peptide))
                matches.put(peptide, uniProtDB.findPeptide(peptide));
        }
    }

    private void reportPeptides(File outputFile){

        FileWriter writer = null;
        try {
            writer = new FileWriter(outputFile);
            if (!header.isEmpty())
                writer.write(header+"\n");

            for (int i=0;i<lines.size();i++) {

                String peptide = peptides.get(i);

                Map<String, List<PeptideUniProtSequenceMatch>> ppms = matches.get(peptide);
                writer.write(lines.get(i)+"\t"+ppms.keySet().toString()+"\n");
                writer.flush();
            }

            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void readPeptidesTxt(File peptideFile) {
        try {
            String line = "";
            BufferedReader br = new BufferedReader(new FileReader(peptideFile));

            line = br.readLine();
            String[] fields = line.split("\t");

            int index = -1;
            int i = 0;
            for (String field : fields) {
                if (field.equalsIgnoreCase(columnName)) {
                    index = i;
                    break;
                }
                i++;
            }

            if (index==-1) { // not header
                String peptide = fields[0];
                peptides.add(peptide);
                lines.add(line);
                header = "";

            } else {
                header = line+"\tFastaAnnot";
            }

            while ((line = br.readLine()) != null) {
                if (line.trim().isEmpty()) continue;

                fields = line.split("\t");
                String peptide = fields[index<0?0:index];
                peptides.add(peptide);
                lines.add(line);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {

        AnnotatePeptides annotatePeptides =  new AnnotatePeptides("Version 1.0.0");
        try {
            annotatePeptides.init(args).parseOptions(args).run();
        } catch (MissingOptionException e) {
        } catch (ParseException e) {
            annotatePeptides.printOptions(args, e.getMessage());
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
