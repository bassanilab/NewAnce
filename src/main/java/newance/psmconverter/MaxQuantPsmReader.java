/*
Copyright (C) SIB - Swiss Institute of Bioinformatics, Lausanne, Switzerland
Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
*/

package newance.psmconverter;

import com.google.common.base.Optional;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import newance.mzjava.mol.AminoAcid;
import newance.mzjava.mol.modification.ModAttachment;
import newance.mzjava.mol.modification.Modification;
import newance.mzjava.mol.modification.ModificationResolver;
import newance.mzjava.mol.modification.unimod.UnimodModificationResolver;
import newance.proteinmatch.SequenceVariant;
import newance.proteinmatch.VariantProtDB;
import newance.psmcombiner.GroupedFDRCalculator;
import org.apache.commons.io.FilenameUtils;


import java.io.File;
import java.io.UnsupportedEncodingException;
import java.net.URLEncoder;
import java.sql.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author Markus Müller
 */

public class MaxQuantPsmReader {

    private final GroupedFDRCalculator groupedFDRCalculator;
    private final ModificationResolver modResolver;
    private final Map<String,String> peptidesLeadingProteinMap; // for decoys
    private final Map<String,String> peptidesMutationMap; // for decoys

    public MaxQuantPsmReader() {

        this.groupedFDRCalculator = null;
        this.modResolver = makeDefaultModResolver();
        this.peptidesLeadingProteinMap = new HashMap<>();
        this.peptidesMutationMap = new HashMap<>();
    }

    public MaxQuantPsmReader(GroupedFDRCalculator groupedFDRCalculator) {

        this.groupedFDRCalculator = groupedFDRCalculator;
        this.modResolver = makeDefaultModResolver();
        this.peptidesLeadingProteinMap = new HashMap<>();
        this.peptidesMutationMap = new HashMap<>();
    }

    private static UnimodModificationResolver makeDefaultModResolver() {

        UnimodModificationResolver modResolver = new UnimodModificationResolver();
        modResolver.putTranslate("de", "Deamidated");
        modResolver.putTranslate("ox", "Oxidation");
        modResolver.putTranslate("ac", "Acetyl");
        modResolver.putTranslate("ph", "Phospho");
        return modResolver;
    }

    private void extractPeptideProteinMap(File peptidesFile) {

        peptidesLeadingProteinMap.clear();
        peptidesMutationMap.clear();

        try {
            // load the driver into memory
            Class.forName("org.relique.jdbc.csv.CsvDriver");

            Properties props = new Properties();
            props.put("fileExtension", "." + FilenameUtils.getExtension(peptidesFile.getName()));
            // create a connection. The first command line parameter is assumed to
            //  be the directory in which the .csv files are held
            Connection conn = DriverManager.getConnection("jdbc:relique:csv:" + peptidesFile.getParent() + "?separator=" + URLEncoder.encode("\t", "UTF-8"), props);

            // create a Statement object to execute the query with
            Statement stmt = conn.createStatement();

            // Select the ID and NAME columns from sample.csv
            ResultSet results = stmt.executeQuery("SELECT * FROM " + FilenameUtils.removeExtension(peptidesFile.getName()));
            boolean hasMutations = false;
            try {
                hasMutations = results.findColumn("Mutated")>0;
            } catch (SQLException e) {
                System.out.println("No mutations were used in this MaxQuant search.");
            }

            while (results.next()) {
                String leadingProteinStr = results.getString("Leading razor protein");
                if (leadingProteinStr.startsWith("REV__")) {
                    leadingProteinStr = leadingProteinStr.replace("REV__","DECOY_");
                }

                String peptideSeq = results.getString("Sequence");
                peptidesLeadingProteinMap.put(peptideSeq,leadingProteinStr);

                if (hasMutations) {
                    String isVariant = results.getString("Mutated");

                    if (isVariant.equals("Yes")) {
                        peptidesMutationMap.put(peptideSeq, results.getString("Mutation names"));
                    }
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

    public void parse(File file, PeptideSpectrumMatchList callback) {

        if (peptidesLeadingProteinMap.isEmpty()) {

            File peptidesFile = new File(file.getParent()+File.separator+"peptides.txt");
            if (peptidesFile.exists()) {
                extractPeptideProteinMap(peptidesFile);
            }
            else {
                System.out.println("WARNING: peptides file "+peptidesFile.getAbsolutePath()+" not found !");
            }
        }

        parse(file, "\t", callback);
    }

    public void parse(File file, String delimiter, PeptideSpectrumMatchList peptideSpectrumMatchList) {

        try {
            // load the driver into memory
            Class.forName("org.relique.jdbc.csv.CsvDriver");

            Properties props = new Properties();
            props.put("fileExtension", "." + FilenameUtils.getExtension(file.getName()));
            // create a connection. The first command line parameter is assumed to
            //  be the directory in which the .csv files are held
            Connection conn = DriverManager.getConnection("jdbc:relique:csv:" + file.getParent() + "?separator="
                    + URLEncoder.encode(delimiter, "UTF-8"), props);

            // create a Statement object to execute the query with
            Statement stmt = conn.createStatement();

            // Select the ID and NAME columns from sample.csv
            ResultSet results = stmt.executeQuery("SELECT * FROM " + FilenameUtils.removeExtension(file.getName()));

            while (results.next()) {

                SpectrumInfo spectrumInfo = getSpectrumInfo(results);
                if (!peptideSpectrumMatchList.isValidSpectrum(spectrumInfo)) continue;

                List<PeptideMatchDataWrapper> psms = makePeptideMatches(results);
                for (PeptideMatchDataWrapper psm : psms) {

                    if (peptideSpectrumMatchList.isValidProtein(psm.getProteins())) peptideSpectrumMatchList.resultRead(spectrumInfo, psm);
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


    protected int parseModification(int idx, int pos, char[] chars, ListMultimap<Object, Modification> modMatchMap) {

        int stackCnt = 1;
        int j=idx+1;
        int modifNameEndIdx = j;
        while (stackCnt>0 && j<chars.length) {
            if (chars[j]=='(') stackCnt++;
            if (chars[j]==')') stackCnt--;

            if (chars[j]==' ' || chars[j]=='(' || stackCnt==0) {
                if (modifNameEndIdx==idx+1) modifNameEndIdx = j;
            }

            j++;
        }

        Optional<Modification> modOpt = modResolver.resolve(new String(chars,idx+1,modifNameEndIdx-idx-1));
        if ('_'==chars[idx-1]) {

            if (modOpt.isPresent()) {
                modMatchMap.put(idx == 1 ? ModAttachment.N_TERM : ModAttachment.C_TERM, modOpt.get());
            }
        } else {

            if (modOpt.isPresent()) {
                modMatchMap.put(pos, modOpt.get());
            }
        }

        return j;
    }

    protected PeptideMatchDataWrapper makeModifiedPeptideMatch(String modifiedSequence) throws SQLException {

        List<AminoAcid> sequence = new ArrayList<>();
        ListMultimap<Object, Modification> modMatchMap = ArrayListMultimap.create();

        char[] chars = modifiedSequence.toCharArray();
        int pos = 0;
        int i = 0;
        while (i<chars.length) {
            if (chars[i]=='_') {
                i++;
                continue;
            }

            if (chars[i]=='(') {
                i = parseModification(i, pos-1, chars, modMatchMap);
            } else {
                sequence.add(AminoAcid.valueOf(chars[i]));
                pos += 1;
                i++;
            }
        }

        PeptideMatchDataWrapper peptideMatch = new PeptideMatchDataWrapper(sequence);
        for (Object key : modMatchMap.keySet()) {

            for (Modification mod : modMatchMap.get(key)) {

                if (key instanceof Integer) {

                    int position = (Integer) key;
                    peptideMatch.addModificationMatch(position, mod);
                } else if (key instanceof ModAttachment) {

                    ModAttachment modAttachment = (ModAttachment) key;
                    peptideMatch.addModificationMatch(modAttachment, mod);
                }
            }
        }

        return peptideMatch;
    }

    public PeptideMatchDataWrapper makeFirstModifiedPeptideMatch(ResultSet results) throws SQLException {

        String peptide = results.getString("Modified sequence");
        PeptideMatchDataWrapper psm = makeModifiedPeptideMatch(peptide);
        String peptideSeq = results.getString("Sequence");
        boolean isDecoy = true;
        boolean isVariant = peptidesMutationMap.containsKey(peptideSeq);

        psm.setVariant(isVariant);

        List<String> proteins = new ArrayList<>();
        for (String ac : getAccessionCode(results)) {

            boolean decoy =  ac.equals("DECOY");

            // other protein ids than uniprot
            String leadingProtein = peptidesLeadingProteinMap.get(peptideSeq);
            psm.setLeadingProtein(leadingProtein);

            if (decoy) {
                if (peptidesLeadingProteinMap.containsKey(peptideSeq))
                    proteins.add(peptidesLeadingProteinMap.get(peptideSeq));
                else
                    System.out.println("WARNING: peptide seq "+peptideSeq+" not found in peptides.txt.");
            } else {
                isDecoy = false;
                proteins.add(ac);
                if (isVariant) psm.setVariant(true);
            }
        }

        psm.setProteins(proteins);
        psm.setDecoy(isDecoy);

        setValuesFirst(psm, results);

        return psm;
    }

    protected PeptideMatchDataWrapper makeSecondModifiedPeptideMatch(ResultSet results) throws SQLException {

        String[] allSequences = results.getString("All modified sequences").split(";");

        if (allSequences.length<=1) return null;

        String peptide = allSequences[1];
        PeptideMatchDataWrapper psm = makeModifiedPeptideMatch(peptide);

        Set<String> proteins = new HashSet<>();
        proteins.add("unknown second hit");

        setValuesSecond(psm, results);

        return psm;
    }


    protected List<PeptideMatchDataWrapper> makePeptideMatches(ResultSet results) throws SQLException {

        List<PeptideMatchDataWrapper> psms = new ArrayList<>();
        psms.add(makeFirstModifiedPeptideMatch(results));
//        PeptideMatchDataWrapper psm = makeSecondModifiedPeptideMatch(results);
//        if (psm!=null) psms.add(psm);

        return psms;
    }

    protected SpectrumInfo getSpectrumInfo(ResultSet results) throws SQLException {

        int charge = results.getInt("Charge");
        int scanNumber = results.getInt("Scan number");
        int scanIndex = results.getInt("Scan index");
        String rawFile = results.getString("Raw file");

        SpectrumInfo spectrumInfo = new SpectrumInfo(rawFile + "." + scanNumber + "." + scanNumber + "." + charge);
        spectrumInfo.setCharge(charge);
        spectrumInfo.setIndex(scanIndex);
        spectrumInfo.setScanNumber(scanNumber);
        spectrumInfo.setPrecursorNeutralMass(results.getDouble("Mass"));
        spectrumInfo.setRetentionTime(results.getDouble("Retention time"));

        return spectrumInfo;
    }

    protected Collection<String> getAccessionCode(ResultSet results) throws SQLException {

        String proteins = results.getString("Proteins");
        String reverse = results.getString("Reverse");
        if (proteins != null) {

            if (reverse.equals("+")) {
                return Lists.newArrayList("DECOY");
            }

            String[] accessions = proteins.split(";");
            return Lists.newArrayList(accessions);
        } else {

            return Collections.emptyList();
        }
    }

    protected void setValuesFirst(PeptideMatchDataWrapper peptideMatch, ResultSet results) throws SQLException {

        peptideMatch.setNumMissedCleavages(results.getInt("Missed cleavages"));
        peptideMatch.addScore("Number of Matches", results.getDouble("Number of Matches"));
        peptideMatch.addScore("Score", results.getDouble("Score"));
        peptideMatch.addScore("Delta score", results.getDouble("Delta score"));
        peptideMatch.addScore("Mass Error [ppm]", results.getDouble("Mass Error [ppm]"));
        peptideMatch.addScore("Intensity coverage", results.getDouble("Intensity coverage"));
        peptideMatch.addScore("Localization prob", results.getDouble("Localization prob"));
        peptideMatch.addScore("PEP", results.getDouble("PEP"));
        peptideMatch.setRank(1);
    }

    protected void setValuesSecond(PeptideMatchDataWrapper peptideMatch, ResultSet results) throws SQLException {

        String[] allScores = results.getString("All scores").split(";");

        double score = Double.parseDouble(allScores[1]);

        double deltaScore = score;
        if (allScores.length==3)
            deltaScore = score - Double.parseDouble(allScores[2]);

        peptideMatch.addScore("Score", score);
        peptideMatch.addScore("Delta score", deltaScore);
        peptideMatch.setRank(2);
    }
}
