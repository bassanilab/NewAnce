/**
 * Copyright (c) 2010, SIB. All rights reserved.
 *
 * SIB (Swiss Institute of Bioinformatics) - http://www.isb-sib.ch Host -
 * http://mzjava.expasy.org
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer. Redistributions in binary
 * form must reproduce the above copyright notice, this list of conditions and
 * the following disclaimer in the documentation and/or other materials provided
 * with the distribution. Neither the name of the SIB/GENEBIO nor the names of
 * its contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL SIB/GENEBIO BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package newance.psmconverter;

import com.google.common.base.Optional;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import org.apache.commons.io.FilenameUtils;
import org.expasy.mzjava.proteomics.mol.AminoAcid;
import org.expasy.mzjava.proteomics.mol.modification.ModAttachment;
import org.expasy.mzjava.proteomics.mol.modification.Modification;
import org.expasy.mzjava.proteomics.mol.modification.ModificationResolver;
import org.expasy.mzjava.proteomics.mol.modification.unimod.UnimodModificationResolver;

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

public class MaxQuantPsmReader2 {

    private final ModificationResolver modResolver;
    private final Map<String,String> peptidesProteinMap; // for decoys
    private final Map<String,String> peptidesMutationMap; // for decoys

    private final Pattern aaPattern = Pattern.compile("(.)\\((..)\\)|([A-Z])");

    public MaxQuantPsmReader2() {

        this.modResolver = makeDefaultModResolver();
        this.peptidesProteinMap = new HashMap<>();
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

        peptidesProteinMap.clear();

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

            while (results.next()) {
                String proteinStr = results.getString("Leading razor protein");
                if (proteinStr.startsWith("REV__")) {
                    proteinStr = proteinStr.replace("REV__","DECOY_");
                }

                String peptideSeq = results.getString("Sequence");
                peptidesProteinMap.put(peptideSeq,proteinStr);

                if (results.findColumn("Mutated")>0) {
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

        if (peptidesProteinMap.isEmpty()) {

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
            Connection conn = DriverManager.getConnection("jdbc:relique:csv:" + file.getParent() + "?separator=" + URLEncoder.encode(delimiter, "UTF-8"), props);

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

    public PeptideMatchDataWrapper makeModifiedPeptideMatch(String modifiedSequence) throws SQLException {

        Matcher matcher = aaPattern.matcher(modifiedSequence);

        List<AminoAcid> sequence = new ArrayList<>();
        ListMultimap<Object, Modification> modMatchMap = ArrayListMultimap.create();

        int size = 0;
        while (matcher.find()) {

            String modAA = matcher.group(1);
            String mod = matcher.group(2);
            String unModAA = matcher.group(3);

            if (modAA != null) {

                Optional<Modification> modOpt = modResolver.resolve(mod);
                if ("_".equals(modAA)) {

                    if (modOpt.isPresent()) {
                        modMatchMap.put(size == 0 ? ModAttachment.N_TERM : ModAttachment.C_TERM, modOpt.get());
                    }
                } else {

                    sequence.add(AminoAcid.valueOf(modAA));
                    if (modOpt.isPresent()) {
                        modMatchMap.put(size, modOpt.get());
                    }
                    size += 1;
                }
            } else {

                sequence.add(AminoAcid.valueOf(unModAA));
                size += 1;
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

        Set<String> proteins = new HashSet<>();
        for (String ac : getAccessionCode(results)) {

            boolean decoy =  ac.equals("DECOY");

            // other protein ids than uniprot
            if (decoy) {
                if (peptidesProteinMap.containsKey(peptideSeq))
                    proteins.add(peptidesProteinMap.get(peptideSeq));
                else
                    System.out.println("WARNING: peptide seq "+peptideSeq+" not found in peptides.txt.");
            } else {
                isDecoy = false;
                if (isVariant)
                    proteins.add(ac+" "+peptidesMutationMap.get(peptideSeq));
                else
                    proteins.add(ac);
            }
        }

        psm.setProteins(proteins);
        psm.setDecoy(isDecoy);
        psm.setVariant(isVariant);

        setValuesFirst(psm, results);

        return psm;
    }

   public PeptideMatchDataWrapper makeSecondModifiedPeptideMatch(ResultSet results) throws SQLException {

        String[] allSequences = results.getString("All modified sequences").split(";");

        if (allSequences.length<=1) return null;

        String peptide = allSequences[1];
        PeptideMatchDataWrapper psm = makeModifiedPeptideMatch(peptide);

        Set<String> proteins = new HashSet<>();
        proteins.add("unknown second hit");

        setValuesSecond(psm, results);

        return psm;
    }


    public List<PeptideMatchDataWrapper> makePeptideMatches(ResultSet results) throws SQLException {

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
