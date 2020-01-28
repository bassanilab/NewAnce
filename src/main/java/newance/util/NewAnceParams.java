/**
 * Copyright (C) 2019, SIB/LICR. All rights reserved
 *
 * SIB, Swiss Institute of Bioinformatics
 * Ludwig Institute for Cancer Research (LICR)
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer. Redistributions in binary
 * form must reproduce the above copyright notice, this list of conditions and
 * the following disclaimer in the documentation and/or other materials provided
 * with the distribution. Neither the name of the SIB/LICR nor the names of
 * its contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL SIB/LICR BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

package newance.util;


import org.expasy.mzjava.proteomics.mol.modification.Modification;

import java.io.*;
import java.nio.file.InvalidPathException;
import java.util.*;
import java.util.regex.Pattern;

/**
 * @author Markus Müller
 */

public class NewAnceParams implements Serializable {

    private static NewAnceParams instance = null;

    // When reading PSM results files the modifications are mapped to the modifications defined here. Modifications are mapped by means of there mass only
    // (this will change in future versions, where the mapping will be done by mass and amino acid)

    private String modificationStr = "Cysteinyl:C3H5NO2S,Oxidation:O,Carbamidomethyl:C2H3NO,Phospho:HO3P,Deamidated:H-1N-1O," +
            "Ammonia-loss:H-3N-1,Acetyl:C2H2O,Methyl:CH2,Dimethyl:C2H4,Trimethyl:C3H6,Amidated:H-1,Glu->pyro-Glu:H-2O-1," +
            "Propionyl:C3H4O,Pyro-carbamidomethyl:C2O";

    private Set<Modification> modifications = null;

    // Mass tolerance for the modification mapping. This value is independent of instrument precision, but refers to the precision of the modification mass values
    // indicated in the PSM result files.
    private double modifMatchMassTol =  0.015;

    // Decoy prefix in comet
    private String cometDecoyProtPrefix = "DECOY_";

    // Remove excluded proteins from list of psm proteins for every PSM. remove PSM if there are only excluded proteins.
    private Pattern excludedProtPattern = null;

    // Predicate to filter PeptideSpectrumMatch objects before writing them to hadoop file
    private int minCharge = 1;
    private int maxCharge = 5;
    private int minPeptideLength = 8;
    private int maxPeptideLength = 25;
    private int maxRank = 1;

    private int nrThreads = 1;

    private double fdrCometThreshold = 0.03;

    private String protCodingGroup = "";
    private String noncanonicalGroup = "";
    private Pattern spectrumRegExp = null;
    private Pattern codingProtRegExp = null;
    private Set<String> forcedNoncanonicalProts = new HashSet<>();

    private int minNrPsmsPerHisto = 100000;

    private String outputPrefix = "";

    private int nrBins3D = 40;
    private double minXCorr = 0.0;
    private double maxXCorr = 5.0;
    private int nrXCorrBins = nrBins3D;
    private double minDeltaCn = 0.0;
    private double maxDeltaCn = 1.0;
    private int nrDeltaCnBins = nrBins3D;
    private double minSpScore = 0.0;
    private double maxSpScore = 2500.0;
    private int nrSpScoreBins = nrBins3D;

    private int smoothDegree = 1;

    private String fdrControlMethod = "global";

    private String version = "version 1.4.0";

    private String cometPsmDir = "";
    private String maxquantPsmDir = "";
    private Pattern cometPsmRegExp = null;
    private Pattern maxquantPsmRegExp = null;
    private boolean includeMaxQuant = true;
    private String outputDir = "";
    private boolean reportHistos = false;
    private String readHistos = "";
    private String searchFastaFile = "";
    private String uniprotFastaFile= "";
    private boolean doPeptideProteinGrouping = false;
    private String writeParamsFile = "";
    private String readParamsFile = "";

    private String maxQuantMainScore = "Score";
    private String cometMainScore = "xcorr";
    private double maxQuantMainScoreMinValue = 10f;
    private double cometMainScoreMinValue = 1f;

    private final Map<String,String> variableValueMap;

    public static NewAnceParams getInstance() {
        if (instance==null) {
            instance = new NewAnceParams();
        }

        return instance;
    }

    private NewAnceParams() {
        variableValueMap = new HashMap<>();
        maxquantPsmRegExp = Pattern.compile("^msms.txt$");

        modifications = new HashSet<>();
        for (String modif : modificationStr.split(",")){
            modifications.add(Modification.parseModification(modif));
        }
    }

    public String toString() {
        String res = "";

        res +=  "modifications="+Arrays.asList(modifications).toString()+"\n";
        res +=  "modifMatchMassTol="+modifMatchMassTol+"\n";
        res +=  "cometDecoyProtPrefix="+ cometDecoyProtPrefix +"\n";
        res +=  "excludedProtPattern="+excludedProtPattern+"\n";
        res +=  "minCharge="+minCharge+"\n";
        res +=  "maxCharge="+maxCharge+"\n";
        res +=  "minPeptideLength="+minPeptideLength+"\n";
        res +=  "maxPeptideLength="+maxPeptideLength+"\n";
        res +=  "maxRank="+maxRank+"\n";
        res +=  "smoothDegree="+smoothDegree+"\n";
        res +=  "nrThreads="+nrThreads+"\n";
        res +=  "fdrCometThreshold="+fdrCometThreshold+"\n";
        res +=  "protCodingGroup="+ protCodingGroup +"\n";
        res +=  "noncanonicalGroup="+ noncanonicalGroup +"\n";
        res +=  "spectrumRegExp="+spectrumRegExp+"\n";
        res +=  "codingProtRegExp="+ codingProtRegExp +"\n";
        res +=  "forcedNoncanonicalProts="+ forcedNoncanonicalProts +"\n";
        res +=  "minNrPsmsPerHisto="+minNrPsmsPerHisto+"\n";
        res +=  "outputPrefix="+outputPrefix+"\n";
        res +=  "minXCorr="+minXCorr+"\n";
        res +=  "maxXCorr="+maxXCorr+"\n";
        res +=  "nrXCorrBins="+nrXCorrBins+"\n";
        res +=  "minDeltaCn="+minDeltaCn+"\n";
        res +=  "maxDeltaCn="+maxDeltaCn+"\n";
        res +=  "nrDeltaCnBins="+nrDeltaCnBins+"\n";
        res +=  "minSpScore="+minSpScore+"\n";
        res +=  "maxSpScore="+maxSpScore+"\n";
        res +=  "nrSpScoreBins="+nrSpScoreBins+"\n";
        res +=  "fdrControlMethod="+fdrControlMethod+"\n";
        res +=  "version="+version+"\n";
        res +=  "cometPsmDir="+cometPsmDir+"\n";
        res +=  "maxquantPsmDir="+maxquantPsmDir+"\n";
        res +=  "cometPsmRegExp="+cometPsmRegExp+"\n";
        res +=  "maxquantPsmRegExp="+maxquantPsmRegExp+"\n";
        res +=  "includeMaxQuant="+includeMaxQuant+"\n";
        res +=  "outputDir="+outputDir+"\n";
        res +=  "reportHistos="+reportHistos+"\n";
        res +=  "readHistos="+readHistos+"\n";
        res +=  "searchFastaFile="+searchFastaFile+"\n";
        res +=  "uniprotFastaFile="+uniprotFastaFile+"\n";
        res +=  "doPeptideProteinGrouping="+doPeptideProteinGrouping+"\n";
        res +=  "writeParamsFile="+writeParamsFile+"\n";
        res +=  "readParamsFile="+readParamsFile+"\n";

        return res;
    }

    public void read(String paramsFileName) {
        try {
            BufferedReader reader = new BufferedReader(new FileReader(paramsFileName));
            String line;

            while ((line = reader.readLine()) != null) {

                line = line.trim();
                if (line.isEmpty()) continue;
                if (line.startsWith("#")) continue;
                if (!line.startsWith("#") && !line.contains("=")) continue;

                String[] fields = line.split("=");
                if (fields.length==2) variableValueMap.put(fields[0].trim(),fields[1].trim());
                else variableValueMap.put(fields[0].trim(),"");
            }

            reader.close();

        } catch (IOException e) {
            System.out.println("Cannot write NewAnceParams to file: "+paramsFileName);
        }

        setReadVariables();
    }

    private void setReadVariables() {

        if (variableValueMap.containsKey("modifications")) {
            String buf = "";
            try {
                for (String modif : getSetValue("modifications", variableValueMap.get("modifications"))) {
                    buf = modif;
                    modifications.add(Modification.parseModification(modif));
                }
            } catch(IllegalArgumentException e) {
                new RuntimeException("Invalid modification string "+buf+" provided. Abort");
            }
        }

        if (variableValueMap.containsKey("modifMatchMassTol")) {
            modifMatchMassTol = getDoubleValue("modifMatchMassTol", variableValueMap.get("modifMatchMassTol"));
        }

        if (variableValueMap.containsKey("cometDecoyProtPrefix")) {
            cometDecoyProtPrefix = getStringValue("cometDecoyProtPrefix",variableValueMap.get("cometDecoyProtPrefix"));
        }

        if (variableValueMap.containsKey("excludedProtPattern")) {
            excludedProtPattern = getPatternValue("excludedProtPattern",variableValueMap.get("excludedProtPattern"));
        }

        if (variableValueMap.containsKey("minCharge")) {
            minCharge = getIntegerValue("minCharge",variableValueMap.get("minCharge"));
        }

        if (variableValueMap.containsKey("maxCharge")) {
            maxCharge = getIntegerValue("maxCharge",variableValueMap.get("maxCharge"));
        }

        if (variableValueMap.containsKey("minPeptideLength")) {
            minPeptideLength = getIntegerValue("minPeptideLength",variableValueMap.get("minPeptideLength"));
        }

        if (variableValueMap.containsKey("maxPeptideLength")) {
            maxPeptideLength = getIntegerValue("maxPeptideLength",variableValueMap.get("maxPeptideLength"));
        }

        if (variableValueMap.containsKey("maxRank")) {
            maxRank = getIntegerValue("maxRank",variableValueMap.get("maxRank"));
        }

        if (variableValueMap.containsKey("nrThreads")) {
            nrThreads = getIntegerValue("nrThreads",variableValueMap.get("nrThreads"));
        } else {
            int nrProc = Runtime.getRuntime().availableProcessors();
            nrThreads = (nrProc>2)?nrProc-2:1;
        }

        if (variableValueMap.containsKey("smoothDegree")) {
            smoothDegree = getIntegerValue("smoothDegree",variableValueMap.get("smoothDegree"));
        }

        if (variableValueMap.containsKey("fdrCometThreshold")) {
            fdrCometThreshold = getDoubleValue("fdrCometThreshold",variableValueMap.get("fdrCometThreshold"));
        }

        if (variableValueMap.containsKey("codingProtRegExp")) {
            codingProtRegExp = getPatternValue("codingProtRegExp",variableValueMap.get("codingProtRegExp"));
        }

        if (codingProtRegExp !=null) {
            if (variableValueMap.containsKey("protCodingGroup")) {
                protCodingGroup = getStringValue("protCodingGroup",variableValueMap.get("protCodingGroup"));
            } else {
                throw new RuntimeException("No name for protein coding group provided. Abort");
            }

            if (variableValueMap.containsKey("noncanonicalGroup")) {
                noncanonicalGroup = getStringValue("noncanonicalGroup",variableValueMap.get("noncanonicalGroup"));
            } else {
                throw new RuntimeException("No name for non-canonical group provided. Abort");
            }
        }

        if (variableValueMap.containsKey("spectrumRegExp")) {
            spectrumRegExp = getPatternValue("spectrumRegExp",variableValueMap.get("spectrumRegExp"));
        }

        if (variableValueMap.containsKey("forcedNoncanonicalProts")) {
            forcedNoncanonicalProts = getSetValue("forcedNoncanonicalProts", variableValueMap.get("forcedNoncanonicalProts"));
        }

        if (variableValueMap.containsKey("minNrPsmsPerHisto")) {
            minNrPsmsPerHisto = getIntegerValue("minNrPsmsPerHisto",variableValueMap.get("minNrPsmsPerHisto"));
        }

        if (variableValueMap.containsKey("outputPrefix")) {
            outputPrefix = variableValueMap.get("outputPrefix"); // can be empty
        }

        if (variableValueMap.containsKey("minXCorr")) {
            minXCorr = getDoubleValue("minXCorr",variableValueMap.get("minXCorr"));
        }

        if (variableValueMap.containsKey("maxXCorr")) {
            maxXCorr = getDoubleValue("maxXCorr",variableValueMap.get("maxXCorr"));
        }

        if (variableValueMap.containsKey("nrXCorrBins")) {
            nrXCorrBins = getIntegerValue("nrXCorrBins",variableValueMap.get("nrXCorrBins"));
        }

        if (variableValueMap.containsKey("minDeltaCn")) {
            minDeltaCn = getDoubleValue("minDeltaCn",variableValueMap.get("minDeltaCn"));
        }

        if (variableValueMap.containsKey("maxDeltaCn")) {
            maxDeltaCn = getDoubleValue("maxDeltaCn",variableValueMap.get("maxDeltaCn"));
        }

        if (variableValueMap.containsKey("nrDeltaCnBins")) {
            nrDeltaCnBins = getIntegerValue("nrDeltaCnBins",variableValueMap.get("nrDeltaCnBins"));
        }

        if (variableValueMap.containsKey("minSpScore")) {
            minSpScore = getDoubleValue("minSpScore",variableValueMap.get("minSpScore"));
        }

        if (variableValueMap.containsKey("maxSpScore")) {
            maxSpScore = getDoubleValue("maxSpScore",variableValueMap.get("maxSpScore"));
        }

        if (variableValueMap.containsKey("nrSpScoreBins")) {
            nrSpScoreBins = getIntegerValue("nrSpScoreBins",variableValueMap.get("nrSpScoreBins"));
        }

        if (variableValueMap.containsKey("fdrControlMethod")) {
            fdrControlMethod = getStringValue("fdrControlMethod",variableValueMap.get("fdrControlMethod").trim().toLowerCase(),
            new HashSet<>(Arrays.asList(new String[]{"global","groupwise"})));
        }

        if (variableValueMap.containsKey("cometPsmDir")) {
            cometPsmDir = getDirectoryValue("cometPsmDir",variableValueMap.get("cometPsmDir"));
        } else {
            throw new RuntimeException("No valid value for Comet directory provided. Abort");
        }

        if (variableValueMap.containsKey("cometPsmRegExp")) {
            cometPsmRegExp = getPatternValue("cometPsmRegExp",variableValueMap.get("cometPsmRegExp"));
        } else {
            throw new RuntimeException("No valid value for Comet pep.xml file regexp provided. Abort");
        }

        if (variableValueMap.containsKey("includeMaxQuant")) {
            includeMaxQuant = getBooleanValue("includeMaxQuant",variableValueMap.get("includeMaxQuant"));
        }

        if (includeMaxQuant) {

            if (variableValueMap.containsKey("maxquantPsmDir")) {
                maxquantPsmDir = getDirectoryValue("maxquantPsmDir", variableValueMap.get("maxquantPsmDir"));
            } else {
                throw new RuntimeException("No valid value for MaxQuant directory provided. Abort");
            }

            if (variableValueMap.containsKey("maxquantPsmRegExp")) {
                maxquantPsmRegExp = getPatternValue("maxquantPsmRegExp", variableValueMap.get("maxquantPsmRegExp"));
            }
        }

        if (variableValueMap.containsKey("outputDir")) {
            outputDir = getNewDirectoryValue("outputDir",variableValueMap.get("outputDir"));
        } else {
            throw new RuntimeException("No valid value for output directory provided. Abort");
        }

        if (variableValueMap.containsKey("reportHistos")) {
            reportHistos = getBooleanValue("reportHistos",variableValueMap.get("reportHistos"));
        }

        if (variableValueMap.containsKey("readHistos")) {
            readHistos = getDirectoryValue("readHistos",variableValueMap.get("readHistos"));
        }

        if (variableValueMap.containsKey("searchFastaFile")) {
            searchFastaFile = getFileValue("searchFastaFile",variableValueMap.get("searchFastaFile"));
        }

        if (variableValueMap.containsKey("uniprotFastaFile")) {
            uniprotFastaFile = getFileValue("uniprotFastaFile",variableValueMap.get("uniprotFastaFile"));
        } else {
            throw new RuntimeException("No valid value for Uniprot fasta file provided. Abort");
        }

        if (variableValueMap.containsKey("doPeptideProteinGrouping")) {
            doPeptideProteinGrouping = getBooleanValue("doPeptideProteinGrouping",variableValueMap.get("doPeptideProteinGrouping"));
        }

        if (variableValueMap.containsKey("writeParamsFile")) {
            writeParamsFile = getNewFileValue("writeParamsFile",variableValueMap.get("writeParamsFile"));
        }

        if (variableValueMap.containsKey("readParamsFile")) {
            readParamsFile = getNewFileValue("readParamsFile",variableValueMap.get("readParamsFile"));
        }

        if (variableValueMap.containsKey("maxQuantMainScoreMinValue")) {
            maxSpScore = getDoubleValue("maxQuantMainScoreMinValue",variableValueMap.get("maxQuantMainScoreMinValue"));
        }

        if (variableValueMap.containsKey("cometMainScoreMinValue")) {
            maxSpScore = getDoubleValue("cometMainScoreMinValue",variableValueMap.get("cometMainScoreMinValue"));
        }

        if (variableValueMap.containsKey("maxQuantMainScore")) {
            outputPrefix = variableValueMap.get("maxQuantMainScore");
        }

        if (variableValueMap.containsKey("cometMainScore")) {
            outputPrefix = variableValueMap.get("cometMainScore");
        }

        checkVariableValues();

        if (!writeParamsFile.isEmpty()) {
            write(writeParamsFile);
        }
    }

    private boolean checkVariableValues() {

        if (minCharge<0 || minCharge>maxCharge) throw new RuntimeException("Error in minCharge value: "+minCharge+". Abort.");
        if (minPeptideLength<0 || minPeptideLength>maxPeptideLength) throw new RuntimeException("Error in maxPeptideLength value: "+minPeptideLength+". Abort.");
        if (modifMatchMassTol<0.0) throw new RuntimeException("Error in modifMatchMassTol value: "+modifMatchMassTol+". Abort.");
        if (maxRank<=0) throw new RuntimeException("Error in maxRank value: "+maxRank+". Abort.");
        if (nrThreads<=0) throw new RuntimeException("Error in nrThreads value: "+nrThreads+". Abort.");
        if (fdrCometThreshold<0.0 || fdrCometThreshold>1.0) throw new RuntimeException("Error in fdrCometThreshold value: "+fdrCometThreshold+". Abort.");
        if (maxRank<=0) throw new RuntimeException("Error in maxRank value: "+maxRank+". Abort.");

        if (minNrPsmsPerHisto<=0) throw new RuntimeException("Error in minNrPsmsPerHisto value: "+minNrPsmsPerHisto+". Abort.");
        if (minXCorr>=maxXCorr) throw new RuntimeException("Error in minXCorr value: "+minXCorr+". Abort.");
        if (minDeltaCn>=maxDeltaCn) throw new RuntimeException("Error in minDeltaCn value: "+minDeltaCn+". Abort.");
        if (minSpScore>=maxSpScore) throw new RuntimeException("Error in minSpScore value: "+minSpScore+". Abort.");
        if (nrXCorrBins<=0) throw new RuntimeException("Error in nrXCorrBins value: "+nrXCorrBins+". Abort.");
        if (nrDeltaCnBins<=0) throw new RuntimeException("Error in nrDeltaCnBins value: "+nrDeltaCnBins+". Abort.");
        if (nrSpScoreBins<=0) throw new RuntimeException("Error in nrSpScoreBins value: "+nrSpScoreBins+". Abort.");
        if (maxRank<=0) throw new RuntimeException("Error in maxRank value: "+maxRank+". Abort.");

        if (doPeptideProteinGrouping && !(new File(searchFastaFile)).exists()) throw new RuntimeException("Error in doPeptideProteinGrouping value: "+doPeptideProteinGrouping+". Abort.");

        return true;
    }

    public void add(String variable, String value) {

        if (value.trim().isEmpty()) return;

        variableValueMap.put(variable, value);
    }

    public void finalize() {
        setReadVariables();
    }

    public void write(String paramsFileName) {
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(paramsFileName));
            writer.write(toString());
            writer.close();
        } catch (IOException e) {
            System.out.println("Cannot write NewAnceParams to file: "+paramsFileName);
        }
    }

    private Set<Modification> getModifications(String variable, String value) {

        Set<String> modifStrs = getSetValue(variable, value);
        Set<Modification> modifications = new HashSet<>();

        for (String modifStr : modifStrs) {
            try {
                modifications.add(Modification.parseModification(modifStr.trim()));
            } catch(IllegalArgumentException e) {

                throw new RuntimeException("Invalid value "+value+" for variable "+variable+". Value is not a modification.");
            }
        }

        return modifications;
    }

    private String getStringValue(String variable, String value, Set<String> allowedValues) {

        if (allowedValues.contains(value)) {
            return value;
        } else {
            throw new RuntimeException("Invalid value "+value+" for variable "+variable+". Value has to be one of "+allowedValues.toString());
        }
    }

    private String getStringValue(String variable, String value) {

        if (!value.isEmpty()) {
            return value;
        } else {
            throw new RuntimeException("Invalid value "+value+" for variable "+variable+". Value is empty");
        }
    }

    private Pattern getPatternValue(String variable, String value) {

        try {
            return Pattern.compile(value);
        } catch (NumberFormatException e) {
            throw new RuntimeException("Invalid value "+value+" for variable "+variable+". Not a valid regular expression.");
        }
    }

    private boolean getBooleanValue(String variable, String value) {

        if (value.toLowerCase().equals("true")) {
            return true;
        } else if (value.toLowerCase().equals("false")) {
            return false;
        } else {
            throw new RuntimeException("Invalid value "+value+" for variable "+variable+". Value has to be either true or false.");
        }
    }

    private String getDirectoryValue(String variable, String value) {

        File dir = new File(value);
        if (!dir.exists() && !dir.isDirectory()) {
            throw new InvalidPathException("Invalid value "+value+" for variable "+variable+"."," Directory does not exist or is not a directory.");
        }

        return value;
    }

    private String getNewDirectoryValue(String variable, String value) {

        File dir = new File(value);

        if (!dir.exists()) {
            dir.mkdirs();
        }

        if (!dir.isDirectory()) {
            throw new InvalidPathException("Invalid value "+value+" for variable "+variable+"."," Is not a directory.");
        }

        return value;
    }

    private String getFileValue(String variable, String value) {

        File file = new File(value);
        if (!file.exists()) {
            throw new InvalidPathException("Invalid value "+value+" for variable "+variable+"."," File does not exist.");
        }

        return value;
    }

    private String getNewFileValue(String variable, String value) {

        File file = new File(value);
        try {
            if (file.exists()) file.delete();
            if (file.createNewFile()) return value;
        } catch (IOException e) {
            throw new InvalidPathException("Invalid value "+value+" for variable "+variable+"."," File cannot be created.");
        }

        return value;
    }

    private Double getDoubleValue(String variable, String value) {

        try {
            return Double.valueOf(value);
        } catch (NumberFormatException e) {
            throw new RuntimeException("Invalid value "+value+" for variable "+variable+". Not a double value.");

        }
    }

    private Integer getIntegerValue(String variable, String value) {

        try {
            return Integer.valueOf(value);
        } catch (NumberFormatException e) {
            throw new RuntimeException("Invalid value "+value+" for variable "+variable+". Not an integer value.");
        }
    }

    private Set<String> getSetValue(String variable, String value) {

        if (value.equals("[]")) {
            return new HashSet<>();
        } else if (value.startsWith("[") && value.endsWith("]") ) {
            String s = value.substring(1,value.length()-1);
            return new HashSet<>(Arrays.asList(s.split(",")));
        } else {
            throw new RuntimeException("Invalid value "+value+" for variable "+variable+". Not a comma separated list value.");
        }
    }

    public String getModificationStr() {
        return modificationStr;
    }

    public Set<Modification> getModifications() {
        return modifications;
    }

    public double getModifMatchMassTol() {
        return modifMatchMassTol;
    }

    public String getCometDecoyProtPrefix() {
        return cometDecoyProtPrefix;
    }

    public Pattern getExcludedProtPattern() {
        return excludedProtPattern;
    }

    public int getMinCharge() {
        return minCharge;
    }

    public int getMaxCharge() {
        return maxCharge;
    }

    public int getMinPeptideLength() {
        return minPeptideLength;
    }

    public int getMaxPeptideLength() {
        return maxPeptideLength;
    }

    public int getMaxRank() {
        return maxRank;
    }

    public int getNrThreads() {
        return nrThreads;
    }

    public double getFdrCometThreshold() {
        return fdrCometThreshold;
    }

    public String getProtCodingGroup() {
        return protCodingGroup;
    }

    public String getNoncanonicalGroup() {
        return noncanonicalGroup;
    }

    public Pattern getSpectrumRegExp() {
        return spectrumRegExp;
    }

    public Pattern getCodingProtRegExp() {
        return codingProtRegExp;
    }

    public Set<String> getForcedNoncanonicalProts() {
        return forcedNoncanonicalProts;
    }

    public int getMinNrPsmsPerHisto() {
        return minNrPsmsPerHisto;
    }

    public String getOutputPrefix() {
        return outputPrefix;
    }

    public int getNrBins3D() {
        return nrBins3D;
    }

    public double getMinXCorr() {
        return minXCorr;
    }

    public double getMaxXCorr() {
        return maxXCorr;
    }

    public int getNrXCorrBins() {
        return nrXCorrBins;
    }

    public double getMinDeltaCn() {
        return minDeltaCn;
    }

    public double getMaxDeltaCn() {
        return maxDeltaCn;
    }

    public int getNrDeltaCnBins() {
        return nrDeltaCnBins;
    }

    public double getMinSpScore() {
        return minSpScore;
    }

    public double getMaxSpScore() {
        return maxSpScore;
    }

    public int getNrSpScoreBins() {
        return nrSpScoreBins;
    }

    public String getFdrControlMethod() {
        return fdrControlMethod;
    }

    public String getVersion() {
        return version;
    }

    public String getCometPsmDir() {
        return cometPsmDir;
    }

    public String getMaxquantPsmDir() {
        return maxquantPsmDir;
    }

    public Pattern getCometPsmRegExp() {
        return cometPsmRegExp;
    }

    public Pattern getMaxquantPsmRegExp() {
        return maxquantPsmRegExp;
    }

    public boolean isIncludeMaxQuant() {
        return includeMaxQuant;
    }

    public String getOutputDir() {
        return outputDir;
    }

    public boolean isReportHistos() {
        return reportHistos;
    }

    public String getSearchFastaFile() {
        return searchFastaFile;
    }

    public String getUniprotFastaFile() {
        return uniprotFastaFile;
    }

    public boolean isDoPeptideProteinGrouping() {
        return doPeptideProteinGrouping;
    }

    public String getWriteParamsFile() {
        return writeParamsFile;
    }

    public String getReadParamsFile() {
        return readParamsFile;
    }

    public int getSmoothDegree() {
        return smoothDegree;
    }

    public String getReadHistos() {
        return readHistos;
    }

    public String getMaxQuantMainScore() {
        return maxQuantMainScore;
    }

    public String getCometMainScore() {
        return cometMainScore;
    }

    public double getMaxQuantMainScoreMinValue() {
        return maxQuantMainScoreMinValue;
    }

    public double getCometMainScoreMinValue() {
        return cometMainScoreMinValue;
    }
}
