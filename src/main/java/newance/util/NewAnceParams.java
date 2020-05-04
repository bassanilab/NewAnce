/*
Copyright (C) SIB - Swiss Institute of Bioinformatics, Lausanne, Switzerland
Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
*/

package newance.util;



import newance.mzjava.mol.modification.Modification;

import java.io.*;
import java.nio.BufferUnderflowException;
import java.nio.file.InvalidPathException;
import java.util.*;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

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
    private boolean debug = false;

    private double fdrCometThreshold = 0.03;

    private String groupingMethod = "none";
    private List<String> groupNames = new ArrayList<>();
    private List<Pattern> groupRegExs = new ArrayList<>();
    private String protCodingGroup = "prot";
    private String noncanonicalGroup = "nonc";
    private Pattern spectrumRegExp = null;
    private Pattern codingProtRegExp = null;
    private String proteinGroupMapFile = "";
    private Map<String,Set<String>> proteinGroupMap = new HashMap<>();
    private boolean writeCometTabFile = false;

    private int minNrPsmsPerHisto = 100000;

    private String outputTag = "";

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

    private String fdrControlMethod = "combined";

    private String version = "version 1.6.0";

    private String cometPsmDir = "";
    private String maxquantPsmDir = "";
    private Pattern cometPsmRegExp = null;
    private Pattern maxquantPsmRegExp = null;
    private boolean includeMaxQuant = false;
    private String outputDir = "";
    private boolean reportHistos = false;
    private String readHistos = "";
    private boolean forceHistos = false;
    private String searchFastaFile = "";
    private String uniprotFastaFile= "";
    private boolean doPeptideProteinGrouping = false;
    private boolean writeParamsFile = false;
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
        res +=  "excludedProtPattern="+((excludedProtPattern==null)?"":excludedProtPattern.toString())+"\n";
        res +=  "minCharge="+minCharge+"\n";
        res +=  "maxCharge="+maxCharge+"\n";
        res +=  "minPeptideLength="+minPeptideLength+"\n";
        res +=  "maxPeptideLength="+maxPeptideLength+"\n";
        res +=  "maxRank="+maxRank+"\n";
        res +=  "debug="+debug+"\n";
        res +=  "smoothDegree="+smoothDegree+"\n";
        res +=  "nrThreads="+nrThreads+"\n";
        res +=  "fdrCometThreshold="+fdrCometThreshold+"\n";
        res +=  "protCodingGroup="+ protCodingGroup +"\n";
        res +=  "noncanonicalGroup="+ noncanonicalGroup +"\n";
        res +=  "spectrumRegExp="+spectrumRegExp+"\n";
        res +=  "codingProtRegExp="+ codingProtRegExp +"\n";
        res +=  "proteinGroupMapFile="+ proteinGroupMapFile +"\n";
        res +=  "minNrPsmsPerHisto="+minNrPsmsPerHisto+"\n";
        res +=  "outputTag="+ outputTag +"\n";
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
        res +=  "forceHistos="+forceHistos+"\n";
        res +=  "searchFastaFile="+searchFastaFile+"\n";
        res +=  "uniprotFastaFile="+uniprotFastaFile+"\n";
        res +=  "doPeptideProteinGrouping="+doPeptideProteinGrouping+"\n";
        res +=  "maxQuantMainScore="+maxQuantMainScore+"\n";
        res +=  "cometMainScore="+cometMainScore+"\n";
        res +=  "maxQuantMainScoreMinValue="+maxQuantMainScoreMinValue+"\n";
        res +=  "cometMainScoreMinValue="+cometMainScoreMinValue+"\n";
        res +=  "groupingMethod="+groupingMethod+"\n";
        res +=  "groupNames="+groupNames+"\n";
        res +=  "groupRegExs="+groupRegExs+"\n";
        res +=  "writeCometTabFile="+writeCometTabFile+"\n";

        res +=  "writeParamsFile="+writeParamsFile+"\n";

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
                if (fields.length==2 && !fields[1].trim().equals("null") && !fields[1].trim().isEmpty()) variableValueMap.put(fields[0].trim(),fields[1].trim());
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
            Set<Modification> modifs = new HashSet<>();
            try {
                for (String modif : getSetValue("modifications", variableValueMap.get("modifications"))) {
                    buf = modif;
                    if (!modif.trim().isEmpty()) modifs.add(Modification.parseModification(modif));
                }
            } catch(IllegalArgumentException e) {
                throw new RuntimeException("Invalid modification string "+buf+" provided. Abort");
            }

            if (!modifs.isEmpty()) {
                modifications.clear();
                modifications.addAll(modifs);
            }
        }

        if (variableValueMap.containsKey("modifMatchMassTol")) {
            modifMatchMassTol = getDoubleValue("modifMatchMassTol", variableValueMap.get("modifMatchMassTol"), 0.0, Double.MAX_VALUE);
        }

        if (variableValueMap.containsKey("cometDecoyProtPrefix")) {
            cometDecoyProtPrefix = getStringValue("cometDecoyProtPrefix",variableValueMap.get("cometDecoyProtPrefix"));
        }

        if (variableValueMap.containsKey("excludedProtPattern")) {
            excludedProtPattern = getPatternValue("excludedProtPattern",variableValueMap.get("excludedProtPattern"));
        }

        if (variableValueMap.containsKey("minCharge")) {
            minCharge = getIntegerValue("minCharge",variableValueMap.get("minCharge"), 1, 10);
        }

        if (variableValueMap.containsKey("maxCharge")) {
            maxCharge = getIntegerValue("maxCharge",variableValueMap.get("maxCharge"), 1, 10);
        }

        if (variableValueMap.containsKey("minPeptideLength")) {
            minPeptideLength = getIntegerValue("minPeptideLength",variableValueMap.get("minPeptideLength"), 4, 1000);
        }

        if (variableValueMap.containsKey("maxPeptideLength")) {
            maxPeptideLength = getIntegerValue("maxPeptideLength",variableValueMap.get("maxPeptideLength"), 4, 1000);
        }

        if (variableValueMap.containsKey("maxRank")) {
            maxRank = getIntegerValue("maxRank",variableValueMap.get("maxRank"), 1, 100);
        }

        if (variableValueMap.containsKey("nrThreads")) {
            nrThreads = getIntegerValue("nrThreads",variableValueMap.get("nrThreads"), 1, 256);
        } else {
            int nrProc = Runtime.getRuntime().availableProcessors();
            nrThreads = (nrProc>2)?nrProc-2:1;
            nrThreads = (nrThreads>256)?256:nrThreads;
        }

        if (variableValueMap.containsKey("smoothDegree")) {
            smoothDegree = getIntegerValue("smoothDegree",variableValueMap.get("smoothDegree"), 0, 100);
        }

        if (variableValueMap.containsKey("fdrCometThreshold")) {
            fdrCometThreshold = getDoubleValue("fdrCometThreshold",variableValueMap.get("fdrCometThreshold"), 0.0, Double.MAX_VALUE);
        }

        if (variableValueMap.containsKey("groupingMethod")) {
            groupingMethod = getStringValue("groupingMethod",variableValueMap.get("groupingMethod"),new HashSet<>(Arrays.asList(new String[]{"fasta","modif","none"})));
        }

        if (variableValueMap.containsKey("groupNames")) {
            groupNames = getListValue("groupNames",variableValueMap.get("groupNames"));

            if (hasDuplicates(groupNames)) throw new RuntimeException("Duplicates in groupNames: "+groupNames+". Abort.");
            if (hasNA(groupNames)) throw new RuntimeException("Empty items in groupNames: "+groupNames+". Abort.");
        }

        if (variableValueMap.containsKey("groupRegExs")) {
            List<String> groupRegExStr = getListValue("groupRegExs",variableValueMap.get("groupRegExs"));
            if (hasDuplicates(groupRegExStr)) throw new RuntimeException("Duplicates in groupRegExs: "+groupRegExStr+". Abort.");
            if (hasNA(groupRegExStr)) throw new RuntimeException("Empty items in groupRegExs: "+groupRegExStr+". Abort.");

            groupRegExs = new ArrayList<>();
            for (String re : groupRegExStr) {
                try {
                    groupRegExs.add(Pattern.compile(re));
                } catch(PatternSyntaxException e) {
                    throw new RuntimeException("Item in groupRegExs not valid RegExp: "+re+". Abort.");
                }
            }

        }

        if (variableValueMap.containsKey("writeCometTabFile")) {
            writeCometTabFile = getBooleanValue("writeCometTabFile",variableValueMap.get("writeCometTabFile"));
        }

        if (variableValueMap.containsKey("codingProtRegExp")) {
            codingProtRegExp = getPatternValue("codingProtRegExp",variableValueMap.get("codingProtRegExp"));
            groupRegExs.add(codingProtRegExp);
            groupingMethod = "fasta";
        }

        // leave this for back compatibility
        if (codingProtRegExp !=null) {
            if (variableValueMap.containsKey("protCodingGroup")) {
                protCodingGroup = getStringValue("protCodingGroup",variableValueMap.get("protCodingGroup"));
                groupNames.add(protCodingGroup);
            } else {
                throw new RuntimeException("No name for protein coding group provided. Abort");
            }

            if (variableValueMap.containsKey("noncanonicalGroup")) {
                noncanonicalGroup = getStringValue("noncanonicalGroup",variableValueMap.get("noncanonicalGroup"));
                groupNames.add(noncanonicalGroup);
            } else {
                throw new RuntimeException("No name for non-canonical group provided. Abort");
            }
        }

        if (variableValueMap.containsKey("spectrumRegExp")) {
            spectrumRegExp = getPatternValue("spectrumRegExp",variableValueMap.get("spectrumRegExp"));
        }

        if (variableValueMap.containsKey("proteinGroupMapFile")) {
            proteinGroupMapFile = getFileValue("proteinGroupMapFile", variableValueMap.get("proteinGroupMapFile"));
            proteinGroupMap = readProteinGroupMapFile(proteinGroupMapFile);
        }

        if (variableValueMap.containsKey("minNrPsmsPerHisto")) {
            minNrPsmsPerHisto = getIntegerValue("minNrPsmsPerHisto",variableValueMap.get("minNrPsmsPerHisto"), 0, Integer.MAX_VALUE);
        }

        if (variableValueMap.containsKey("outputTag")) {
            outputTag = variableValueMap.get("outputTag"); // can be empty
        }

        if (variableValueMap.containsKey("minXCorr")) {
            minXCorr = getDoubleValue("minXCorr",variableValueMap.get("minXCorr"), 0.0, Double.MAX_VALUE);
        }

        if (variableValueMap.containsKey("maxXCorr")) {
            maxXCorr = getDoubleValue("maxXCorr",variableValueMap.get("maxXCorr"), 0.0, Double.MAX_VALUE);
        }

        if (variableValueMap.containsKey("nrXCorrBins")) {
            nrXCorrBins = getIntegerValue("nrXCorrBins",variableValueMap.get("nrXCorrBins"), 1, 100000);
        }

        if (variableValueMap.containsKey("minDeltaCn")) {
            minDeltaCn = getDoubleValue("minDeltaCn",variableValueMap.get("minDeltaCn"), 0.0, Double.MAX_VALUE);
        }

        if (variableValueMap.containsKey("maxDeltaCn")) {
            maxDeltaCn = getDoubleValue("maxDeltaCn",variableValueMap.get("maxDeltaCn"), 0.0, Double.MAX_VALUE);
        }

        if (variableValueMap.containsKey("nrDeltaCnBins")) {
            nrDeltaCnBins = getIntegerValue("nrDeltaCnBins",variableValueMap.get("nrDeltaCnBins"), 1, 100000);
        }

        if (variableValueMap.containsKey("minSpScore")) {
            minSpScore = getDoubleValue("minSpScore",variableValueMap.get("minSpScore"), 0.0, Double.MAX_VALUE);
        }

        if (variableValueMap.containsKey("maxSpScore")) {
            maxSpScore = getDoubleValue("maxSpScore",variableValueMap.get("maxSpScore"), 0.0, Double.MAX_VALUE);
        }

        if (variableValueMap.containsKey("nrSpScoreBins")) {
            nrSpScoreBins = getIntegerValue("nrSpScoreBins",variableValueMap.get("nrSpScoreBins"), 1, 100000);
        }

        if (variableValueMap.containsKey("fdrControlMethod")) {
            fdrControlMethod = getStringValue("fdrControlMethod",variableValueMap.get("fdrControlMethod").trim().toLowerCase(),
            new HashSet<>(Arrays.asList(new String[]{"combined","separate"})));
        }

        if (variableValueMap.containsKey("cometPsmDir")) {
            cometPsmDir = getDirectoryValue("cometPsmDir",variableValueMap.get("cometPsmDir"));
        }

        if (variableValueMap.containsKey("cometPsmRegExp")) {
            cometPsmRegExp = getPatternValue("cometPsmRegExp",variableValueMap.get("cometPsmRegExp"));
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
        }

        if (variableValueMap.containsKey("reportHistos")) {
            reportHistos = getBooleanValue("reportHistos",variableValueMap.get("reportHistos"));
        }

        if (variableValueMap.containsKey("debug")) {
            debug = getBooleanValue("debug",variableValueMap.get("debug"));
        }

        if (variableValueMap.containsKey("readHistos")) {
            readHistos = getDirectoryValue("readHistos",variableValueMap.get("readHistos"));
        }

        if (variableValueMap.containsKey("forceHistos")) {
            forceHistos = getBooleanValue("forceHistos",variableValueMap.get("forceHistos"));
        }

        if (variableValueMap.containsKey("searchFastaFile")) {
            searchFastaFile = getFileValue("searchFastaFile",variableValueMap.get("searchFastaFile"));
        }

        if (variableValueMap.containsKey("uniprotFastaFile")) {
            uniprotFastaFile = getFileValue("uniprotFastaFile",variableValueMap.get("uniprotFastaFile"));
        }

        if (variableValueMap.containsKey("doPeptideProteinGrouping")) {
            doPeptideProteinGrouping = getBooleanValue("doPeptideProteinGrouping",variableValueMap.get("doPeptideProteinGrouping"));
        }

        if (variableValueMap.containsKey("writeParamsFile")) {
            writeParamsFile = getBooleanValue("writeParamsFile",variableValueMap.get("writeParamsFile"));
        }

        if (variableValueMap.containsKey("readParamsFile")) {
            readParamsFile = getFileValue("readParamsFile",variableValueMap.get("readParamsFile"));
        }

        if (variableValueMap.containsKey("maxQuantMainScoreMinValue")) {
            maxQuantMainScoreMinValue = getDoubleValue("maxQuantMainScoreMinValue",variableValueMap.get("maxQuantMainScoreMinValue"), Double.MIN_VALUE, Double.MAX_VALUE);
        }

        if (variableValueMap.containsKey("cometMainScoreMinValue")) {
            cometMainScoreMinValue = getDoubleValue("cometMainScoreMinValue",variableValueMap.get("cometMainScoreMinValue"), Double.MIN_VALUE, Double.MAX_VALUE);
        }

        if (variableValueMap.containsKey("maxQuantMainScore")) {
            maxQuantMainScore = getStringValue("maxQuantMainScore",variableValueMap.get("maxQuantMainScore"),new HashSet<>(Arrays.asList(new String[]{"Score","Delta score","Localization prob","PEP"})));
        }

        if (variableValueMap.containsKey("cometMainScore")) {
            cometMainScore = getStringValue("cometMainScore",variableValueMap.get("cometMainScore"),new HashSet<>(Arrays.asList(new String[]{"xcorr","spscore","deltacn","expect"})));
        }

        checkVariableValues();
    }

    private boolean checkVariableValues() {

        if (minXCorr>=maxXCorr) throw new RuntimeException("Error in minXCorr value: "+minXCorr+". Abort.");
        if (minDeltaCn>=maxDeltaCn) throw new RuntimeException("Error in minDeltaCn value: "+minDeltaCn+". Abort.");
        if (minSpScore>=maxSpScore) throw new RuntimeException("Error in minSpScore value: "+minSpScore+". Abort.");

        if (doPeptideProteinGrouping && !(new File(searchFastaFile)).exists()) throw new RuntimeException("Error in doPeptideProteinGrouping value: "+doPeptideProteinGrouping+". Abort.");

        if (groupingMethod.equals("fasta")) {

            if (groupNames.size()!=groupRegExs.size()+1) throw new RuntimeException("Incompatible size of groupNames and groupRegExs : "+groupNames+" - "+groupRegExs+". Abort.");

        }

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

    public static Set<Modification> getModifications(String variable, String value) {

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

    public static String getStringValue(String variable, String value, Set<String> allowedValues) {

        if (allowedValues.contains(value)) {
            return value;
        } else {
            throw new RuntimeException("Invalid value "+value+" for variable "+variable+". Value has to be one of "+allowedValues.toString());
        }
    }

    public static String getStringValue(String variable, String value) {

        if (!value.isEmpty()) {
            return value;
        } else {
            throw new RuntimeException("Invalid value "+value+" for variable "+variable+". Value is empty");
        }
    }

    public static Pattern getPatternValue(String variable, String value) {

        if (value.isEmpty()) return null;

        try {
            return Pattern.compile(value);
        } catch (PatternSyntaxException e) {
            throw new RuntimeException("Invalid value "+value+" for variable "+variable+". Not a valid regular expression.");
        }
    }

    public static boolean getBooleanValue(String variable, String value) {

        if (value.toLowerCase().equals("true")) {
            return true;
        } else if (value.toLowerCase().equals("false")) {
            return false;
        } else {
            throw new RuntimeException("Invalid value "+value+" for variable "+variable+". Value has to be either true or false.");
        }
    }

    public static String getDirectoryValue(String variable, String value) {

        File dir = new File(value);
        if (!dir.exists() && !dir.isDirectory()) {
            throw new InvalidPathException("Invalid value "+value+" for variable "+variable+"."," Directory does not exist or is not a directory.");
        }

        return value;
    }

    public static String getNewDirectoryValue(String variable, String value) {

        File dir = new File(value);

        if (!dir.exists()) {
            dir.mkdirs();
        }

        if (!dir.isDirectory()) {
            throw new InvalidPathException("Invalid value "+value+" for variable "+variable+"."," Is not a directory.");
        }

        return value;
    }

    public static String getFileValue(String variable, String value) {

        File file = new File(value);
        if (!file.exists()) {
            throw new InvalidPathException("Invalid value "+value+" for variable "+variable+"."," File does not exist.");
        }

        return value;
    }

    public static Double getDoubleValue(String variable, String value, double min, double max) {

        Double d;
        try {
            d = Double.valueOf(value);

            if (d<min || d>max) {
                throw new RuntimeException("Invalid value "+value+" for variable "+variable+". Value out of range.");
            }

            return d;
        } catch (NumberFormatException e) {
            throw new RuntimeException("Invalid value "+value+" for variable "+variable+". Not a double value.");

        }
    }

    public static Integer getIntegerValue(String variable, String value, int min, int max) {

        Integer i;
        try {
            i = Integer.valueOf(value);

            if (i<min || i>max) {
                throw new RuntimeException("Invalid value "+value+" for variable "+variable+". Value out of range.");
            }

            return i;
        } catch (NumberFormatException e) {
            throw new RuntimeException("Invalid value "+value+" for variable "+variable+". Not an integer value.");
        }
    }

    public static Set<String> getSetValue(String variable, String value) {

        int i = 0;
        while (value.charAt(i)=='[' && i<value.length()) i++;
        int j = value.length()-1;
        while (value.charAt(j)==']' && j>=0) j--;

        int l = (i<value.length()-1-j)?i:value.length()-1-j;
        i = l;
        j = value.length()-l;

        String listStr = value.substring(i,j);

        if (listStr.isEmpty()) {
            return new HashSet<>();
        } else {
            return new HashSet<>(Arrays.asList(listStr.split(",")));
        }
    }

    public static List<String> getListValue(String variable, String value) {

        int i = 0;
        while (value.charAt(i)=='[' && i<value.length()) i++;
        int j = value.length()-1;
        while (value.charAt(j)==']' && j>=0) j--;

        int l = (i<value.length()-1-j)?i:value.length()-1-j;
        i = l;
        j = value.length()-l;

        String listStr = value.substring(i,j);

        if (listStr.isEmpty()) {
            return new ArrayList<>();
        } else {
            return new ArrayList<>(Arrays.asList(listStr.split(",")));
        }
    }

    protected static boolean hasDuplicates(List<String> list) {


        for (int i=0;i<list.size()-1;i++) {

            String first = list.get(i);
            for (int j=i+1;i<list.size();i++) {
                if (first.equals(list.get(j))) return true;
            }
        }

        return false;
    }


    protected static boolean hasNA(List<String> list) {

        for (int i=0;i<list.size();i++) {
            if (list.get(i).isEmpty()) return true;
        }

        return false;
    }

    protected static Map<String,Set<String>> readProteinGroupMapFile(String fileName) {

        Map<String,Set<String>> proteinGroupMap = new HashMap<>();

        try {
            BufferedReader reader = new BufferedReader(new FileReader(new File(fileName)));

            String line = "";

            while ((line = reader.readLine())!=null) {
                String[] fields = line.split("\t");

                proteinGroupMap.putIfAbsent(fields[1],new HashSet<>());
                proteinGroupMap.get(fields[1]).add(fields[0]);
            }


        } catch (FileNotFoundException e) {

        } catch (IOException e) {

        }

        return proteinGroupMap;
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

    public Map<String,Set<String>> getProteinGroupMap() {
        return proteinGroupMap;
    }

    public int getMinNrPsmsPerHisto() {
        return minNrPsmsPerHisto;
    }

    public String getOutputTag() {
        return outputTag;
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

    public boolean isWriteParamsFile() {
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

    public boolean isForceHistos() {
        return forceHistos;
    }

    public boolean isDebug() {
        return debug;
    }

    public String getGroupingMethod() {
        return groupingMethod;
    }

    public List<String> getGroupNames() {
        return groupNames;
    }

    public List<Pattern> getGroupRegExs() {
        return groupRegExs;
    }

    public boolean isWriteCometTabFile() {
        return writeCometTabFile;
    }
}
