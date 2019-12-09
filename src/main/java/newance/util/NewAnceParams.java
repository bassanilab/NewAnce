package newance.util;


import java.io.Serializable;
import java.util.Set;
import java.util.regex.Pattern;

/**
 * @author Markus Muller
 * @version 0.0
 */

public class NewAnceParams implements Serializable {

    private static NewAnceParams instance = null;

    // When reading PSM results files the modifications are mapped to the modifications defined here. Modifications are mapped by means of there mass only
    // (this will change in future versions, where the mapping will be done by mass and amino acid)

    private String modificationStr = "Cysteinyl:C3H5NO2S,Oxidation:O,Carbamidomethyl:C2H3NO,Phospho:HO3P,Deamidated:H-1N-1O," +
            "Ammonia-loss:H-3N-1,Acetyl:C2H2O,Methyl:CH2,Dimethyl:C2H4,Trimethyl:C3H6,Amidated:H-1,Glu->pyro-Glu:H-2O-1," +
            "Propionyl:C3H4O,Pyro-carbamidomethyl:C2O";

    // Mass tolerance for the modification mapping. This value is independent of instrument precision, but refers to the precision of the modification mass values
    // indicated in the PSM result files.
    private double modifMatchMassTol =  0.015;

    // Remove decoy proteins from list of psm proteins for every PSM. remove PSM if there are only decoy proteins.
    private String decoyProtPattern = "DECOY_";

    // Remove excluded proteins from list of psm proteins for every PSM. remove PSM if there are only excluded proteins.
    private String excludedProtPattern = "";

    // Predicate to filter PeptideMatchData objects before writing them to hadoop file
    private int minCharge = 1;
    private int maxCharge = 5;
    private int minPeptideLength = 8;
    private int maxPeptideLength = 25;
    private int maxRank = 2;

    private int nrThreads = 1;

    private double fdrCometThreshold = 0.03;

    private String standartGroup;
    private String crypticGroup;
    private Pattern spectrumRegExp;
    private Pattern stdProtRegExp;
    private Set<String> forcedCrypticProts;

    private int minNrPsmsPerHisto = 50000;

    private String outputPrefix;

    private final int nrBins3D = 40;
    private double minXCorr = 0.0;
    private double maxXCorr = 5.0;
    private int nrXCorrBins = nrBins3D;
    private double minDeltaCn = 0.0;
    private double maxDeltaCn = 1.0;
    private int nrDeltaCnBins = nrBins3D;
    private double minSpScore = 0.0;
    private double maxSpScore = 2500.0;
    private int nrSpScoreBins = nrBins3D;

    private String fdrControlMethod;

    private String version = "version 1.3.0";


    public static NewAnceParams getInstance() {
        if (instance==null) {
            instance = new NewAnceParams();
        }

        return instance;
    }

    private NewAnceParams() {
    }

    public String toString() {
        String res = "";

        res +=  "modificationStr="+modificationStr+"\n";
        res +=  "modifMatchMassTol="+modifMatchMassTol+"\n";

        res +=  "decoyProtPattern="+decoyProtPattern+"\n";
        res +=  "excludedProtPattern="+excludedProtPattern+"\n";
        res +=  "minCharge="+minCharge+"\n";
        res +=  "maxCharge="+maxCharge+"\n";
        res +=  "minPeptideLength="+minPeptideLength+"\n";
        res +=  "maxPeptideLength="+maxPeptideLength+"\n";
        res +=  "maxRank="+maxRank+"\n";
        res +=  "nrThreads="+nrThreads+"\n";
        res +=  "fdrCometThreshold="+fdrCometThreshold+"\n";
        res +=  "standartGroup="+standartGroup+"\n";
        res +=  "crypticGroup="+crypticGroup+"\n";
        res +=  "spectrumRegExp="+spectrumRegExp+"\n";
        res +=  "stdProtRegExp="+stdProtRegExp+"\n";
        res +=  "forcedCrypticProts="+forcedCrypticProts+"\n";
        res +=  "minNrPsmsPerHisto="+minNrPsmsPerHisto+"\n";
        res +=  "outputPrefix="+outputPrefix+"\n";
        res +=  "nrBins3D="+nrBins3D+"\n";
        res +=  "minXCorr="+minXCorr+"\n";
        res +=  "maxXCorr="+maxXCorr+"\n";
        res +=  "nrXCorrBins="+nrXCorrBins+"\n";
        res +=  "minDeltaCn="+minDeltaCn+"\n";
        res +=  "maxDeltaCn="+maxDeltaCn+"\n";
        res +=  "nrDeltaCnBins="+nrDeltaCnBins+"\n";
        res +=  "minSpScore="+minSpScore+"\n";
        res +=  "maxSpScore="+maxSpScore+"\n";
        res +=  "nrSpScoreBins="+nrSpScoreBins+"\n";

        return res;
    }

    public String getDecoyProtPattern() {return decoyProtPattern;}

    public String getExcludedProtPattern() {return excludedProtPattern;}

    public String getModificationStr() {
        return modificationStr;
    }

    public double getModifMatchMassTol() {
        return modifMatchMassTol;
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

    public void setModificationStr(String modificationStr) {
        this.modificationStr = modificationStr;
    }

    public void setModifMatchMassTol(double modifMatchMassTol) {
        this.modifMatchMassTol = modifMatchMassTol;
    }

    public void setDecoyProtPattern(String decoyProtPattern) {
        this.decoyProtPattern = decoyProtPattern;
    }

    public void setExcludedProtPattern(String excludedProtPattern) {
        this.excludedProtPattern = excludedProtPattern;
    }

    public void setMinCharge(int minCharge) {
        this.minCharge = minCharge;
    }

    public void setMaxCharge(int maxCharge) {
        this.maxCharge = maxCharge;
    }

    public void setMinPeptideLength(int minPeptideLength) {
        this.minPeptideLength = minPeptideLength;
    }

    public void setMaxPeptideLength(int maxPeptideLength) {
        this.maxPeptideLength = maxPeptideLength;
    }

    public void setMaxRank(int maxRank) {
        this.maxRank = maxRank;
    }

    public double getFdrCometThreshold() {
        return fdrCometThreshold;
    }

    public void setFdrCometThreshold(double fdrCometThreshold) {
        this.fdrCometThreshold = fdrCometThreshold;
    }

    public String getStandartGroup() {
        return standartGroup;
    }

    public void setStandartGroup(String standartGroup) {
        this.standartGroup = standartGroup;
    }

    public String getCrypticGroup() {
        return crypticGroup;
    }

    public void setCrypticGroup(String crypticGroup) {
        this.crypticGroup = crypticGroup;
    }

    public Pattern getSpectrumRegExp() {
        return spectrumRegExp;
    }

    public void setSpectrumRegExp(Pattern spectrumRegExp) {
        this.spectrumRegExp = spectrumRegExp;
    }

    public Pattern getStdProtRegExp() {
        return stdProtRegExp;
    }

    public void setStdProtRegExp(Pattern stdProtRegExp) {
        this.stdProtRegExp = stdProtRegExp;
    }

    public Set<String> getForcedCrypticProts() {
        return forcedCrypticProts;
    }

    public void setForcedCrypticProts(Set<String> excludeHumanProts) {
        this.forcedCrypticProts = excludeHumanProts;
    }

    public int getNrThreads() {
        return nrThreads;
    }

    public void setNrThreads(int nrThreads) {
        this.nrThreads = nrThreads;
    }

    public int getNrBins3D() {
        return nrBins3D;
    }

    public double getMinXCorr() {
        return minXCorr;
    }

    public void setMinXCorr(double minXCorr) {
        this.minXCorr = minXCorr;
    }

    public double getMaxXCorr() {
        return maxXCorr;
    }

    public void setMaxXCorr(double maxXCorr) {
        this.maxXCorr = maxXCorr;
    }

    public int getNrXCorrBins() {
        return nrXCorrBins;
    }

    public void setNrXCorrBins(int nrXCorrBins) {
        this.nrXCorrBins = nrXCorrBins;
    }

    public double getMinDeltaCn() {
        return minDeltaCn;
    }

    public void setMinDeltaCn(double minDeltaCn) {
        this.minDeltaCn = minDeltaCn;
    }

    public double getMaxDeltaCn() {
        return maxDeltaCn;
    }

    public void setMaxDeltaCn(double maxDeltaCn) {
        this.maxDeltaCn = maxDeltaCn;
    }

    public int getNrDeltaCnBins() {
        return nrDeltaCnBins;
    }

    public void setNrDeltaCnBins(int nrDeltaCnBins) {
        this.nrDeltaCnBins = nrDeltaCnBins;
    }

    public double getMinSpScore() {
        return minSpScore;
    }

    public void setMinSpScore(double minSpScore) {
        this.minSpScore = minSpScore;
    }

    public double getMaxSpScore() {
        return maxSpScore;
    }

    public void setMaxSpScore(double maxSpScore) {
        this.maxSpScore = maxSpScore;
    }

    public int getNrSpScoreBins() {
        return nrSpScoreBins;
    }

    public void setNrSpScoreBins(int nrSpScoreBins) {
        this.nrSpScoreBins = nrSpScoreBins;
    }

    public String getOutputPrefix() {
        return outputPrefix;
    }

    public void setOutputPrefix(String outputPrefix) {
        this.outputPrefix = outputPrefix;
    }

    public int getMinNrPsmsPerHisto() {
        return minNrPsmsPerHisto;
    }

    public void setMinNrPsmsPerHisto(int minNrPsmsPerHisto) {
        this.minNrPsmsPerHisto = minNrPsmsPerHisto;
    }

    public String getFdrControlMethod() {
        return fdrControlMethod;
    }

    public void setFdrControlMethod(String fdrControlMethod) {
        this.fdrControlMethod = fdrControlMethod;
    }

    public String getVersion() {
        return version;
    }

    public void setVersion(String version) {
        this.version = version;
    }
}
