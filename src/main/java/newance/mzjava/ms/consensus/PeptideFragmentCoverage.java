package newance.mzjava.ms.consensus;

import newance.mzjava.mol.IonType;
import newance.mzjava.mol.Peptide;
import newance.mzjava.ms.peaklist.PeakList;
import newance.mzjava.ms.spectrum.PepFragAnnotation;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by markusmueller on 05.05.21.
 */
public class PeptideFragmentCoverage {

    protected PeakList<PepFragAnnotation> spectrum;
    protected Peptide peptide;
    protected String peptideFragInfo = null;
    protected String peptideFragInfoShort = null;
    protected String peptideFragInfoLong = null;
    protected double sequenceCoverage = -1.0;
    protected double spectrumCoverage = -1.0;

    public PeptideFragmentCoverage(PeakList<PepFragAnnotation> spectrum, Peptide peptide) {

        this.peptide = peptide;
        this.spectrum = spectrum;
    }

    public String getPeptideFragInfo() {
        if (peptideFragInfo==null) calcPeptideFragInfo();

        return peptideFragInfo;
    }

    public String getPeptideFragInfoShort() {
        if (peptideFragInfoShort==null) calcPeptideFragInfo();

        return peptideFragInfoShort;
    }

    public String getPeptideFragInfoLong() {
        if (peptideFragInfoLong==null) calcPeptideFragInfo();

        return peptideFragInfoLong;
    }

    public double getSequenceCoverage() {
        if (sequenceCoverage < 0) calcPeptideFragInfo();

        return sequenceCoverage;
    }

    public double getSpectrumCoverage() {
        if (spectrumCoverage < 0) calcPeptideFragInfo();

        return spectrumCoverage;
    }

    protected void calcPeptideFragInfo() {


        double totIntensity = spectrum.getTotalIonCurrent();
        double maxIntensity = spectrum.getIntensity(spectrum.getMostIntenseIndex());
        double totMatchedIntensity = 0.0;
        List<List<PepFragAnnotation>> peptideFrags = new ArrayList<>();
        for (int i=0;i<=peptide.size();i++){
            peptideFrags.add(new ArrayList<>());
        }
        List<Double> matchedIntensitiesY = new ArrayList<>();
        for (int i=0;i<=peptide.size();i++){
            matchedIntensitiesY.add(0.0);
        }
        List<Double> matchedIntensitiesB = new ArrayList<>();
        for (int i=0;i<=peptide.size();i++){
            matchedIntensitiesB.add(0.0);
        }

        peptideFragInfoLong = "";
        for (int i=0;i<spectrum.size();i++) {
            if (spectrum.hasAnnotationsAt(i)) {
                List<PepFragAnnotation> annots = spectrum.getAnnotations(i);
                if (!peptideFragInfoLong.isEmpty()) peptideFragInfoLong += ",";
                peptideFragInfoLong +=
                        String.format("(%.5f,%.5f: ",spectrum.getMz(i),spectrum.getIntensity(i)/maxIntensity);
                String annotList = "";
                for (PepFragAnnotation annot : annots) {
                    if (!annotList.isEmpty()) annotList += ",";
                    annotList += annot.getIonType()+"|"+annot.getCharge()+"|"+annot.getFragment()+"|"+
                            annot.getNeutralLoss().getFormula();
                }
                peptideFragInfoLong += annotList+")";
            }
        }

        for (int i=0;i<spectrum.size();i++) {
            if (spectrum.hasAnnotationsAt(i)) {

                double h = spectrum.getIntensity(i);
                totMatchedIntensity += h;
                h /= maxIntensity;

                List<PepFragAnnotation> annots = spectrum.getAnnotations(i);
                for (PepFragAnnotation annot : annots) {
                    IonType ionType = annot.getIonType();
                    String fragSeq = annot.getFragment().toString();
                    if (ionType==IonType.a || ionType==IonType.b) {
                        int pos = fragSeq.length()-1;
                        peptideFrags.get(pos).add(annot);
                        if (h>matchedIntensitiesB.get(pos)) matchedIntensitiesB.set(pos,h);
                    }
                    if (ionType==IonType.y) {
                        int pos = peptide.size()-fragSeq.length()-1;
                        peptideFrags.get(pos).add(annot);
                        if (h>matchedIntensitiesY.get(pos)) matchedIntensitiesY.set(pos,h);
                    }
                }
            }
        }

        spectrumCoverage = totMatchedIntensity/totIntensity;

        sequenceCoverage = 0.0;
        peptideFragInfoShort = "";
        peptideFragInfo = "";
        for (int i=0;i<peptideFrags.size()-1;i++) {
            peptideFragInfo += peptide.getSymbol(i);
            peptideFragInfoShort += peptide.getSymbol(i);
            List<PepFragAnnotation> annots = peptideFrags.get(i);
            if (annots.size()>0) {
                sequenceCoverage += 1.0;
                int intensityRangeY = (int)Math.floor(matchedIntensitiesY.get(i)*10);
                int intensityRangeB = (int)Math.floor(matchedIntensitiesB.get(i)*10);
                peptideFragInfo += "(";
                boolean y = false;
                boolean b = false;
                for (int j=0;j<annots.size();j++){
                    IonType ionType = annots.get(j).getIonType();
                    if (ionType==IonType.a || ionType==IonType.b) b = true;
                    if (ionType==IonType.y) y = true;
                }
                if (y) peptideFragInfo += "y"+intensityRangeY;
                if (b) peptideFragInfo += "b"+intensityRangeB;
                peptideFragInfo += ")";
                peptideFragInfoShort += "|";
            }
        }

        sequenceCoverage /= (peptide.size()-1);
    }
}
