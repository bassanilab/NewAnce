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

package newance.psmcombiner;

import newance.psmconverter.PeptideMatchData;
import org.expasy.mzjava.proteomics.mol.Peptide;
import org.expasy.mzjava.proteomics.mol.modification.ModAttachment;
import org.expasy.mzjava.proteomics.mol.modification.Modification;
import org.expasy.mzjava.proteomics.mol.modification.ModificationList;

import java.util.List;
import java.util.function.BiFunction;

/**
 * @author Markus Müller
 */

public class Psm2StringFunction implements BiFunction<String, List<PeptideMatchData>, String> {

    public enum TabStringMode {COMBINED, COMET, MAXQUANT};

    protected final TabStringMode tabStringMode;
    protected final GroupedFDRCalculator groupedFDRCalculator;


    public Psm2StringFunction(TabStringMode tabStringMode, GroupedFDRCalculator groupedFDRCalculator) {
        this.tabStringMode = tabStringMode;
        this.groupedFDRCalculator = groupedFDRCalculator;
    }

    public Psm2StringFunction(TabStringMode tabStringMode) {
        this.tabStringMode = tabStringMode;
        this.groupedFDRCalculator = null;
    }

    @Override
    public String apply(String specID, List<PeptideMatchData> peptideMatchData) {

        String txt = "";
        for (PeptideMatchData psm : peptideMatchData) {
            txt += getTabString(specID, psm);
        }

        return txt;
    }


    public String getHeader() {

        if (tabStringMode == TabStringMode.COMET)
            return "Spectrum\tScanNr\tCharge\tRT\tNeutralMass\tPeptide\tSequence\tPeptideMass\tModifName\tModifPosition\tModifMass\tModifAA\tProteins\t" +
                    "IsVariant\tIsDecoy\tRank\tXCorr\tDeltaCn\tSpScore\tNegLogPv\tmassdiff\ttot_num_ions\tnum_matched_ions\tlFDR";
        if (tabStringMode == TabStringMode.MAXQUANT)
            return "Spectrum\tScanNr\tCharge\tRT\tNeutralMass\tPeptide\tSequence\tPeptideMass\tModifName\tModifPosition\tModifMass\tModifAA\tProteins\t" +
                    "IsVariant\tIsDecoy\tRank\tMass.Error[ppm]\tScore\tDelta.score\tLocalization.prob";
        else
            return "Spectrum\tScanNr\tCharge\tRT\tNeutralMass\tPeptide\tSequence\tPeptideMass\tModifName\tModifPosition\tModifMass\tModifAA\tProteins\t" +
                    "IsVariant\tIsDecoy\tComet.Rank\tComet.XCorr\tComet.DeltaCn\tComet.SpScore\tComet.NegLogPv\tComet.massdiff\tComet.tot_num_ions\t" +
                    "Comet.num_matched_ions\tComet.lFDR\tMaxQuant.Mass.Error[ppm]\tMaxQuant.Score\tMaxQuant.Delta.score\tMaxQuant.Localization.prob";

    }

    public String getTabString(String specID, PeptideMatchData psm) {

        String rt = String.format("%.5f",psm.getScore("rt"));
        String mass = String.format("%.5f",psm.getScore("mass"));

        if (tabStringMode == TabStringMode.COMET)
            return specID+"\t"+(int)psm.getScore("sn")+"\t"+psm.getCharge()+"\t"+rt+"\t"+mass+"\t"+
                    getCometString(psm);
        if (tabStringMode == TabStringMode.MAXQUANT)
            return specID+"\t"+(int)psm.getScore("sn")+"\t"+psm.getCharge()+"\t"+rt+"\t"+mass+"\t"+
                    getMaxQuantString(psm);
        else
            return specID+"\t"+(int)psm.getScore("sn")+"\t"+psm.getCharge()+"\t"+rt+"\t"+mass+"\t"+
                    getCometCombString(psm)+"\t"+getMaxQuantCombString(psm);

    }

    protected String getCometCombString(PeptideMatchData psm) {

        String protACs = psm.getProteinAcc().toString();
        boolean isVariant = (protACs.contains("variant__"));
        int rank = (int) psm.getScore("rank");
        String lfdrStr = (groupedFDRCalculator==null)?"NaN":String.format("%.5f",groupedFDRCalculator.getLocalFDR(psm));
        String expectStr = String.format("%.5f",psm.getScore("neg_log10_p"));
        Peptide peptide = psm.getPeptide();
        String pepMass = String.format("%.5f",peptide.getMolecularMass());
        String modifString = getModifString(peptide);

        return  psm.getPeptide().toString()+"\t"+psm.getPeptide().toSymbolString()+"\t"+pepMass+"\t"+modifString+"\t"+protACs+"\t"+isVariant+"\t"+psm.isDecoy()+"\t"+rank+"\t"+
                psm.getScore("xcorr")+"\t"+psm.getScore("deltacn")+"\t"+psm.getScore("spscore")+"\t"+expectStr+"\t+"+
                psm.getScore("mass_diff")+"\t"+(int)psm.getScore("tot_num_ions")+"\t"+(int)psm.getScore("matched_num_ions")+"\t"+
                lfdrStr;
    }

    protected String getMaxQuantCombString(PeptideMatchData psm) {

        return psm.getScore("mass-error-ppm")+"\t"+psm.getScore("score")+"\t"+psm.getScore("delta-score")+"\t"+
                psm.getScore("mod-pos-prob");
    }

    protected String getCometString(PeptideMatchData psm) {

        String protACs = psm.getProteinAcc().toString();
        boolean isVariant = (protACs.contains("variant__"));
        int rank = (int) psm.getScore("rank");
        String lfdrStr = (groupedFDRCalculator==null)?"NaN":String.format("%.5f",groupedFDRCalculator.getLocalFDR(psm));
        String expectStr = String.format("%.5f",psm.getScore("neg_log10_p"));
        Peptide peptide = psm.getPeptide();
        String pepMass = String.format("%.5f",peptide.getMolecularMass());
        String modifString = getModifString(peptide);

        return  psm.getPeptide().toString()+"\t"+psm.getPeptide().toSymbolString()+"\t"+pepMass+"\t"+modifString+"\t"+protACs+"\t"+isVariant+"\t"+psm.isDecoy()+"\t"+
                rank+"\t"+psm.getScore("xcorr")+"\t"+psm.getScore("deltacn")+"\t"+psm.getScore("spscore")+"\t"+expectStr+"\t+" +
                psm.getScore("mass_diff")+"\t"+(int)psm.getScore("tot_num_ions")+"\t"+(int)psm.getScore("matched_num_ions")+"\t"+
                lfdrStr;
    }

    protected String getMaxQuantString(PeptideMatchData psm) {

        String protACs = psm.getProteinAcc().toString();
        boolean isVariant = (protACs.contains("variant__"));
        int rank = (int) psm.getScore("rank");
        Peptide peptide = psm.getPeptide();
        String pepMass = String.format("%.5f",peptide.getMolecularMass());
        String modifString = getModifString(peptide);

        return  psm.getPeptide().toString()+"\t"+psm.getPeptide().toSymbolString()+"\t"+pepMass+"\t"+modifString+"\t"+protACs+"\t"+isVariant+"\t"+psm.isDecoy()+"\t"+
                rank+"\t"+psm.getScore("mass-error-ppm")+"\t"+psm.getScore("score")+"\t"+psm.getScore("delta-score")+"\t"+
                psm.getScore("mod-pos-prob");
    }

    private String getModifString(Peptide peptide)  {

        if (!peptide.hasModifications()) return "NA\tNA\tNA\tNA";

        String modifNames = "";
        String modifPos = "";
        String modifMass = "";
        String modifAA = "";

        if (peptide.hasModificationAt(ModAttachment.N_TERM)) {
            for (Modification modif : peptide.getModifications(ModAttachment.nTermSet)) {
                modifNames = (modifNames.isEmpty())?modif.getLabel():","+modif.getLabel();
                modifPos = (modifPos.isEmpty())?"0":",0";
                String massStr = String.format("%.5f",modif.getMolecularMass());
                modifMass = (modifMass.isEmpty())?massStr:","+massStr;
                modifAA = (modifAA.isEmpty())?"NT":",NT";
            }
        }

        if (peptide.hasModificationAt(ModAttachment.SIDE_CHAIN)) {
            for (int i : peptide.getModificationIndexes(ModAttachment.sideChainSet)) {

                for (Modification modif : peptide.getModificationsAt(i, ModAttachment.sideChainSet)) {
                    modifNames = (modifNames.isEmpty())?modif.getLabel():","+modif.getLabel();
                    String posStr = String.format("%d",i+1);
                    modifPos = (modifPos.isEmpty())?posStr:","+posStr;
                    String massStr = String.format("%.5f",modif.getMolecularMass());
                    modifMass = (modifMass.isEmpty())?massStr:","+massStr;
                    modifAA = (modifAA.isEmpty())?peptide.getSymbol(i).getSymbol():","+peptide.getSymbol(i).getSymbol();
                }
            }
        }

        if (peptide.hasModificationAt(ModAttachment.C_TERM)) {
            for (Modification modif : peptide.getModifications(ModAttachment.cTermSet)) {
                modifNames = (modifNames.isEmpty())?modif.getLabel():","+modif.getLabel();
                String posStr = String.format("%d",peptide.size()+1);
                modifPos = (modifPos.isEmpty())?posStr:","+posStr;
                String massStr = String.format("%.5f",modif.getMolecularMass());
                modifMass = (modifMass.isEmpty())?massStr:","+massStr;
                modifAA = (modifAA.isEmpty())?"NT":",NT";
            }
        }

        return modifNames+"\t"+modifPos+"\t"+modifMass+"\t"+modifAA;
    }

}
