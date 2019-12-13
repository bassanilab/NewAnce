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

import java.util.List;
import java.util.function.BiFunction;

/**
 * @author Markus Müller
 */

public class Psm2StringFunction implements BiFunction<String, List<PeptideMatchData>, String> {

    public enum TabStringMode {COMBINED, COMET, MAXQUANT};

    protected final TabStringMode tabStringMode;

    final String mode;

    public Psm2StringFunction(String mode, TabStringMode tabStringMode) {
        this.mode = mode;
        this.tabStringMode = tabStringMode;
    }

    public Psm2StringFunction(String mode) {
        this.mode = mode;
        this.tabStringMode = TabStringMode.COMBINED;
    }

    @Override
    public String apply(String specID, List<PeptideMatchData> peptideMatchData) {

        String txt = "";
        for (PeptideMatchData psm : peptideMatchData) {
            if (mode.equals("tab")) txt += getTabString(specID, psm);
            else if (mode.equals("txt")) txt += getTextString(psm);
            else if (mode.equals("pep")) txt += getPeptideString(psm);
            else if (mode.equals("fasta")) txt += getFastaString(psm);
            else return "";
        }

        return txt;
    }


    public String getHeader() {

        if (!mode.equals("tab")) return "";

        if (tabStringMode == TabStringMode.COMET)
            return "Spectrum\tCharge\tRT\tRank\tPeptide\tSequence\tProteins\tisDecoy\tComet.Xcorr\tComet.DeltaCn\tComet.SpScore\tComet.NegLogPv";
        if (tabStringMode == TabStringMode.MAXQUANT)
            return "Spectrum\tCharge\tRT\tRank\tPeptide\tSequence\tProteins\tisDecoy\tMaxQuant.DeltaMass(ppm)\tMaxQuant.Score\tMaxQuant.DeltaScore\tMaxQuant.Localisation.Prob";
        else
            return "Spectrum\tCharge\tRT\tComet.Rank\tComet.Peptide\tComet.Sequence\tComet.Proteins\tComet.isDecoy\tComet.Xcorr\tComet.DeltaCn\tComet.SpScore\tComet.NegLogPv\t" +
                    "MaxQuant.DeltaMass(ppm)\tMaxQuant.Score\tMaxQuant.DeltaScore\tMaxQuant.Localisation.Prob";
    }

    public String getTabString(String specID, PeptideMatchData psm) {

        if (tabStringMode == TabStringMode.COMET)
            return specID+"\t"+psm.getCharge()+"\t"+psm.getScore("rt")+"\t"+String.format("%.0f",psm.getScore("rank"))
                    +"\t"+getCometString(psm);
        if (tabStringMode == TabStringMode.MAXQUANT)
            return specID+"\t"+psm.getCharge()+"\t"+psm.getScore("rt")+"\t"+String.format("%.0f",psm.getScore("rank"))
                    +"\t"+getMaxQuantString(psm);
        else
            return specID+"\t"+psm.getCharge()+"\t"+psm.getScore("rt")+"\t"+String.format("%.0f",psm.getScore("rank"))
                    +"\t"+getCometCombString(psm)+"\t"+getMaxQuantCombString(psm);

    }

    protected String getCometCombString(PeptideMatchData psm) {

        return  psm.getPeptide().toString()+"\t"+psm.getPeptide().toSymbolString()+"\t"+psm.getProteinAcc()+"\t"+psm.isDecoy()+"\t"+
                psm.getScore("xcorr")+"\t"+psm.getScore("deltacn")+"\t"+psm.getScore("spscore")+"\t"+psm.getScore("neg_log10_p");
    }

    protected String getMaxQuantCombString(PeptideMatchData psm) {

        return psm.getScore("mass-error-ppm")+"\t"+psm.getScore("score")+"\t"+psm.getScore("delta-score")+"\t"+psm.getScore("mod-pos-prob");
    }

    protected String getCometString(PeptideMatchData psm) {

        return  psm.getPeptide().toString()+"\t"+psm.getPeptide().toSymbolString()+"\t"+psm.getProteinAcc()+"\t"+psm.isDecoy()+"\t"+
                psm.getScore("xcorr")+"\t"+psm.getScore("deltacn")+"\t"+psm.getScore("spscore")+"\t"+psm.getScore("neg_log10_p");
    }

    protected String getMaxQuantString(PeptideMatchData psm) {

        return  psm.getPeptide().toString()+"\t"+psm.getPeptide().toSymbolString()+"\t"+psm.getProteinAcc()+"\t"+psm.isDecoy()+"\t"+
                psm.getScore("mass-error-ppm")+"\t"+psm.getScore("score")+"\t"+psm.getScore("delta-score")+"\t"+psm.getScore("mod-pos-prob");
    }

    protected String getFastaString(PeptideMatchData psm) {

        return ">"+psm.getPeptide().toSymbolString()+";"+psm.getProteinAcc()+"\n"+psm.toSymbolString();
    }

    protected String getTextString(PeptideMatchData psm) {

        return psm.getPeptide().toString()+"\t"+psm.getProteinAcc();
    }

    protected String getPeptideString(PeptideMatchData psm) {

        return psm.getPeptide().toSymbolString();
    }

}
