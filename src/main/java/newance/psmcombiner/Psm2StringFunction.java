/*
Copyright (C) SIB - Swiss Institute of Bioinformatics, Lausanne, Switzerland
Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
*/

package newance.psmcombiner;

import newance.mzjava.mol.Peptide;
import newance.mzjava.mol.modification.ModAttachment;
import newance.mzjava.mol.modification.Modification;
import newance.proteinmatch.SequenceVariant;
import newance.psmconverter.PeptideSpectrumMatch;

import java.util.List;
import java.util.function.BiFunction;

/**
 * @author Markus Müller
 */

public abstract class Psm2StringFunction implements BiFunction<String, List<PeptideSpectrumMatch>, String> {

    protected abstract String getScoreString(PeptideSpectrumMatch psm);
    protected abstract String getScoreHeader();

    @Override
    public String apply(String specID, List<PeptideSpectrumMatch> peptideSpectrumMatchData) {

        String txt = "";
        for (PeptideSpectrumMatch psm : peptideSpectrumMatchData) {
            if (!txt.isEmpty()) txt += "\n";
            txt += getSpecIDString(specID, psm) + "\t" + getPSMString(psm) + "\t" + getScoreString(psm);
        }

        return txt;
    }

    public String getHeader() {
        return getSpecIDHeader()+"\t"+getPSMHeader()+"\t"+getScoreHeader();
    }

    protected String getSpecIDHeader() {

        return "Spectrum\tScanNr\tCharge\tRT\tNeutralMass";
    }

    protected String getPSMHeader() {
        return "Peptide\tSequence\tPeptideMass\tModifName\tModifPosition\tModifMass\tModifAA\tProteins\tIsVariant\t" +
                "VariantStart\tVariantEnd\tVariantSeq\tVariantID\tWTSeq\tIsDecoy\tRank\tGroup";
    }

    protected String getSpecIDString(String specID, PeptideSpectrumMatch psm) {

        String rt = String.format("%.5f",psm.getRetentionTime());
        String mass = String.format("%.5f",psm.getNeutralPrecMass());

        return specID+"\t"+psm.getScanNr()+"\t"+psm.getCharge()+"\t"+rt+"\t"+mass;

    }

    protected String getPSMString(PeptideSpectrumMatch psm) {

        Peptide peptide = psm.getPeptide();
        String protACs = psm.getProteinIDs().toString();
        int rank = (int) psm.getScore("rank");
        String pepMass = String.format("%.5f",peptide.getMolecularMass());
        String modifString = getModifString(peptide);
        String variantString = getVariantString(psm);

        return  peptide.toString()+"\t"+peptide.toSymbolString()+"\t"+pepMass+"\t"+modifString+"\t"+protACs+"\t"+
                variantString+"\t"+psm.isDecoy()+"\t"+ rank+"\t"+ psm.getGroup();

    }

    private String getModifString(Peptide peptide)  {

        if (!peptide.hasModifications()) return "NA\tNA\tNA\tNA";

        String modifNames = "";
        String modifPos = "";
        String modifMass = "";
        String modifAA = "";

        if (peptide.hasModificationAt(ModAttachment.N_TERM)) {
            for (Modification modif : peptide.getModifications(ModAttachment.nTermSet)) {
                modifNames += (modifNames.isEmpty())?modif.getLabel():","+modif.getLabel();
                modifPos += (modifPos.isEmpty())?"0":",0";
                String massStr = String.format("%.5f",modif.getMolecularMass());
                modifMass += (modifMass.isEmpty())?massStr:","+massStr;
                modifAA += (modifAA.isEmpty())?"NT":",NT";
            }
        }

        if (peptide.hasModificationAt(ModAttachment.SIDE_CHAIN)) {
            for (int i : peptide.getModificationIndexes(ModAttachment.sideChainSet)) {

                for (Modification modif : peptide.getModificationsAt(i, ModAttachment.sideChainSet)) {
                    modifNames += (modifNames.isEmpty())?modif.getLabel():","+modif.getLabel();
                    String posStr = String.format("%d",i+1);
                    modifPos += (modifPos.isEmpty())?posStr:","+posStr;
                    String massStr = String.format("%.5f",modif.getMolecularMass());
                    modifMass += (modifMass.isEmpty())?massStr:","+massStr;
                    modifAA += (modifAA.isEmpty())?peptide.getSymbol(i).getSymbol():","+peptide.getSymbol(i).getSymbol();
                }
            }
        }

        if (peptide.hasModificationAt(ModAttachment.C_TERM)) {
            for (Modification modif : peptide.getModifications(ModAttachment.cTermSet)) {
                modifNames += (modifNames.isEmpty())?modif.getLabel():","+modif.getLabel();
                String posStr = String.format("%d",peptide.size()+1);
                modifPos += (modifPos.isEmpty())?posStr:","+posStr;
                String massStr = String.format("%.5f",modif.getMolecularMass());
                modifMass += (modifMass.isEmpty())?massStr:","+massStr;
                modifAA += (modifAA.isEmpty())?"CT":",CT";
            }
        }

        return modifNames+"\t"+modifPos+"\t"+modifMass+"\t"+modifAA;
    }

    public String getVariantString(PeptideSpectrumMatch psm) {

        String wtSeqStr = psm.getWtSequence();
        if (!psm.isVariant()) return "false\tNA\tNA\tNA\tNA\t"+wtSeqStr;

        String posStartStr = "";
        String posEndStr = "";
        String mutSeqStr = "";
        String annotStr = "";

        for (SequenceVariant variant : psm.getVariants()) {

            posStartStr += (posStartStr.isEmpty())?(variant.getStartWT()+1):","+(variant.getStartWT()+1);
            posEndStr += (posEndStr.isEmpty())?(variant.getEndWT()+1):","+(variant.getEndWT()+1);
            mutSeqStr += (mutSeqStr.isEmpty())?variant.getMutatedSequence():","+variant.getMutatedSequence();
            annotStr += (annotStr.isEmpty())?variant.getInfo():","+variant.getInfo();
        }

        return "true\t"+posStartStr+"\t"+posEndStr+"\t"+mutSeqStr+"\t"+annotStr+"\t"+wtSeqStr;
    }

}
