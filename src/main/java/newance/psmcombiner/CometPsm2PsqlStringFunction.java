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

import java.util.Arrays;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * @author Emma Ricart
 */

public class CometPsm2PsqlStringFunction extends Psm2StringFunction {

    protected final GroupedFDRCalculator groupedFDRCalculator;
    protected final Map<String, Float> grpThresholdMap;


    public CometPsm2PsqlStringFunction(GroupedFDRCalculator groupedFDRCalculator, Map<String, Float> grpThresholdMap) {
        this.groupedFDRCalculator = groupedFDRCalculator;
        this.grpThresholdMap = grpThresholdMap;
    }

    @Override
    public String getScoreString(PeptideSpectrumMatch psm) {

        float lfdr = groupedFDRCalculator.getLocalFDR(psm);
        String lfdrStr = (groupedFDRCalculator ==null)?"":String.format("%.5f",lfdr);
        String expectStr = String.format("%.5f",psm.getScore("expect"));
        String mdStr = String.format("%.5f",psm.getScore("mass_diff"));
        String pass = "NA";
        if (grpThresholdMap!=null) pass = (lfdr<=grpThresholdMap.get(psm.getGroup()))?"true":"false";

        return  psm.getScore("xcorr")+"\t"+psm.getScore("deltacn")+"\t"+
                psm.getScore("spscore")+"\t"+expectStr+"\t+"+mdStr+"\t"+
                (int)psm.getScore("tot_num_ions")+"\t"+(int)psm.getScore("matched_num_ions")+"\t"+
                lfdrStr+"\t"+pass;
    }

    @Override
    public String getScoreHeader() {
        return "XCorr\tDeltaCn\tSpScore\tExpect\tMassdiff\tTot_num_ions\tNum_matched_ions\tComet.lFDR\tComet.passFDR";
    }

    @Override
    protected String getPSMString(PeptideSpectrumMatch psm) {

        Peptide peptide = psm.getPeptide();
        //edited by Emma
        String protACs ="{"+ psm.getProteinIDs().stream()
                .map(proteinID -> "\"" + proteinID + "\"" )
                .collect(Collectors.joining(","))+"}";
        int rank = (int) psm.getScore("rank");
        String pepMass = String.format("%.5f",peptide.getMolecularMass());
        String modifString = getModifString(peptide);
        String variantString = getVariantString(psm);

        return  peptide.toString()+"\t"+peptide.toSymbolString()+"\t"+pepMass+"\t"+modifString+"\t"+protACs+"\t"+
                variantString+"\t"+psm.isDecoy()+"\t"+ rank+"\t"+ psm.getGroup();

    }

    @Override
    protected String getModifString(Peptide peptide)  {

        if (!peptide.hasModifications()) return "NA\tNA\tNA\tNA";

        int numMods=peptide.getModificationCount();

        String[] a_modifNames = new String[numMods];
        int[] a_modifPos = new int [numMods];
        double[] a_modifMass = new double[numMods];
        String[] a_modifAA = new String[numMods];

        int n = 0;


        if (peptide.hasModificationAt(ModAttachment.N_TERM)) {
            for (Modification modif : peptide.getModifications(ModAttachment.nTermSet)) {
                a_modifNames[n] = modif.getLabel();
                a_modifPos[n] = 0;
                a_modifMass[n] = modif.getMolecularMass();
                a_modifAA[n] ="NT";
                n++;
            }
        }

        if (peptide.hasModificationAt(ModAttachment.SIDE_CHAIN)) {
            for (int i : peptide.getModificationIndexes(ModAttachment.sideChainSet)) {
                for (Modification modif : peptide.getModificationsAt(i, ModAttachment.sideChainSet)) {
                    a_modifNames[n] = modif.getLabel();
                    a_modifPos[n] = i+1;
                    a_modifMass[n] = modif.getMolecularMass();
                    a_modifAA[n] =peptide.getSymbol(i).getSymbol();
                    n++;
                }
            }
        }

        if (peptide.hasModificationAt(ModAttachment.C_TERM)) {
            for (Modification modif : peptide.getModifications(ModAttachment.cTermSet)) {

                a_modifNames[n] = modif.getLabel();
                a_modifPos[n] = peptide.size()+1;
                a_modifMass[n] = modif.getMolecularMass();
                a_modifAA[n] ="CT";
                n++;
            }
        }

        String modifNames =strArray2PsqlString(a_modifNames);
        String modifPos = intArray2PsqlString(a_modifPos);
        String modifMass = "{"+ Arrays.stream(a_modifMass)
                .mapToObj(d -> String.format("%.5f",d))
                .collect(Collectors.joining(","))+"}";
        String modifAA =strArray2PsqlString(a_modifAA);

        return modifNames+"\t"+modifPos+"\t"+modifMass+"\t"+modifAA;
    }
    @Override
    public String getVariantString(PeptideSpectrumMatch psm) {

        String wtSeqStr = psm.getWtSequence();
        if (!psm.isVariant()) return "false\tNA\tNA\tNA\tNA\t"+wtSeqStr;

        int numVars=psm.getVariants().size();

        int[] a_posStartStr = new int[numVars];
        int[] a_posEndStr = new int[numVars];
        String[] a_mutSeqStr = new String[numVars];
        String[] a_annotStr = new String[numVars];

        int n=0;

        for (SequenceVariant variant : psm.getVariants()) {
            a_posStartStr[n] = variant.getStartWT()+1;
            a_posEndStr[n] = variant.getEndWT()+1;
            a_mutSeqStr[n] = variant.getMutatedSequence();
            a_annotStr[n] = variant.getInfo();
            n++;
        }

        String posStartStr = intArray2PsqlString(a_posStartStr);
        String posEndStr = intArray2PsqlString(a_posEndStr);
        String mutSeqStr =strArray2PsqlString(a_mutSeqStr);
        String annotStr = strArray2PsqlString(a_annotStr);

        return "true\t"+posStartStr+"\t"+posEndStr+"\t"+mutSeqStr+"\t"+annotStr+"\t"+wtSeqStr;
    }


    private static String intArray2PsqlString (int[] array) {
        return "{"+ Arrays.stream(array)
                .mapToObj(i -> String.format("%d",i))
                .collect(Collectors.joining(","))+"}";

    }

    private static String strArray2PsqlString (String[] array) {
        return "{"+ Arrays.stream(array)
                .map(s -> "\"" + s + "\"")
                .collect(Collectors.joining(","))+"}";
    }


}
