/*
Copyright (C) SIB - Swiss Institute of Bioinformatics, Lausanne, Switzerland
Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
*/

package newance.scripts;

import newance.mzjava.mol.Peptide;
import newance.psmcombiner.Psm2StringFunction;
import newance.psmconverter.PeptideSpectrumMatch;

import java.util.List;

/**
 * @author Markus Müller
 */

public class PepXMLPSm2StringFunction extends Psm2StringFunction {

    public PepXMLPSm2StringFunction() {
    }

    @Override
    public String getScoreString(PeptideSpectrumMatch psm) {

        String hyperscoreStr = String.format("%.5f",psm.getScore("hyperscore"));
        String nextscoreStr = String.format("%.5f",psm.getScore("nextscore"));
        String mdStr = String.format("%.5f",psm.getScore("mass_diff"));
        String expectStr = String.format("%.5f",psm.getScore("expect"));
        int tot_num_ions = (int)psm.getScore("tot_num_ions");
        int matched_num_ions = (int)psm.getScore("matched_num_ions");

        return  hyperscoreStr+"\t"+nextscoreStr+"\t"+expectStr+"\t"+mdStr+"\t"+tot_num_ions+"\t"+matched_num_ions;
    }

    @Override
    public String apply(String specID, List<PeptideSpectrumMatch> peptideSpectrumMatchData) {

        String[] fields  = specID.split("\\.");
        String specID_pepxml = String.format("%s.%05d.%05d.%s", fields[0], Integer.parseInt(fields[1]), Integer.parseInt(fields[2]), fields[3]);
        String txt = "";
        for (PeptideSpectrumMatch psm : peptideSpectrumMatchData) {
            if (!txt.isEmpty()) txt += "\n";
            txt += getSpecIDString(specID_pepxml, psm) + "\t" + getPSMString(psm) + "\t" + getScoreString(psm);
        }

        return txt;
    }

    @Override
    public String getPSMHeader() {
        return "Peptide\tSequence\tPeptideMass\tModifName\tModifPosition\tModifMass\tModifAA\tProteins\tRank";
    }

    @Override
    public String getScoreHeader() {
        return "hyperscore\tnextscore\texpect\tmass_diff\ttot_num_ions\tmatched_num_ions";
    }


    protected String getPSMString(PeptideSpectrumMatch psm) {

        Peptide peptide = psm.getPeptide();
        String protACs = psm.getProteinIDs().toString();
        int rank = (int) psm.getScore("rank");
        rank++;
        String pepMass = String.format("%.5f",peptide.getMolecularMass());
        String modifString = getModifString(peptide);

        return  peptide.toString()+"\t"+peptide.toSymbolString()+"\t"+pepMass+"\t"+modifString+"\t"+protACs+"\t"+rank;

    }

}
