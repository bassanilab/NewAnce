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
import newance.psmconverter.PeptideSpectrumMatch;

import java.util.List;
import java.util.Map;
import java.util.function.BiFunction;

/**
 * @author Markus Müller
 */

public class CometPsm2StringFunction extends Psm2StringFunction {

    protected final GroupedFDRCalculator groupedFDRCalculator;
    protected final Map<String, Float> grpThresholdMap;


    public CometPsm2StringFunction(GroupedFDRCalculator groupedFDRCalculator, Map<String, Float> grpThresholdMap) {
        this.groupedFDRCalculator = groupedFDRCalculator;
        this.grpThresholdMap = grpThresholdMap;
    }

    @Override
    public String getScoreString(PeptideSpectrumMatch psm) {

        String lfdrStr = "NA";
        String pass = "NA";
        if (groupedFDRCalculator != null) {
            float lfdr = groupedFDRCalculator.getLocalFDR(psm);
            lfdrStr = String.format("%.5f",lfdr);

            if (grpThresholdMap!=null) pass = (lfdr<=grpThresholdMap.get(psm.getGroup()))?"true":"false";
        }
        String expectStr = String.format("%.5f",psm.getScore("expect"));
        String mdStr = String.format("%.5f",psm.getScore("mass_diff"));

        return  psm.getScore("xcorr")+"\t"+psm.getScore("deltacn")+"\t"+
                psm.getScore("spscore")+"\t"+expectStr+"\t+"+mdStr+"\t"+
                (int)psm.getScore("tot_num_ions")+"\t"+(int)psm.getScore("matched_num_ions")+"\t"+
                lfdrStr+"\t"+pass;
    }

    @Override
    public String getScoreHeader() {
        return "XCorr\tDeltaCn\tSpScore\tExpect\tMassdiff\tTot_num_ions\tNum_matched_ions\tComet.lFDR\tComet.passFDR";
    }

}
