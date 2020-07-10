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

import newance.psmconverter.PeptideSpectrumMatch;

import java.util.Map;

/**
 * @author Markus Müller
 */

public class MaxQuantPsm2StringFunction extends Psm2StringFunction {

    protected final GroupedFDRCalculator groupedFDRCalculator;
    protected final Map<String, Float> grpThresholdMap;


    public MaxQuantPsm2StringFunction(GroupedFDRCalculator groupedFDRCalculator, Map<String, Float> grpThresholdMap) {
        this.groupedFDRCalculator = groupedFDRCalculator;
        this.grpThresholdMap = grpThresholdMap;
    }

    @Override
    public String getScoreString(PeptideSpectrumMatch psm) {

        float lfdr = groupedFDRCalculator.getLocalFDR(psm);
        String lfdrStr = (groupedFDRCalculator ==null)?"":String.format("%.5f",lfdr);
        String pass = "NA";
        if (grpThresholdMap!=null) pass = (lfdr<=grpThresholdMap.get(psm.getGroup()))?"true":"false";

        return psm.getScore("Mass Error [ppm]")+"\t"+psm.getScore("Score")+"\t"+
                psm.getScore("Delta score")+"\t"+psm.getScore("PEP")+"\t"+
                psm.getScore("Localization prob")+"\t"+lfdrStr+"\t"+psm.getGroup()+"\t"+pass;
    }

    @Override
    public String getScoreHeader() {
        return "Mass.Error[ppm]\tScore\tDelta.score\tPEP\tLocalization.prob\tlFDR\tGroup\tpassFDR";
    }

}
