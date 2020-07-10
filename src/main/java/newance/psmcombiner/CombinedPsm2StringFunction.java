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

/**
 * @author Markus Müller
 */

public class CombinedPsm2StringFunction extends Psm2StringFunction {

    protected final GroupedFDRCalculator cometGroupedFDRCalculator;
    protected final GroupedFDRCalculator maxquantGroupedFDRCalculator;
    protected final Map<String, Float> cometGrpThresholdMap;
    protected final Map<String, Float> maxquantGrpThresholdMap;

    protected final CometPsm2StringFunction cometPsm2StringFunction;
    protected final MaxQuantPsm2StringFunction maxQuantPsm2StringFunction;


    public CombinedPsm2StringFunction(GroupedFDRCalculator cometGroupedFDRCalculator,
                                      GroupedFDRCalculator maxquantGroupedFDRCalculator,
                                      Map<String, Float> cometGrpThresholdMap,
                                      Map<String, Float> maxquantGrpThresholdMap) {

        this.cometGroupedFDRCalculator = cometGroupedFDRCalculator;
        this.maxquantGroupedFDRCalculator = maxquantGroupedFDRCalculator;
        this.cometGrpThresholdMap = cometGrpThresholdMap;
        this.maxquantGrpThresholdMap = maxquantGrpThresholdMap;

        cometPsm2StringFunction = new CometPsm2StringFunction(cometGroupedFDRCalculator, cometGrpThresholdMap);
        maxQuantPsm2StringFunction = new MaxQuantPsm2StringFunction(maxquantGroupedFDRCalculator, maxquantGrpThresholdMap);
    }

    @Override
    public String getScoreString(PeptideSpectrumMatch psm) {
        return cometPsm2StringFunction.getScoreString(psm)+"\t"+maxQuantPsm2StringFunction.getScoreString(psm);
    }

    @Override
    public String getScoreHeader() {
        return cometPsm2StringFunction.getScoreHeader()+"\t"+maxQuantPsm2StringFunction.getScoreHeader();
    }

}
