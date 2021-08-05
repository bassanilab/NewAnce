/*
Copyright (C) SIB - Swiss Institute of Bioinformatics, Lausanne, Switzerland
Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
*/

package newance.util;

import com.google.common.base.Preconditions;
import gnu.trove.map.TObjectDoubleMap;
import newance.psmconverter.PeptideMatchDataWrapper;

import java.io.Serializable;

/**
 * @author Markus Müller
 */

public class PsmPredicate implements Serializable {

    private final NewAnceParams params;

    public PsmPredicate(NewAnceParams params) {
        this.params = params;
    }

    public boolean check(PeptideMatchDataWrapper psm, int charge) {

        TObjectDoubleMap<String> scoreMap = psm.getScoreMap();

        if (psm.getSequence().equals("LSSSSQHGPSY")){
            System.out.println(psm.toPeptide());
        }
        if (psm.getRank() > params.getMaxRank()) return false;

        if (charge<params.getMinCharge()) return false;
        if (charge>params.getMaxCharge()) return false;

        if (scoreMap.containsKey("spscore") && scoreMap.get("spscore")<params.getMinSpScorePSM()) return false;
        if (scoreMap.containsKey("xcorr") && scoreMap.get("xcorr")<params.getMinXCorrPSM()) return false;
        if (scoreMap.containsKey("deltacn") && scoreMap.get("deltacn")<params.getMinDeltaCnPSM()) return false;
        if (scoreMap.containsKey("Score") && scoreMap.get("Score")<params.getMinScorePSM()) return false;
        if (scoreMap.containsKey("Delta score") && scoreMap.get("Delta score")<params.getMinDeltaScorePSM()) return false;
        if (scoreMap.containsKey("PEP") && scoreMap.get("PEP")>params.getMinPEPPSM()) return false;

        String peptide = psm.getSequence();
        if (peptide.length() < params.getMinPeptideLength()) {
            return false;
        } else if (peptide.length() > params.getMaxPeptideLength()) {
            return false;
        }

        return true;
    }

}
