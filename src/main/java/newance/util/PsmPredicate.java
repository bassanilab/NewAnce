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
import newance.psmconverter.PeptideMatchDataWrapper;

import java.io.Serializable;

/**
 * @author Markus Müller
 */

public class PsmPredicate implements Serializable {

    public enum ScoreOrder {LARGER,SMALLER};

    private final int minCharge;
    private final int maxCharge;
    private final int minPeptideLength;
    private final int maxPeptideLength;
    private final String scoreName;
    private final double minPsmScore;
    private final int maxRank;
    private final ScoreOrder order;

    public PsmPredicate(int minCharge, int maxCharge,
                        int minPeptideLength, int maxPeptideLength, int maxRank,
                        String scoreName, double minPsmScore, ScoreOrder order) {
        Preconditions.checkArgument(minPeptideLength <= maxPeptideLength);
        Preconditions.checkArgument(minCharge <= maxCharge);
        Preconditions.checkArgument(maxRank >= 0);

        this.minCharge = minCharge;
        this.maxCharge = maxCharge;
        this.minPeptideLength = minPeptideLength;
        this.maxPeptideLength = maxPeptideLength;
        this.scoreName = scoreName;
        this.minPsmScore = minPsmScore;
        this.maxRank = maxRank;
        this.order = order;
    }

    public boolean check(PeptideMatchDataWrapper psm, int charge) {

        if (!psm.getScoreMap().containsKey(this.scoreName)) return false;

        if (order==ScoreOrder.LARGER) {
            if (psm.getScore(scoreName) < minPsmScore) return false;
        } else {
            if (psm.getScore(scoreName) > minPsmScore) return false;
        }

        // !!!!temporary fix, need to be included in command line options
        if (psm.getScoreMap().containsKey("spscore")) {
            if (psm.getScore("spscore") < 100.0) return false;
        }

        if (psm.getRank() > this.maxRank) return false;

        String peptide = psm.getSequence();
        if (peptide.length() < this.minPeptideLength) {
            return false;
        } else if (peptide.length() > this.maxPeptideLength) {
            return false;
        }

        if (charge<minCharge) return false;
        if (charge>maxCharge) return false;

        return true;
    }

}
