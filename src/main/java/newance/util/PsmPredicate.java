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
