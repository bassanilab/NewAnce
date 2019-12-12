package newance.util;

import com.google.common.base.Preconditions;
import org.expasy.mzjava.proteomics.ms.ident.PeptideMatch;

import java.io.Serializable;

/**
 * Copyright (C) 2019
 * @author Markus Müller
 * @Institutions: SIB, Swiss Institute of Bioinformatics; Ludwig Institute for Cancer Research
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

    public boolean check(PeptideMatch psm, int charge) {

        if (!psm.getScoreMap().containsKey(this.scoreName)) return false;

        if (order==ScoreOrder.LARGER) {
            if (psm.getScore(scoreName) < minPsmScore) return false;
        } else {
            if (psm.getScore(scoreName) > minPsmScore) return false;
        }

        if (psm.getScore("rank") > this.maxRank) return false;

        String peptide = psm.toSymbolString();
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
