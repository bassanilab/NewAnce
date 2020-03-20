package newance.mzjava.ms.peaklist.peaktransformer;

import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TIntObjectHashMap;
import newance.mzjava.ms.peaklist.DelayedPeakProcessor;
import newance.mzjava.ms.peaklist.PeakAnnotation;
import org.apache.commons.math3.exception.MathInternalError;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.stat.ranking.TiesStrategy;
import org.apache.commons.math3.util.FastMath;

import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * Class to rank apply intensities. The highest peak gets rank 1, the runner up 1 - 1/N, where N is the total
 * number of peaks, the third highest peak 1 - 2/N, and the lowest peak is assigned a value of 1/N. The algorithm
 * to calculate the ranks is copied from org.apache.commons.math.stat.ranking.NaturalRanking with slight adaptations.
 *
 * @author  Markus Muller
 * @version 1.0
 */
public class RankNormalizer<A extends PeakAnnotation>  extends DelayedPeakProcessor<A, A> {

    /**
     * Strategy to resolve ties.
     */
    private final TiesStrategy tiesStrategy;
    /** Source of random data - used only when ties strategy is RANDOM */
    private final RandomDataGenerator randomData;

    private final Comparator<double[]> comparator = new Comparator<double[]>() {
        @Override
        public int compare(double[] o1, double[] o2) {

            return Double.compare(o1[1], o2[1]);
        }
    };

    public RankNormalizer() {

        tiesStrategy = TiesStrategy.AVERAGE;
        randomData = null;
    }

    public RankNormalizer(TiesStrategy tiesStrategy) {

        this.tiesStrategy = tiesStrategy;

        randomData = tiesStrategy == TiesStrategy.RANDOM ? new RandomDataGenerator() : null;
    }

    @Override
    protected void processCached(TDoubleArrayList mzList, TDoubleArrayList intensityList, TIntObjectHashMap<List<A>> annotationMap) {

        TDoubleArrayList peakRanks = rank(intensityList);
        int size =  mzList.size();
        for(int i = 0; i <size; i++) {

            List<A> annotations;
            if(annotationMap.contains(i))
                annotations = annotationMap.get(i);
            else
                annotations = Collections.emptyList();

            sink.processPeak(mzList.get(i), peakRanks.get(i) / size, annotations);
        }
    }

    double[][] ranks = new double[0][2];

    //borrowed from org.apache.commons.math.stat.ranking.NaturalRanking
    TDoubleArrayList rank(TDoubleArrayList data) {

        if(ranks.length < data.size()) ranks = new double[data.size()][2];

        for (int i = 0; i < data.size(); i++) {

            ranks[i][0] = i;
            ranks[i][1] = data.get(i);
        }

        // Sort the IntDoublePairs
        Arrays.sort(ranks, 0, data.size(), comparator);

        // Walk the sorted array, filling output array using sorted positions,
        // resolving ties as we go
        int pos = 1;  // position in sorted array
        data.set((int) ranks[0][0], pos);
        TIntArrayList tiesTrace = new TIntArrayList();
        tiesTrace.add((int)ranks[0][0]);
        for (int i = 1; i < data.size(); i++) {

            if (Double.compare(ranks[i][1], ranks[i - 1][1]) > 0) {
                // tie sequence has ended (or had length 1)
                pos = i + 1;
                if (tiesTrace.size() > 1) {  // if seq is nontrivial, resolve
                    resolveTie(data, tiesTrace);
                }
                tiesTrace.resetQuick();
                tiesTrace.add((int)ranks[i][0]);
            } else {
                // tie sequence continues
                tiesTrace.add((int)ranks[i][0]);
            }
            int index = (int) ranks[i][0];

            data.set(index, pos);
        }
        if (tiesTrace.size() > 1) {  // handle tie sequence at end
            resolveTie(data, tiesTrace);
        }
        return data;
    }

    //borrowed from org.apache.commons.math.stat.ranking.NaturalRanking
    private void resolveTie(TDoubleArrayList ranks, TIntArrayList tiesTrace) {

        // constant value of ranks over tiesTrace
        final double c = ranks.get(tiesTrace.get(0));

        // length of sequence of tied ranks
        final int length = tiesTrace.size();

        switch (tiesStrategy) {
            case  AVERAGE:  // Replace ranks with average
                fill(ranks, tiesTrace, (2 * c + length - 1) / 2d);
                break;
            case MAXIMUM:   // Replace ranks with maximum values
                fill(ranks, tiesTrace, c + length - 1);
                break;
            case MINIMUM:   // Replace ties with minimum
                fill(ranks, tiesTrace, c);
                break;
            case RANDOM:    // Fill with random integral values in [c, c + length - 1]

                long f = FastMath.round(c);
                for (int i = 0; i < tiesTrace.size(); i++) {

                    ranks.set(tiesTrace.get(i), randomData.nextLong(f, f + length - 1));
                }
                break;
            case SEQUENTIAL:  // Fill sequentially from c to c + length - 1
                // walk and fill

                f = FastMath.round(c);
                int i = 0;
                for (int j = 0; j < tiesTrace.size(); j++) {
                    ranks.set(tiesTrace.get(j), f + i++);
                }
                break;
            default: // this should not happen unless TiesStrategy enum is changed
                throw new MathInternalError();
        }
    }

    //borrowed from org.apache.commons.math.stat.ranking.NaturalRanking
    private void fill(TDoubleArrayList data, TIntArrayList tiesTrace, double value) {

        for (int i = 0; i < tiesTrace.size(); i++) {

            data.set(tiesTrace.get(i), value);
        }
    }
}
