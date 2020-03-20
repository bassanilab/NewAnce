package newance.mzjava.ms.library;

import newance.mzjava.ms.peaklist.PeakAnnotation;
import newance.mzjava.ms.peaklist.PeakList;
import newance.mzjava.ms.spectrasim.SimFunc;

import java.util.Comparator;
import java.util.PriorityQueue;

/**
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class SimFuncProcedure<LA extends PeakAnnotation, LPL extends PeakList<LA>, QA extends PeakAnnotation, QPL extends PeakList<QA>> implements Procedure<LPL> {

    private final QPL queryPeakList;
    private final SimFunc<LA, QA> simFunc;

    private final PriorityQueue<Result> priorityQueue;

    public SimFuncProcedure(QPL queryPeakList, SimFunc<LA, QA> simFunc, int capacity) {

        this.queryPeakList = queryPeakList;
        this.simFunc = simFunc;
        priorityQueue = new PriorityQueue<>(capacity + 1, new Comparator<Result>() {
            @Override
            public int compare(Result r1, Result r2) {

                return Double.compare(r1.score, r2.score);
            }
        });
    }

    @Override
    public void execute(LPL libPeakList) {

        double score = simFunc.calcSimilarity(libPeakList, queryPeakList);


    }

    private static class Result<PL> {

        final double score;
        final PL peakList;

        public Result(double score, PL peakList) {

            this.score = score;
            this.peakList = peakList;
        }
    }
}
