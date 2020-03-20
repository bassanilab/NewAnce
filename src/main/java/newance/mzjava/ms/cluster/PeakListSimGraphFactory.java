package newance.mzjava.ms.cluster;

import newance.mzjava.mol.modification.Tolerance;
import newance.mzjava.ms.peaklist.Peak;
import newance.mzjava.ms.peaklist.PeakAnnotation;
import newance.mzjava.ms.peaklist.PeakList;
import newance.mzjava.ms.spectrasim.SimFunc;

import java.util.List;

import static com.google.common.base.Preconditions.checkState;

/**
 * A factory for building similarity graphs given a list of PeakLists that are sorted by charge and m/z.
 *
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class PeakListSimGraphFactory {

    private PeakListSimGraphFactory() {
    }

    /**
     * Build a new DenseSimilarityGraph
     *
     * @param spectra            the spectra, sorted by charge and m/z
     * @param simFunc            the sim function to use
     * @param precursorTolerance the precursor tolerance
     * @param scoreThreshold     if the similarity score is >= scoreThreshold a edge is added
     * @param <S>                the type of spectra
     * @return the new graph
     */
    public static <A extends PeakAnnotation, S extends PeakList<A>> SimilarityGraph<S> build(final List<S> spectra, SimFunc<A, A> simFunc,
                                                                                             final Tolerance precursorTolerance, final double scoreThreshold) {

        return build(spectra, new DenseSimilarityGraph.Builder<S>(), simFunc, precursorTolerance, scoreThreshold);
    }

    /**
     * Build a new SimilarityGraph
     *
     * @param spectra            the spectra, sorted by charge and m/z
     * @param builder            the builder to use
     * @param simFunc            the sim function to use
     * @param precursorTolerance the precursor tolerance
     * @param scoreThreshold     if the similarity score is >= scoreThreshold a edge is added
     * @param <S>                the type of spectra
     * @return the new graph
     */
    public static <A extends PeakAnnotation, S extends PeakList<A>, G extends SimilarityGraph<S>> G build(final List<S> spectra, final SimilarityGraphBuilder<S, G> builder,
                                                                                             SimFunc<A, A> simFunc, final Tolerance precursorTolerance, final double scoreThreshold) {

        for (int i = 0; i < spectra.size(); i++) {

            S spectrum1 = spectra.get(i);
            builder.add(spectrum1);
            Peak precursor1 = spectrum1.getPrecursor();

            for (int j = i + 1; j < spectra.size(); j++) {

                S spectrum2 = spectra.get(j);
                Peak precursor2 = spectrum2.getPrecursor();

                if(precursor2.compareTo(precursor1) < 0)
                    throw new IllegalStateException("spectra not sorted");

                if (precursor2.getCharge() > precursor1.getCharge() || precursorTolerance.check(precursor1.getMz(), precursor2.getMz()) == Tolerance.Location.LARGER)
                    break;

                checkState(precursorTolerance.withinTolerance(precursor1.getMz(), precursor2.getMz()));

                double score = simFunc.calcSimilarity(spectrum1, spectrum2);
                if (!Double.isNaN(score) && score >= scoreThreshold)
                    builder.add(spectrum1, spectrum2, score);
            }
        }

        return builder.build();
    }
}
