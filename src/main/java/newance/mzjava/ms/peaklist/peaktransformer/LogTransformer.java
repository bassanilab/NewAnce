package newance.mzjava.ms.peaklist.peaktransformer;


import newance.mzjava.ms.peaklist.AbstractPeakProcessor;
import newance.mzjava.ms.peaklist.PeakAnnotation;

import java.util.List;

/**
 * Performs a base e (natural) logarithm transform on peak intensities + 1.
 *
 * @author Markus
 */
public class LogTransformer<A extends PeakAnnotation> extends AbstractPeakProcessor<A, A> {

    @Override
    public void processPeak(double mz, double intensity, List<A> annotations) {

        double transformed = intensity <= 0 ? 0 : Math.log(intensity + 1);

        sink.processPeak(mz, transformed, annotations);
    }
}
