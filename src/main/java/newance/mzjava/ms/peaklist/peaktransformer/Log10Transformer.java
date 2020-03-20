package newance.mzjava.ms.peaklist.peaktransformer;


import newance.mzjava.ms.peaklist.AbstractPeakProcessor;
import newance.mzjava.ms.peaklist.PeakAnnotation;

import java.util.List;

/**
 * Performs a base 10 logarithm transform on peak intensities + 1.
 * <p/>
 * @author Markus
 */
public class Log10Transformer<A extends PeakAnnotation> extends AbstractPeakProcessor<A, A> {

    @Override
    public void processPeak(double mz, double intensity, List<A> annotations) {

        double transformed = intensity <= 0 ? 0 : Math.log10(intensity + 1);

        sink.processPeak(mz, transformed, annotations);
    }
}
