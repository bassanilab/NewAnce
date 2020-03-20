package newance.mzjava.ms.peaklist.peaktransformer;

import newance.mzjava.ms.peaklist.AbstractPeakProcessor;
import newance.mzjava.ms.peaklist.PeakAnnotation;

import java.util.List;

/**
 * Performs a square root transform on peak intensities.
 * <p/>
 * Created with IntelliJ IDEA.
 * @author Markus
 */
public class SqrtTransformer<A extends PeakAnnotation> extends AbstractPeakProcessor<A, A> {

    @Override
    public void processPeak(double mz, double intensity, List<A> annotations) {

        double transformed = intensity <= 0 ? 0 : Math.sqrt(intensity);

        sink.processPeak(mz, transformed, annotations);
    }
}
