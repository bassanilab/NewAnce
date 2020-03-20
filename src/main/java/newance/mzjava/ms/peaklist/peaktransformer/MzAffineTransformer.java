package newance.mzjava.ms.peaklist.peaktransformer;

import newance.mzjava.ms.peaklist.AbstractPeakProcessor;
import newance.mzjava.ms.peaklist.PeakAnnotation;

import java.util.List;

/**
 * Performs affine transformation (slope*mz+offset) on peak m/z.
 * <p/>
 * @author Markus
 */
public class MzAffineTransformer<A extends PeakAnnotation> extends AbstractPeakProcessor<A, A> {

    private final double slope;
    private final double offset;

    /**
     * Make an instance of affine transformer (slope*mz+offset) on peak intensities.
     *
     * @param slope  slope of affine transformation
     * @param offset offset of affine transformation
     */
    public MzAffineTransformer(double slope, double offset) {

        this.slope = slope;
        this.offset = offset;
    }

    @Override
    public void processPeak(double mz, double intensity, List<A> annotations) {

        sink.processPeak(slope * mz + offset, intensity, annotations);
    }
}
