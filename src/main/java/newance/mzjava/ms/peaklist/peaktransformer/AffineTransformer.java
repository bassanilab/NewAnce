package newance.mzjava.ms.peaklist.peaktransformer;


import newance.mzjava.ms.peaklist.AbstractPeakProcessor;
import newance.mzjava.ms.peaklist.PeakAnnotation;

import java.util.List;

/**
 * Performs affine transformation (slope*Intensity+offset) on peak intensities.
 * <p/>
 * @author Markus
 */
public class AffineTransformer<A extends PeakAnnotation> extends AbstractPeakProcessor<A, A> {

    private final double slope;
    private final double offset;

    /**
     * Make an instance of affine transformer (slope*Intensity+offset) on peak intensities.
     *
     * @param slope  slope of affine transformation
     * @param offset offset of affine transformation
     */
    public AffineTransformer(double slope, double offset) {

        this.slope = slope;
        this.offset = offset;
    }

    @Override
    public void processPeak(double mz, double intensity, List<A> annotations) {

        sink.processPeak(mz, slope * intensity + offset, annotations);
    }
}
