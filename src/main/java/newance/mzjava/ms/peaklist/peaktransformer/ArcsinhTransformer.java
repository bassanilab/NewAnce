package newance.mzjava.ms.peaklist.peaktransformer;

import newance.mzjava.ms.peaklist.AbstractPeakProcessor;
import newance.mzjava.ms.peaklist.PeakAnnotation;

import java.util.List;

/**
 * Following Huber et al. Bioinformatics, 18, S96-S104, 2003 we apply the arcsinh transformation t,
 * t: intensity -> gamma*arcsinh(a + b*intensity), where arcsinh(x) = log(x + sqrt(x^2 + 1)),
 * which stabilizes the variance of the peak intensities, i.e. after this transformation all peak intensities
 * have the same variance independent of their height. This transformation is equal to the log transformation
 * for large intensities. According to Huber et al, this transformation is valid for error models of the form:
 * var(I) = (c1*I+c2)^2 + c3 i.e. the variance of the intensity I depend quadratically on the intensity itself.
 * The parameters gamma, a and b are related to the ci's as follows: gamma = 1/c1, a=c2/sqrt(c3), b=c1/sqrt(c3),
 * or may be estimated otherwise.
 *
 * Created with IntelliJ IDEA.
 * User: Markus
 * Date: 28.06.12
 * Time: 18:06
 * To change this template use File | Settings | File Templates.
 */
public class ArcsinhTransformer<A extends PeakAnnotation> extends AbstractPeakProcessor<A, A> {
    private final double a;
    private final double b;
    private final double gamma;

    /**
     * Constructor
     * @param gamma gamma in gamma*arcsinh(a + b*intensity)
     * @param a a in gamma*arcsinh(a + b*intensity)
     * @param b b in gamma*arcsinh(a + b*intensity)
     */
    public ArcsinhTransformer(double gamma, double a, double b)
    {
        this.a = a;
        this.b = b;
        this.gamma = gamma;
    }

    @Override
    public void processPeak(double mz, double intensity, List<A> annotations) {
        intensity = a+b*intensity;
        intensity = gamma*Math.log(intensity+Math.sqrt(intensity*intensity+1.0));
        sink.processPeak(mz, intensity, annotations);
    }
}
