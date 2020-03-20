package newance.mzjava.ms.peaklist.peaktransformer;


import newance.mzjava.ms.peaklist.AbstractPeakProcessor;
import newance.mzjava.ms.peaklist.PeakAnnotation;

import java.util.List;

/**
 * Following Huber et al. Bioinformatics, 18, S96-S104, 2003 the arcsinh transformation t,
 * t: intensity -> gamma*arcsinh(a + b*intensity), where arcsinh(x) = log(x + sqrt(x^2 + 1)),
 * stabilizes the variance of the peak intensities, i.e. after this transformation all peak intensities
 * have the same variance independent of their height. This transformation is equal to the log transformation
 * for large intensities. According to Huber et al, this is valid for error models of the form:
 * var(I) = (c1*I+c2)^2 + c3 i.e. the variance of the intensity I depend quadratically on the intensity itself.
 * The parameters gamma, a and b are related to the ci's as follows: gamma = 1/c1, a=c2/sqrt(c3), b=c1/sqrt(c3),
 * or may be estimated otherwise.
 *
 * This class inverses the arcsinh transformation
 *
 * @author Markus Muller
 * @version 0.0
 */
public class InverseArcsinhTransformer<A extends PeakAnnotation> extends AbstractPeakProcessor<A, A> {
        private final double a;
        private final double b;
        private final double gamma;

        /**
         * Constructor
         * @param gamma gamma in gamma*arcsinh(a + b*intensity)
         * @param a a in gamma*arcsinh(a + b*intensity)
         * @param b b in gamma*arcsinh(a + b*intensity)
         */
        public InverseArcsinhTransformer(double gamma, double a, double b)
        {
            this.a = a;
            this.b = b;
            this.gamma = gamma;
        }

        @Override
        public void processPeak(double mz, double intensity, List<A> annotations) {
            intensity = intensity/gamma;
            intensity = (Math.sinh(intensity)-a)/b;
            sink.processPeak(mz, intensity, annotations);
        }
    }
