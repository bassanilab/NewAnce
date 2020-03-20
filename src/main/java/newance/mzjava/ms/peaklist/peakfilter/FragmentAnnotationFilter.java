package newance.mzjava.ms.peaklist.peakfilter;

import com.google.common.base.Predicate;
import newance.mzjava.ms.peaklist.AbstractPeakProcessor;
import newance.mzjava.ms.peaklist.PeakAnnotation;

import java.util.List;

import static com.google.common.base.Preconditions.checkNotNull;

/**
 * Retain annotated peaks which the given predicate is true.
 *
 * @author fnikitin
 * Date: 6/11/13
 */
public class FragmentAnnotationFilter<A extends PeakAnnotation> extends AbstractPeakProcessor<A, A> {

    private final Predicate<A> predicate;

    public FragmentAnnotationFilter(Predicate<A> predicate) {

        checkNotNull(predicate);

        this.predicate = predicate;
    }

    @Override
    public void processPeak(double mz, double intensity, List<A> annotations) {

        for (A annotation : annotations) {

            if (predicate.apply(annotation)) {

                sink.processPeak(mz, intensity, annotations);
                break;
            }
        }
    }
}
