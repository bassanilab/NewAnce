package newance.mzjava.ms.peaklist;

import gnu.trove.list.TDoubleList;
import gnu.trove.map.TIntObjectMap;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static com.google.common.base.Preconditions.checkNotNull;

/**
 * A list of peak processors.
 * <p>
 * <b>Note that this implementation is not synchronized.</b> If multiple threads access a PeakProcessorChain
 * instance concurrently the process method must be synchronized externally.
 *
 * @author Oliver Horlacher
 * @version 1.0
 */
public class PeakProcessorChain<A extends PeakAnnotation> {

    private final List<PeakProcessor<A, A>> list = new ArrayList<PeakProcessor<A, A>>();
    private PeakSink<A> sink;

    public PeakProcessorChain() {}

    public PeakProcessorChain(PeakProcessor<A, A> peakProcessor) {

        checkNotNull(peakProcessor);

        list.add(peakProcessor);
    }

    public PeakProcessorChain<A> add(PeakProcessor<A, A> peakProcessor) {

        checkNotNull(peakProcessor);

        list.add(peakProcessor);

        return this;
    }

    public void process(PeakList<A> peakList, PeakSink<A> peakSink) {

        checkNotNull(peakList);
        checkNotNull(peakSink);

        sink = peakSink;

        PeakSink<A> firstSink = initialize();

        firstSink.start(peakList.size());
        for(int i = 0; i < peakList.size(); i++) {

            firstSink.processPeak(peakList.getMz(i), peakList.getIntensity(i), peakList.getAnnotations(i));
        }
        firstSink.end();
    }

    public void process(TDoubleList mzList, TDoubleList intensityList, TIntObjectMap<List<A>> annotationMap, PeakSink<A> peakSink) {

        checkNotNull(mzList);
        checkNotNull(intensityList);
        checkNotNull(annotationMap);
        checkNotNull(peakSink);

        sink = peakSink;

        PeakSink<A> firstSink = initialize();

        int size = Math.min(mzList.size(), intensityList.size());
        firstSink.start(size);
        for(int i = 0; i < size; i++) {

            List<A> annotations = annotationMap.get(i);
            if(annotations == null) annotations = Collections.emptyList();

            firstSink.processPeak(mzList.get(i), intensityList.get(i), annotations);
        }
        firstSink.end();
    }

    public void process(Number[] mzArray, Number[] intensityArray, TIntObjectMap<List<A>> annotationMap, PeakSink<A> peakSink) {

        checkNotNull(mzArray);
        checkNotNull(intensityArray);
        checkNotNull(annotationMap);
        checkNotNull(peakSink);

        sink = peakSink;

        PeakSink<A> firstSink = initialize();

        int size = Math.min(mzArray.length, intensityArray.length);
        firstSink.start(size);
        for(int i = 0; i < size; i++) {

            List<A> annotations = annotationMap.get(i);
            if(annotations == null) annotations = Collections.emptyList();

            firstSink.processPeak(mzArray[i].doubleValue(), intensityArray[i].doubleValue(), annotations);
        }
        firstSink.end();
    }

    private PeakSink<A> initialize() {

        for(int i = 0; i < list.size() - 1; i++) {

            list.get(i).setSink(list.get(i + 1));
        }

        PeakSink<A> firstSink;
        if (!list.isEmpty()) {

            list.get(list.size() - 1).setSink(this.sink);
            firstSink = list.get(0);
        } else {

            firstSink = this.sink;
        }

        return firstSink;
    }

    public boolean isEmpty() {

        return list.isEmpty();
    }
}
