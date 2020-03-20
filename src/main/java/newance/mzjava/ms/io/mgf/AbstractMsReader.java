package newance.mzjava.ms.io.mgf;


import newance.mzjava.ms.peaklist.PeakAnnotation;
import newance.mzjava.ms.peaklist.PeakList;
import newance.mzjava.ms.peaklist.PeakProcessorChain;
import newance.mzjava.ms.peaklist.UnsortedPeakListException;
import newance.mzjava.ms.spectrum.Spectrum;

import java.io.IOException;
import java.io.LineNumberReader;
import java.io.Reader;
import java.net.URI;

import static com.google.common.base.Preconditions.checkNotNull;

public abstract class AbstractMsReader<A extends PeakAnnotation, S extends Spectrum<A>> implements IterativeReader<S> {

    // a parsed object waiting for a call to T next())
    private S waitingForDelivery;

    // potential next parsing thrown exception goes there
    private IOException exception;

    private final ParseContext context;

    protected final PeakList.Precision precision;
    protected boolean acceptUnsortedSpectra;
    private final String header;
    private final PeakProcessorChain<A> processorChain;

    protected AbstractMsReader(Reader reader, URI source, PeakList.Precision precision, PeakProcessorChain<A> processorChain) throws IOException {

        context = new ParseContext(reader, source);
        init(context);

        this.precision = precision;
        header = parseHeader(context);

        this.processorChain = processorChain;
    }

    protected void init(ParseContext parseContext) throws IOException {

    }

    /**
     * @return the next T-object or null if no more object to read.
     */
    protected abstract S parseNextEntry(ParseContext context) throws IOException;

    /**
     * @return the context
     */
    public ParseContext getContext() {

        return context;
    }

    @Override
    public boolean hasNext() {

        if (context.isEndOfParsing()) {

            return false;
        }

        // if no object awaiting -> parse next entry
        if (waitingForDelivery == null) {

            exception = parseNextEntryAndWaitForDelivery();
        }

        // if no object yet awaiting -> nothing more to parse
        if (waitingForDelivery == null && exception == null) {

            context.setEndOfParsing(true);

            return false;
        }

        return true;
    }

    @Override
    public S next() throws IOException {

        // parse the next entry if not yet done
        if (exception == null && waitingForDelivery == null) {

            exception = parseNextEntryAndWaitForDelivery();
        }

        if (exception != null) {

            // cannot make a new T
            waitingForDelivery = null;

            throw exception;
        }

        // store the next entry
        S spectrum = waitingForDelivery;

        // ready for next parsing
        waitingForDelivery = null;

        if (!processorChain.isEmpty())
            spectrum.apply(processorChain);

        return spectrum;
    }

    @Override
    public void close() throws IOException {

        context.getCurrentReader().close();
    }

    /**
     * Get the next entry and store it before next() is called.
     */
    private IOException parseNextEntryAndWaitForDelivery() {

        try {

            waitingForDelivery = parseNextEntry(context);

            if (waitingForDelivery != null) {

                context.incrementParsedEntryNumber();
            }
        } catch (Exception e) {

            // no need to continue the parsing
            waitingForDelivery = null;

            return new IOException(e.getMessage() +
                    ": cannot parse entry, context=[" + context + "]", e);
        }

        return null;
    }

    public void acceptUnsortedSpectra() {

        acceptUnsortedSpectra = true;
    }

    /**
     * This is called on reader initialization.  Subclasses should parse any information that occurs
     * before the first spectrum.  The information can returned as a string
     *
     * @param context the parser context
     * @return the data that is before the first spectrum
     * @throws IOException
     */
    protected abstract String parseHeader(ParseContext context) throws IOException;

    @SuppressWarnings("UnusedDeclaration")
    public String getHeader() {

        return header;
    }

    /**
     * Potentially non-order add peaks to spectrum.
     *
     * @param spectrum    the spectrum to add peaks into.
     * @param mzs         the peak's m/zs.
     * @param intensities the peak's intensities.
     * @param size        the spectrum size
     * @throws UnsortedPeakListException if reader is not allowed to add unsorted peaklists.
     */
    protected void addPeaksToSpectrum(S spectrum, double[] mzs, double[] intensities, int size) {

        checkNotNull(spectrum);

        if (!acceptUnsortedSpectra) {

            // catch lower-level runtime exception to rethrow appropriate higher-level runtime exception
            try {

                spectrum.addSorted(mzs, intensities, size);
            } catch (UnsortedPeakListException e) {

                throw new UnsortedPeakListException(context + ": cannot read unsorted spectrum. " +
                        "Call method 'acceptUnsortedSpectra()' if your reader has to deal with unsorted peaks (will be sorted internally).\"", e);
            }
        } else {

            // add peaks separately (potentially unordered)
            for (int i = 0; i < size; i++) {

                spectrum.add(mzs[i], intensities[i]);
            }
        }
    }

    /**
     * Potentially non-order add peak to spectrum.
     *
     * @param spectrum  the spectrum to add peak into.
     * @param mz        the peak's m/z.
     * @param intensity the peak's intensity.
     * @param context   the parse context
     * @return the index at which the peak was added
     * @throws UnsortedPeakListException if reader is not allowed to add unsorted peaklists.
     */
    protected int addPeakToSpectrum(Spectrum<A> spectrum, double mz, double intensity, ParseContext context) {

        checkNotNull(spectrum);

        if (!acceptUnsortedSpectra) {

            if (spectrum.isEmpty() || mz >= spectrum.getMz(spectrum.size() - 1)) {

                return spectrum.add(mz, intensity);
            } else {

                throw new UnsortedPeakListException(context + ": cannot read unsorted spectrum!", spectrum.size() - 1);
            }
        } else {

            return spectrum.add(mz, intensity);
        }

    }

    public static final class ParseContext {

        private LineNumberReader currentReader;

        private int numberOfEntryParsed;
        private boolean isEOPReached;
        private int columnNumber = -1;

        private final URI source;

        public ParseContext(URI source) {

            checkNotNull(source);

            this.source = source;
            isEOPReached = false;
        }

        public ParseContext(Reader reader, URI source) {

            this(source);

            if (reader != null) {

                if (reader instanceof LineNumberReader) {

                    currentReader = (LineNumberReader) reader;
                } else {

                    currentReader = new LineNumberReader(reader);
                }
            }
        }

        public LineNumberReader getCurrentReader() {

            return currentReader;
        }

        public URI getSource() {

            return source;
        }

        // @return the lineNumber or -1 if unavailable
        public final int getLineNumber() {

            if (currentReader == null) {

                return -1;
            }

            return currentReader.getLineNumber() + 1;
        }

        public final int getColumnNumber() {

            return columnNumber;
        }

        public void setColumnNumber(int columnNumber) {

            this.columnNumber = columnNumber;
        }

        public int getNumberOfParsedEntry() {

            return numberOfEntryParsed;
        }

        private void incrementParsedEntryNumber() {

            this.numberOfEntryParsed++;
        }

        public boolean isEndOfParsing() {

            return isEOPReached;
        }

        public void setEndOfParsing(boolean bool) {

            isEOPReached = bool;
        }

        public String toString() {

            StringBuilder sb = new StringBuilder();

            sb.append("source=");
            sb.append(source);

            if (getLineNumber() != -1) {
                sb.append(", line=");
                sb.append(getLineNumber());
            }
            if (getColumnNumber() != -1) {

                sb.append(", column=");
                sb.append(getColumnNumber());
            }

            return sb.toString();
        }
    }
}
