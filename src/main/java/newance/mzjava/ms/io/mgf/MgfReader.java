package newance.mzjava.ms.io.mgf;

import newance.mzjava.ms.peaklist.PeakAnnotation;
import newance.mzjava.ms.peaklist.PeakList;
import newance.mzjava.ms.peaklist.PeakProcessorChain;
import newance.mzjava.ms.spectrum.MsnSpectrum;
import newance.mzjava.ms.spectrum.RetentionTimeList;
import newance.mzjava.ms.spectrum.ScanNumberList;
import newance.mzjava.ms.spectrum.URIBuilder;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.net.URI;
import java.util.ArrayList;
import java.util.List;

import static com.google.common.base.Preconditions.checkNotNull;

/**
 * @author fnikitin
 * @author Oliver Horlacher
 */
public class MgfReader extends AbstractMgfReader<PeakAnnotation, MsnSpectrum> {

    private final List<TitleParser> titleParsers;
    private TitleParser currentTitleParser;

   public MgfReader(File file, TitleParser titleParser) throws IOException {

        super(new FileReader(file), file.toURI(), PeakList.Precision.DOUBLE, new PeakProcessorChain<>());
        titleParsers = new ArrayList<>();
        titleParsers.add(titleParser);
        currentTitleParser = titleParsers.get(0);
    }

    public MgfReader(Reader reader, TitleParser titleParser) throws IOException {

        super(reader, URIBuilder.UNDEFINED_URI, PeakList.Precision.DOUBLE, new PeakProcessorChain<>());
        titleParsers = new ArrayList<>();
        titleParsers.add(titleParser);
        currentTitleParser = titleParsers.get(0);
    }

    public MgfReader(File file, TitleParser titleParser, PeakProcessorChain<PeakAnnotation> processorChain) throws IOException {

        super(new FileReader(file), file.toURI(), PeakList.Precision.DOUBLE, processorChain);
        titleParsers = new ArrayList<>();
        titleParsers.add(titleParser);
        currentTitleParser = titleParsers.get(0);
    }

    public void addTitleParser(TitleParser parser) {

        checkNotNull(parser);
        titleParsers.add(parser);
        currentTitleParser = parser;
    }

    @Override
    protected MsnSpectrum newSpectrum(ParseContext context, PeakList.Precision precision) {

        MsnSpectrum spectrum = new MsnSpectrum(precision);
        spectrum.setSpectrumIndex(context.getNumberOfParsedEntry());
        spectrum.setSpectrumSource(context.getSource());

        return spectrum;
    }

    @Override
    protected void setRetentionTimes(MsnSpectrum spectrum, RetentionTimeList retentionTimeList) {

        spectrum.addRetentionTimes(retentionTimeList);
    }

    @Override
    protected void setScanNumbers(MsnSpectrum spectrum, ScanNumberList scanNumbers) {

        spectrum.addScanNumbers(scanNumbers);
    }

    protected boolean parseTitleTag(String value, MsnSpectrum spectrum) {

        boolean parsed = currentTitleParser.parseTitle(value, spectrum);

        if (!parsed) {

            for (TitleParser titleParser : titleParsers) {

                if (titleParser.parseTitle(value, spectrum)) {

                    currentTitleParser = titleParser;
                    break;
                }
            }
        }

        spectrum.setComment(value);

        return parsed;
    }

    /**
     * Overridden so that the compiler knows that MsnSpectra are returned by next().
     * If this method is missing the compiler thinks Spectrum are returned
     *
     * @inheritDoc
     */
    @Override
    public MsnSpectrum next() throws IOException {

        return super.next();
    }
}
