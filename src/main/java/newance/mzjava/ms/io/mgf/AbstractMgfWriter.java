/**
 * Copyright (c) 2010, SIB. All rights reserved.
 *
 * SIB (Swiss Institute of Bioinformatics) - http://www.isb-sib.ch Host -
 * http://mzjava.expasy.org
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer. Redistributions in binary
 * form must reproduce the above copyright notice, this list of conditions and
 * the following disclaimer in the documentation and/or other materials provided
 * with the distribution. Neither the name of the SIB/GENEBIO nor the names of
 * its contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL SIB/GENEBIO BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package newance.mzjava.ms.io.mgf;

import newance.mzjava.ms.peaklist.Peak;
import newance.mzjava.ms.peaklist.PeakAnnotation;
import newance.mzjava.ms.peaklist.PeakList;
import newance.mzjava.ms.peaklist.Polarity;
import newance.mzjava.ms.spectrum.*;

import java.io.IOException;
import java.io.Writer;
import java.text.NumberFormat;
import java.util.List;

import static com.google.common.base.Preconditions.checkNotNull;

/**
 * Abstract implementation of an MGF writer. This class takes care of writing the standard
 * structure of an MGF file. For each spectrum the following tags are written
 *
 * BEGIN IONS
 * TITLE=
 * PEPMASS=
 * CHARGE=
 * RTINSECONDS=
 * peak mz\tpeak intensity
 * END IONS
 *
 * Subclasses have to implement methods to extract the charge list, retention time, scan number
 * and title to fill in the values for the tags.
 *
 * By default this class does not write any ms1 header data, however subclasses can generate a ms1 header by
 * using a Parameters implementation that contains data for the header and overriding writeMs1Header(P, Writer)
 * method which is called by the constructor.
 *
 * To write a PeakList and associated meta data subclasses have to call writeMs2(PeakList peakList, M metaData).
 *
 * @author Oliver Horlacher
 * @author nikitin
 * @version sqrt -1
 *
 * @param <A> the annotations of the peaks in the peak list
 * @param <M> the object that holds the meta data, can be the same object as the peak list
 * @param <P> the object that holds the parameters for initialising this mgf writer
 */
public abstract class AbstractMgfWriter<A extends PeakAnnotation, M, P extends AbstractMgfWriter.Parameters> {

    public static interface Parameters {

        NumberFormat getMzFormat();

        NumberFormat getIntensityFormat();

        NumberFormat getRetentionTimeFormat();
    }

    private final NumberFormat retentionTimeFormat;
    private final NumberFormat mzFormat;
    private final NumberFormat intensityFormat;

    private final Writer writer;

    public AbstractMgfWriter(Writer writer, P parameters) throws IOException {

        checkNotNull(writer);
        checkNotNull(parameters);

        mzFormat = parameters.getMzFormat();
        intensityFormat = parameters.getIntensityFormat();
        retentionTimeFormat = parameters.getRetentionTimeFormat();

        checkNotNull(mzFormat);
        checkNotNull(intensityFormat);
        checkNotNull(retentionTimeFormat);

        this.writer = writer;

        writeMs1Header(parameters, writer);
    }

    protected void writeMs1Header(P parameters, Writer writer) throws IOException {

        //By default do nothing
    }

    public void writeComment(String comment) throws IOException {

        writer.write("# ");
        writer.write(comment);
        writer.write("\n");
    }

    protected void writeMs2(PeakList<? extends A> peakList, M metaData) throws IOException {

        writer.write("BEGIN IONS\n");
        writeMs2Header(peakList, metaData, writer);
        writePeaks(peakList, writer);
        writer.write("END IONS\n");
        writer.write("\n");
    }

    protected void writeMs2Header(PeakList peakList, M metaData, Writer writer) throws IOException {

        writeTitle(metaData, writer);
        writePepmass(peakList, writer);
        writeCharge(peakList, metaData, writer);
        writeScans(metaData, writer);
        writeRtInSeconds(metaData, writer);
    }

    protected void writeTitle(M metaData, Writer writer) throws IOException {

        String comment = extractTitle(metaData);
        if (comment.length() > 0) {
            writer.write("TITLE=");
            writer.write(comment);
            writer.write("\n");
        }
    }

    protected abstract String extractTitle(M metaData);

    protected void writePepmass(PeakList peakList, Writer writer) throws IOException {

        Peak precursor = peakList.getPrecursor();

        writer.write("PEPMASS=");

        writer.write(mzFormat.format(precursor.getMz()));
        writer.write("\t");

        writer.write(intensityFormat.format(precursor.getIntensity()));
        writer.write("\n");
    }

    protected void writeScans(M metaData, Writer writer) throws IOException {

        List<ScanNumber> scanNumList = extractScanNumbers(metaData);
        if (!scanNumList.isEmpty()) {

            writer.write("SCANS=");
            boolean first = true;
            for (ScanNumber sn : scanNumList) {

                if (first) {

                    first = false;
                } else {

                    writer.write(",");
                }
                if (sn instanceof ScanNumberDiscrete) {

                    writer.write(Integer.toString(sn.getValue()));
                } else if (sn instanceof ScanNumberInterval) {

                    writer.write(Integer.toString(sn.getMinScanNumber()));
                    writer.write('-');
                    writer.write(Integer.toString(sn.getMaxScanNumber()));
                }
            }
            writer.write('\n');
        }
    }

    protected abstract List<ScanNumber> extractScanNumbers(M metaData);

    protected void writeRtInSeconds(M metaData, Writer writer) throws IOException {

        List<RetentionTime> retentionTimeList = extractRetentionTime(metaData);
        if (!retentionTimeList.isEmpty()) {

            writer.write("RTINSECONDS=");
            boolean first = true;
            for (RetentionTime rt : retentionTimeList) {

                if (first) {

                    first = false;
                } else {

                    writer.write(",");
                }

                if (rt instanceof RetentionTimeDiscrete) {

                    writer.write(retentionTimeFormat.format(rt.getTime()));
                } else {

                    writer.write(retentionTimeFormat.format(rt.getMinRetentionTime()));
                    writer.write('-');
                    writer.write(retentionTimeFormat.format(rt.getMaxRetentionTime()));
                }
            }
            writer.write('\n');
        }
    }

    protected abstract List<RetentionTime> extractRetentionTime(M metaData);

    protected void writePeaks(PeakList<? extends A> peakList, Writer writer) throws IOException {

        for (int i = 0; i < peakList.size(); i++) {

            double mz = peakList.getMz(i);
            double intensity = peakList.getIntensity(i);

            writer.write(mzFormat.format(mz));
            writer.write("\t");
            writer.write(intensityFormat.format(intensity));
            writePeakAnnotations(peakList.getAnnotations(i), writer);
            writer.write("\n");
        }
    }

    /**
     * Subclasses that want to add annotation information to a peak can override this method.
     *
     * Before this method is called the peak mz and intensity were written. Writing \n is
     * taken care of in the method that calls this.
     *
     * @param annotations the annotations for the peak being written
     * @param writer the writer
     */
    @SuppressWarnings("UnusedParameters")
    protected void writePeakAnnotations(List<? extends A> annotations, Writer writer) {}

    protected void writeCharge(PeakList peakList, M metaData, Writer writer) throws IOException {

        int[] charges = extractChargeList(peakList, metaData);

        if (charges.length > 0) {

            writer.write("CHARGE=");

            Polarity polarity = peakList.getPrecursor().getPolarity();
            boolean first = true;
            for (int charge : charges) {

                if (first) {

                    first = false;
                } else {

                    writer.write(", ");
                }
                writer.write(String.valueOf(Math.abs(charge)));

                if (polarity == Polarity.POSITIVE) {
                    writer.write('+');
                } else if (polarity == Polarity.NEGATIVE) {

                    writer.write('-');
                }
            }
            writer.write('\n');
        }
    }

    protected abstract int[] extractChargeList(PeakList peakList, M metaData);

    public void close() throws IOException {

        writer.close();
    }
}
