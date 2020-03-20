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


import com.google.common.base.Optional;
import newance.mzjava.ms.peaklist.PeakAnnotation;
import newance.mzjava.ms.peaklist.PeakList;
import newance.mzjava.ms.spectrum.MsnSpectrum;
import newance.mzjava.ms.spectrum.RetentionTime;
import newance.mzjava.ms.spectrum.ScanNumber;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.List;


/**
 * MGF Writer
 *
 * @author nikitin
 * @author Oliver Horlacher
 *
 * @version 1.0
 */
public class MgfWriter extends AbstractMgfWriter<PeakAnnotation, MsnSpectrum, MgfParameters> {

    public MgfWriter(Writer writer, PeakList.Precision precision) throws IOException {

        super(writer, new MgfParameters(precision));
    }

    public MgfWriter(File file, PeakList.Precision precision) throws IOException {

        this(new FileWriter(file, false), precision);
    }

    public MgfWriter(Writer writer, MgfParameters parameters) throws IOException {

        super(writer, parameters);
    }

    public MgfWriter(File file, MgfParameters parameters) throws IOException {

        this(new FileWriter(file, false), parameters);
    }

    public void write(MsnSpectrum peakList) throws IOException {

        writeMs2(peakList, peakList);
    }

    @Override
    protected void writeMs1Header(MgfParameters parameters, Writer writer) throws IOException {

        Optional<MsnSpectrum> headerData = parameters.getHeaderData();
        if (headerData.isPresent()) {

            MsnSpectrum spectrum = headerData.get();
            writeCharge(spectrum, spectrum, writer);
            writePeaks(spectrum, writer);
            writer.write("\n");
        }
    }

    @Override
    protected String extractTitle(MsnSpectrum spectrum) {

        return spectrum.getComment();
    }

    @Override
    protected List<ScanNumber> extractScanNumbers(MsnSpectrum spectrum) {

        return spectrum.getScanNumbers();
    }

    @Override
    protected List<RetentionTime> extractRetentionTime(MsnSpectrum spectrum) {

        return spectrum.getRetentionTimes();
    }

    @Override
    protected int[] extractChargeList(PeakList peakList, MsnSpectrum metaData) {

        return peakList.getPrecursor().getChargeList();
    }
}
