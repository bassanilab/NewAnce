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
package newance.mzjava.ms.spectrum;

import newance.mzjava.mol.Peptide;
import newance.mzjava.ms.peaklist.PeakList;
import newance.mzjava.ms.peaklist.peaktransformer.IdentityPeakProcessor;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

/**
 * @author Oliver Horlacher
 * @version 1.0
 */
public class PeptideSpectrumTest {

    @Test
    public void testEquals() throws Exception {

        PeptideSpectrum spectrum1 = new PeptideSpectrum(Peptide.parse("CERVILAS"), 1, PeakList.Precision.FLOAT);

        PeptideSpectrum spectrum2 = new PeptideSpectrum(Peptide.parse("CERVILAS"), 1, PeakList.Precision.FLOAT);

        PeptideSpectrum spectrum3 = new PeptideSpectrum(Peptide.parse("PEPTIDE"), 2, PeakList.Precision.FLOAT);

        assertEquals(true, spectrum1.equals(spectrum2));
        assertEquals(false, spectrum2.equals(spectrum3));
    }

    @Test
    public void testHashCode() throws Exception {

        PeptideSpectrum spectrum1 = new PeptideSpectrum(Peptide.parse("CERVILAS"), 1, PeakList.Precision.FLOAT);

        PeptideSpectrum spectrum2 = new PeptideSpectrum(Peptide.parse("CERVILAS"), 1, PeakList.Precision.FLOAT);

        PeptideSpectrum spectrum3 = new PeptideSpectrum(Peptide.parse("PEPTIDE"), 2, PeakList.Precision.FLOAT);

        assertEquals(spectrum1.hashCode(), spectrum2.hashCode());
        assertFalse("Spectra have different score values therefore they should have different hash codes",
                spectrum2.hashCode() == spectrum3.hashCode());
    }

    @Test
    public void testCopy() throws Exception {

        Peptide peptide = Peptide.parse("CERVILAS");
        PeptideSpectrum spectrum = new PeptideSpectrum(peptide, 3, PeakList.Precision.FLOAT);

        PeptideSpectrum copy = spectrum.copy(new IdentityPeakProcessor<PepFragAnnotation>());

        assertEquals(peptide, copy.getPeptide());
    }
}
