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
import newance.mzjava.ms.peaklist.PeakProcessor;
import newance.mzjava.ms.peaklist.PeakProcessorChain;
import newance.mzjava.ms.peaklist.peaktransformer.IdentityPeakProcessor;

import java.util.*;

import static com.google.common.base.Preconditions.checkArgument;
import static com.google.common.base.Preconditions.checkNotNull;

/**
 * @author Oliver Horlacher
 * @version 1.0
 */
public class PeptideSpectrum extends Spectrum<PepFragAnnotation> {

    private final Peptide peptide;
    private final int charge;
    private final Set<String> proteinAccessionNumbers = new LinkedHashSet<String>();

    public PeptideSpectrum(Peptide peptide, int charge) {

        this(peptide, charge, 0, Precision.DOUBLE);
    }

    public PeptideSpectrum(Peptide peptide, int charge, Precision precision) {

        this(peptide, charge, 0, precision);
    }

    public PeptideSpectrum(Peptide peptide, int charge, int initialCapacity, Precision precision) {

        super(initialCapacity, precision);

        checkNotNull(peptide);
        checkArgument(charge != 0);

        this.peptide = peptide;
        this.charge = charge;

        setMsLevel(2);
    }

    public PeptideSpectrum(Peptide peptide, int charge, int initialCapacity, double constantIntensity, Precision precision) {

        super(initialCapacity, constantIntensity, precision);

        checkNotNull(peptide);
        this.peptide = peptide;
        this.charge = charge;

        setMsLevel(2);
    }

    protected PeptideSpectrum(PeptideSpectrum src, PeakProcessor<PepFragAnnotation,PepFragAnnotation> peakProcessor) {

        super(src,peakProcessor);

        peptide = src.peptide;
        charge = src.charge;
    }

    protected PeptideSpectrum(PeptideSpectrum src, PeakProcessorChain<PepFragAnnotation> peakProcessorChain) {

        super(src,peakProcessorChain);

        peptide = src.peptide;
        charge = src.charge;
    }

    public PeptideSpectrum(Spectrum<PepFragAnnotation> src, Peptide peptide, int charge) {

        super(src, new IdentityPeakProcessor<PepFragAnnotation>());
        this.peptide = peptide;
        this.charge = charge;

        setMsLevel(2);
    }

    @Override
    public boolean equals(Object o) {

        if (this == o) return true;
        if (!(o instanceof PeptideSpectrum)) return false;
        if (!super.equals(o)) return false;

        PeptideSpectrum that = (PeptideSpectrum) o;

        return !(peptide != null ? !peptide.equals(that.peptide) : that.peptide != null) && charge == that.charge;

    }

    @Override
    public int hashCode() {

        int result = super.hashCode();
        result = 31 * result + (peptide != null ? peptide.hashCode() : 0);
        result = 31 * result + charge;
        return result;
    }

    @Override
    public PeptideSpectrum copy(PeakProcessor<PepFragAnnotation, PepFragAnnotation> peakProcessor) {

        return new PeptideSpectrum(this, peakProcessor);
    }

    @Override
    public PeptideSpectrum copy(PeakProcessorChain<PepFragAnnotation> peakProcessorChain) {

        return new PeptideSpectrum(this,peakProcessorChain);
    }

    public Peptide getPeptide() {

        return peptide;
    }

    public int getCharge() {

        return charge;
    }

    public void addProteinAccessionNumbers(String... acc) {

        Collections.addAll(proteinAccessionNumbers, acc);
    }

    public void addProteinAccessionNumbers(Collection<String> accessionNumbers) {

        proteinAccessionNumbers.addAll(accessionNumbers);
    }

    public List<String> getProteinAccessionNumbers() {

        return Collections.unmodifiableList(new ArrayList<String>(proteinAccessionNumbers));
    }
}
