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
package newance.mzjava.ms.consensus;

import newance.mzjava.ms.peaklist.*;
import newance.mzjava.ms.peaklist.peakfilter.AbstractMergePeakFilter;
import newance.mzjava.ms.spectrum.LibPeakAnnotation;
import newance.mzjava.ms.spectrum.Spectrum;

import java.util.*;

/**
 * @author Oliver Horlacher
 * @author Markus Muller
 * @version sqrt -1
 */
public class ConsensusSpectrum<A extends LibPeakAnnotation> extends Spectrum<A> {

    protected String name = "";
    protected double simScoreMean, simScoreStdev;
    protected double precursorMzMean, precursorMzStdev;

    protected Set<UUID> memberIds;

    public ConsensusSpectrum(int initialCapacity, Precision precision, Set<UUID> memberIds) {

        super(initialCapacity, precision);
        this.memberIds = Collections.unmodifiableSet(new HashSet<>(memberIds));
    }

    public ConsensusSpectrum(ConsensusSpectrum<A> consensusSpectrum, PeakProcessor<A, A> peakProcessor) {

        super(consensusSpectrum, peakProcessor);

        copyFields(consensusSpectrum);
    }

    public ConsensusSpectrum(ConsensusSpectrum<A> consensusSpectrum, PeakProcessorChain<A> peakProcessorChain) {

        super(consensusSpectrum, peakProcessorChain);

        copyFields(consensusSpectrum);
    }

    private void copyFields(ConsensusSpectrum<A> consensusSpectrum) {

        memberIds = Collections.unmodifiableSet(new HashSet<>(consensusSpectrum.memberIds));
        name = consensusSpectrum.name;
        simScoreMean = consensusSpectrum.simScoreMean;
        simScoreStdev = consensusSpectrum.simScoreStdev;
        precursorMzMean = consensusSpectrum.precursorMzMean;
        precursorMzStdev = consensusSpectrum.precursorMzStdev;
    }

    public void addMemberIds(Collection<UUID> memberIds) {

        this.memberIds = Collections.unmodifiableSet(new HashSet<>(memberIds));
    }

    public Set<UUID> getMemberIds() {

        return memberIds;
    }

    public int nrOfMembers() {

        return memberIds.size();
    }

    public void setName(String name) {

        this.name = name;
    }

    public String getName() {

        return name;
    }

    public ConsensusSpectrum setScoreStats(double simScoreMean, double simScoreStdev) {

        this.simScoreMean = simScoreMean;
        this.simScoreStdev = simScoreStdev;

        return this;
    }

    public ConsensusSpectrum setPrecursorStats(double precursorMean, double precursorStdev) {

        this.precursorMzMean = precursorMean;
        this.precursorMzStdev = precursorStdev;

        return this;
    }

    public double getSimScoreMean() {

        return simScoreMean;
    }

    public double getSimScoreStdev() {

        return simScoreStdev;
    }

    public double getPrecursorMzMean() {

        return precursorMzMean;
    }

    public double getPrecursorMzStdev() {

        return precursorMzStdev;
    }

    @Override
    public ConsensusSpectrum<A> copy(PeakProcessor<A, A> peakProcessor) {

        return new ConsensusSpectrum<>(this, peakProcessor);
    }

    @Override
    public ConsensusSpectrum<A> copy(PeakProcessorChain<A> peakProcessorChain) {

        return new ConsensusSpectrum<>(this, peakProcessorChain);
    }

    @Override
    public boolean equals(Object o) {

        if (this == o) return true;
        if (!(o instanceof ConsensusSpectrum)) return false;
        if (!super.equals(o)) return false;

        ConsensusSpectrum that = (ConsensusSpectrum) o;

        return Double.compare(that.precursorMzMean, precursorMzMean) == 0 &&
                Double.compare(that.precursorMzStdev, precursorMzStdev) == 0 &&
                Double.compare(that.simScoreMean, simScoreMean) == 0 &&
                Double.compare(that.simScoreStdev, simScoreStdev) == 0 &&
                name.equals(that.name);
    }

    @Override
    public int hashCode() {

        int result = super.hashCode();
        long temp;
        result = 31 * result + name.hashCode();
        temp = Double.doubleToLongBits(simScoreMean);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(simScoreStdev);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(precursorMzMean);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(precursorMzStdev);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        return result;
    }

    /**
     * @author Markus Muller
     * @version 0.0
     */
    public static class Builder extends AbstractConsensusSpectrumBuilder<LibPeakAnnotation, Builder, ConsensusSpectrum<LibPeakAnnotation>> {

        protected Builder() {

        }

        protected void calculatePrecursor(ConsensusSpectrum<? extends LibPeakAnnotation> consensus) {

            Peak precursor;

            double meanMz = 0.0;
            double meanIntensity = 0.0;
            double sigmaMz = 0.0;

            if (spectra.size() > 1) {

                int n = 0;
                Set<Integer> charges = new HashSet<>();


                for (PeakList pl : spectra) {
                    n++;
                    double delta = pl.getPrecursor().getMz() - meanMz;
                    meanMz += delta / n;
                    sigmaMz += delta * (pl.getPrecursor().getMz() - meanMz);
                    delta = pl.getPrecursor().getIntensity() - meanIntensity;
                    meanIntensity += delta / n;
                    int[] chargeArray = pl.getPrecursor().getChargeList();
                    for (int aChargeArray : chargeArray) {
                        charges.add(aChargeArray);
                    }
                }

                if (n > 1) sigmaMz /= n - 1;

                int[] chargeArray = new int[charges.size()];
                int i = 0;
                for (Integer ch : charges) {
                    chargeArray[i++] = ch;
                }

                precursor = new Peak(meanMz, meanIntensity, chargeArray);

            } else {
                precursor = spectra.get(0).getPrecursor();

            }

            consensus.setPrecursor(precursor);
            consensus.setPrecursorStats(meanMz, Math.sqrt(sigmaMz));
        }

        public static Builder getBuilder() {

            return new Builder();
        }

        public ConsensusSpectrum<LibPeakAnnotation> build() {

            Set<UUID> memberIDs = new HashSet<>();
            if (spectra.isEmpty()) return new ConsensusSpectrum<>(0, Precision.DOUBLE, memberIDs);

            for (PeakList<PeakAnnotation> spectrum : spectra) {
                memberIDs.add(spectrum.getId());
            }
            ConsensusSpectrum<LibPeakAnnotation> consensus = new ConsensusSpectrum<>(initialCapacity, Precision.DOUBLE, memberIDs);
            calculatePrecursor(consensus);
            consensus.setMsLevel(msLevel);

            MergePeakFilter mergePeakFilter = new MergePeakFilter(fragMzTolerance, maxMzClusterWidth, intCombMeth, spectra.size());

            if (spectra.size() > 1) {
                PeakListMerger<PeakAnnotation> peakListMerger = new PeakListMerger<>();
                peakListMerger
                        .setSink(mergePeakFilter)
                        .setSink(new ConsensusPeakSink<>(consensus, peakFraction, minPeakCount));
                peakListMerger.merge(spectra);
            } else {
                PeakList spectrum = spectra.get(0);
                consensus.ensureCapacity(spectrum.size());
                for (int i = 0; i < spectrum.size(); i++) {
                    consensus.add(spectrum.getMz(i), spectrum.getIntensity(i), mergePeakFilter.createPeakAnnotation(1, 0.0, 0.0));
                }
                consensus.trimToSize();
            }

            return consensus;
        }

        private class MergePeakFilter extends AbstractMergePeakFilter<PeakAnnotation, LibPeakAnnotation> {

            public MergePeakFilter(double mzMaxDiff, double maxMzClusterWidth, IntensityMode intCombMeth, int totSpectrumCount) {

                super(mzMaxDiff, maxMzClusterWidth, intCombMeth, totSpectrumCount);
            }

            @Override
            protected LibPeakAnnotation createPeakAnnotation(int peakCount, double mzStd, double intensityStd) {

                return new LibPeakAnnotation(peakCount, mzStd, intensityStd);
            }
        }
    }
}