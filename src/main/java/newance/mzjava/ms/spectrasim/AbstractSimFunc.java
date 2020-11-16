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
package newance.mzjava.ms.spectrasim;


import newance.mzjava.mol.modification.Tolerance;
import newance.mzjava.ms.peaklist.PeakAnnotation;
import newance.mzjava.ms.peaklist.PeakList;
import newance.mzjava.ms.spectrasim.peakpairprocessor.*;

/**
 * Abstract similarity function. This class takes care of vectorizing the PeakList pairs.
 *
 * Vectorizing aligns the two PeakLists and then processes the aligned peak pairs using the
 * PeakPairProcessor chain.
 *
 * @author Oliver Horlacher
 * @version 1.0
 */
public abstract class AbstractSimFunc<X extends PeakAnnotation, Y extends PeakAnnotation> implements SimFunc<X, Y>, PeakPairSink<X, Y> {

    private final PeakPairProcessorChain<X, Y> peakPairProcessorChain;

    /**
     * Creates a AbstractSimFunc that uses a DefaultPeakListAligner and clear processor chain.
     *
     * @param tolerance the tolerance to use when aligning the two PeakLists
     */
    protected AbstractSimFunc(Tolerance tolerance) {

        this.peakPairProcessorChain = new PeakPairProcessorChain<X, Y>(new DefaultPeakListAligner<X, Y>(tolerance));
    }

    /**
     * Creates a AbstractSimFunc that uses a DefaultPeakListAligner and clear processor chain.
     *
     * @param peakListAligner object to align the two PeakLists
     */
    protected AbstractSimFunc(PeakListAligner<X, Y> peakListAligner) {

        this.peakPairProcessorChain = new PeakPairProcessorChain<X, Y>(peakListAligner);
    }

    /**
     * Creates a AbstractSimFunc that uses the supplied PeakListAligner and PeakPairProcessor chain.
     *
     * @param peakListAligner the object that is used to align the PeakList pairs
     * @param chain PeakPairProcessor chain that is used to pre-process the aligned peaks before calculating the similarity
     */
    protected AbstractSimFunc(PeakListAligner<X, Y> peakListAligner, PeakPairProcessor<X, Y>... chain) {

        this.peakPairProcessorChain = new PeakPairProcessorChain<X, Y>(peakListAligner, chain);
    }

    /**
     * Aligns the two PeakLists and then processes the aligned peak tupplets using the
     * PeakPairProcessor chain.
     *
     * @param peakListX the first peak list
     * @param peakListY the second peak list
     */
    public void vectorize(PeakList<X> peakListX, PeakList<Y> peakListY) {

        peakPairProcessorChain.process(this, peakListX, peakListY);
    }

    @Override
    public void begin(PeakList<X> xPeakList, PeakList<Y> yPeakList) {}

    @Override
    public void end() {}
}
