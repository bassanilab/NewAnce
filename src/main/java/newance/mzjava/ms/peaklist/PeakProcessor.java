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
package newance.mzjava.ms.peaklist;

/**
 * A peak sink that is one element in a peak processor chain.
 * <p>
 * The peak processor is used to create a processor chain to process a PeakList.
 * For example, to extract the top 100 peaks from a PeakList and then normalize to the
 * second most intense peak a processor chain is created as follows:
 * <p>
 * <code>
 *      NPeaksFilter topPeaksFilter = new TopPeakFilter(100);
 *      NthPeakNormalizer scaleToPeakProcessor = new NthPeakNormalizer(2);
 *
 *      PeakProcessorChain processorChain = new PeakProcessorChain()
 *                                      .add(topPeaksFilter)
 *                                      .add(scaleToPeakProcessor);
 *
 *      //To process in place
 *      peakList.apply(processorChain);
 *
 *      //To create a copy that is processed
 *      peakList.copy(processorChain);
 * </code>
 *
 * @author Oliver Horlacher
 * @version 1.0
 */
public interface PeakProcessor<IN extends PeakAnnotation, OUT extends PeakAnnotation> extends PeakSource<OUT>, PeakSink<IN> {
}

