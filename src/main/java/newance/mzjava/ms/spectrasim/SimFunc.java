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

import newance.mzjava.ms.peaklist.PeakAnnotation;
import newance.mzjava.ms.peaklist.PeakList;

/**
 * Interface for any class that can calculate the similarity between two PeakLists
 *
 * @author Oliver Horlacher
 * @version 1.0
 */
public interface SimFunc<X extends PeakAnnotation, Y extends PeakAnnotation> {

    /**
     * Calculate the similarity score of the two vectors.
     *
     * Returns the similarity score or NaN if the score cannot be calculated
     *
     * @param plX the fist peak list
     * @param plY the second peak list
     * @return the similarity score or NaN if the score cannot be calculated
     */
    double calcSimilarity(PeakList<X> plX, PeakList<Y> plY);

    /**
     * Returns the total number peaks in plX and plY from the last call to calcSimilarity
     *
     * @return the total number peaks in plX and plY from the last call to calcSimilarity
     */
    int getTotalPeakCount();

    /**
     * Return the best score that can be returned by calls to calcSimilarity.
     * For similarity functions for which the best score has no upper limit
     * Double.POSITIVE_INFINITY is returned and where the best score has no lower
     * limit Double.NEGATIVE_INFINITY is returned.
     *
     * @return the best score that can be returned by calls to calcSimilarity
     */
    double getBestScore();

    /**
     * Return the worst the score that can be returned by calls to calcSimilarity.
     * For similarity functions for which the worst score has no upper limit
     * Double.POSITIVE_INFINITY is returned and where the worst score has no lower
     * limit Double.NEGATIVE_INFINITY is returned.
     *
     * @return the worst score that can be returned by calls to calcSimilarity
     */
    double getWorstScore();
}
