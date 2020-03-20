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

/**
 * A scan number can either be a discrete number or a number interval.
 * <p/>
 * This interface serves to hide this complexity from code for which is not relevant.
 *
 * @author Oliver Horlacher
 * @version 1.0
 */
public interface ScanNumber {

    public static final ScanNumber SCAN_NUMBER_UNKNOWN = new ScanNumberDiscrete(-1);

    /**
     * Return the scan number.
     *
     * @return the scan number
     */
    int getValue();

    /**
     * Return the smallest scan number
     *
     * @return the smallest scan number
     */
    int getMinScanNumber();

    /**
     * Return the largest scan number
     *
     * @return the largest scan number
     */
    int getMaxScanNumber();

    /**
     * Return true if this scan number contains the <code>scanNumber</code>, false otherwise
     *
     * @param scanNumber the scan number to test
     * @return true if this scan number contains the <code>scanNumber</code>, false otherwise
     */
    boolean contains(ScanNumber scanNumber);
}
