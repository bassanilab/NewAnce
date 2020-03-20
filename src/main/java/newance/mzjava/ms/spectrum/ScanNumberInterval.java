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
 * @author Oliver Horlacher
 * @version 1.0
 */
public class ScanNumberInterval implements ScanNumber {
    
    private final int minScanNumber, maxScanNumber;

    public ScanNumberInterval(int minScanNumber, int maxScanNumber) {

        if(minScanNumber > maxScanNumber) throw new IllegalArgumentException("Min scan number has to be <= than max scan number.  minScanNumber was " + minScanNumber + " and maxScanNumber was " + maxScanNumber);

        this.minScanNumber = minScanNumber;
        this.maxScanNumber = maxScanNumber;
    }

    @Override
    public int getValue() {
        
        return minScanNumber + (maxScanNumber - minScanNumber)/2;
    }

    public int getMinScanNumber() {

        return minScanNumber;
    }

    public int getMaxScanNumber() {

        return maxScanNumber;
    }

    @Override
    public boolean contains(ScanNumber scanNumber) {

        return scanNumber.getMinScanNumber() >= minScanNumber && scanNumber.getMaxScanNumber() <= maxScanNumber;
    }

    @Override
    public boolean equals(Object o) {

        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        ScanNumberInterval that = (ScanNumberInterval) o;

        return maxScanNumber == that.maxScanNumber &&
                minScanNumber == that.minScanNumber;

    }

    @Override
    public int hashCode() {

        int result = minScanNumber;
        result = 31 * result + maxScanNumber;
        return result;
    }

    @Override
    public String toString() {
        return "ScanNumberInterval{" +
                "scanNumber=" + minScanNumber +
                "-" + maxScanNumber +
                '}';
    }
}
