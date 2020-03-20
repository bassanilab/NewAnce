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


import newance.mzjava.ms.peaklist.Copyable;

import java.util.ArrayList;
import java.util.Collections;

/**
 * @author Oliver Horlacher
 * @version 1.0
 */
public class ScanNumberList extends ArrayList<ScanNumber> implements Copyable<ScanNumberList> {

    /**
     * Construct an empty ScanNumberList
     */
    public ScanNumberList() {
    }

    /**
     * Copy constructor
     *
     * @param scanNumbers the scan number list to copy
     */
    public ScanNumberList(ScanNumberList scanNumbers) {

        super(scanNumbers.size());

        for(ScanNumber scanNumber : scanNumbers) {

            add(scanNumber);
        }
    }

    /**
     * Construct a ScanNumberList that contains a ScanNumberDiscrete for each entry in <code>scanNumbers</code>
     *
     * @param scanNumbers the initial content of this ScanNumberList
     */
    public ScanNumberList(int... scanNumbers) {

        super(scanNumbers.length);

        for(int scanNumber : scanNumbers) {

            add(new ScanNumberDiscrete(scanNumber));
        }
    }

    /**
     * Construct a ScanNumberList that contains the ScanNumber's in <code>scanNumbers</code>
     *
     * @param scanNumbers the array whose elements are to be placed in this ScanNumberList
     */
    public ScanNumberList(ScanNumber... scanNumbers) {

        Collections.addAll(this, scanNumbers);
    }

    /**
     * Return the first scan number in this ScanNumberList
     *
     * @return the first scan number in this ScanNumberList
     * @throws IndexOutOfBoundsException if this ScanNumberList is empty
     */
    public ScanNumber getFirst(){

        return get(0);
    }

    /**
     * Return the last scan number in this ScanNumberList
     *
     * @return the last scan number in this ScanNumberList
     * @throws IndexOutOfBoundsException if this ScanNumberList is empty
     */
    public ScanNumber getLast() {

        return get(size() - 1);
    }

    @Override
    public ScanNumberList copy() {

        return new ScanNumberList(this);
    }

    /**
     * Add a ScanNumberInterval that starts at <code>minScanNumber</code> and ends at <code>maxScanNumber</code> to
     * this ScanNumberList
     *
     * @param minScanNumber the min scan number
     * @param maxScanNumber the max scan number
     */
    public void add(int minScanNumber, int maxScanNumber) {

        add(new ScanNumberInterval(minScanNumber, maxScanNumber));
    }

    /**
     * Add a ScanNumberDiscrete to this ScanNumberList
     *
     * @param scanNumber the scan number of the ScanNumberDiscrete that is to be added to this ScanNumberList
     */
    public void add(int scanNumber) {

        add(new ScanNumberDiscrete(scanNumber));
    }

    @Override
    public boolean contains(Object o) {

        if(!(o instanceof ScanNumber)) return false;

        ScanNumber nbr = (ScanNumber)o;
        for(ScanNumber scanNumber : this) {

            if(scanNumber.contains(nbr)) return true;
        }

        return false;
    }
}
