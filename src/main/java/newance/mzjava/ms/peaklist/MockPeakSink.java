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

import gnu.trove.list.array.TDoubleArrayList;

import java.util.List;

/**
 * @author Oliver Horlacher
 * @version 1.0
 */
public class MockPeakSink<A extends PeakAnnotation> implements PeakSink<A> {

    private final TDoubleArrayList mzList = new TDoubleArrayList();
    private final TDoubleArrayList intensityList = new TDoubleArrayList();

    private int startCalled = 0;
    private int endCalled = 0;

    @Override
    public void processPeak(double mz, double intensity, List<A> annotations) {

        mzList.add(mz);
        intensityList.add(intensity);
    }

    public void end() {

        endCalled++;
    }

    public void start(int size) {

        mzList.ensureCapacity(size);
        intensityList.ensureCapacity(size);

        startCalled++;
    }

    public double[] getMzList() {

        return mzList.toArray();
    }

    public double[] getIntensityList() {

        return intensityList.toArray();
    }

    public int getStartCalled() {

        return startCalled;
    }

    public int getEndCalled() {

        return endCalled;
    }
}

