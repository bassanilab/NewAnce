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
package newance.mzjava.mol;

/**
 * {@code Counter} is a simple counter and a better alternative than {@code
 * Integer} as a new instance is created each time you want to set the value.
 *
 * @author nikitin
 * @author Oliver Horlacher
 *
 * @version 1.0
 *
 */
public class Counter {

    private int count;

    public Counter(){

        count = 0;
    }

    public Counter(int start) {

        count = start;
    }

    /**
     * Reset the count to 0
     */
    public void reset() {

        count = 0;
    }

    @Override
    public boolean equals(Object o) {

        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Counter counter = (Counter) o;

        return count == counter.count;

    }

    @Override
    public int hashCode() {

        return count;
    }

    /**
     * Increment counter by 1.
     */
    public void increment() {

        count++;
    }

    /**
     * Increment counter by {@code inc}.
     */
    public void increment(int inc) {

        count += inc;
    }

    /**
     * Decrement counter by 1.
     */
    public void decrement() {

        count--;
    }

    /**
     * Decrement counter by {@code dec}.
     */
    public void decrement(int dec) {

        count -= dec;
    }

    /**
     * @return counter value.
     */
    public int getCount() {

        return count;
    }

    public String toString() {

        return Integer.toString(count);
    }
}

