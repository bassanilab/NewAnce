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
package newance.mzjava.ms.io.mgf;

import java.io.IOException;

/**
 * A reader that can be iterated over
 *
 * @author nikitin
 * @author Oliver Horlacher
 * @version 1.0
 */
public interface IterativeReader<T> {

    /**
     * The empty IterativeReader.
     */
    static final IterativeReader<Object> EMPTY_READER = new IterativeReader<Object>() {
        @Override
        public boolean hasNext() {

            return false;
        }

        @Override
        public Object next() throws IOException {

            return null;
        }

        @Override
        public void close() throws IOException {

            // Do nothing because this reader is empty
        }
    };

    /**
     * Returns {@code true} if the reader has more elements.
     * (In other words, returns {@code true} if {@link #next} would
     * return an element rather than throwing an exception.)
     *
     * @return {@code true} if the reader has more elements
     */
    boolean hasNext();

    /**
     * Returns the next element in the reader.
     *
     * @return the next element in the iteration
     * @throws IOException if there is an io exception
     */
    T next() throws IOException;

    /**
     * Closes the reader and releases any system resources associated with
     * it.  Once the reader has been closed, further hasNext() or next()
     * invocations will throw an IOException.
     * Closing a previously closed reader has no effect.
     *
     * @exception  IOException  If an I/O error occurs
     */
    void close() throws IOException;
}
