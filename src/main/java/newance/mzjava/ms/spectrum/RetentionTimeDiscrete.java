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
 * The retention time is the time it takes for a particular
 * analyte to pass through the liquid chromatography system (from the column
 * inlet to the detector) under set conditions.
 *
 * Retention times are stored in seconds
 *
 * @author masselot
 * @author nikitin
 * @author Oliver Horlacher
 *
 * @version 1.0
 */
public class RetentionTimeDiscrete implements RetentionTime {

    private final double time;

    /**
     * Constructor for a retention time.  The <code>time</code> is
     * converted to seconds.
     *
     * @param time the time for this retention time
     * @param unit the unit for the <code>time</code>
     */
    public RetentionTimeDiscrete(final double time, final TimeUnit unit) {

        this.time = unit.convert(time, TimeUnit.SECOND);
    }

    /**
     * Copy constructor
     *
     * @param src the RetentionTime to copy
     */
    public RetentionTimeDiscrete(RetentionTimeDiscrete src) {

        this.time = src.time;
    }

    @Override
    public double getTime() {

        return time;
    }

    @Override
    public double getMinRetentionTime() {
        
        return time;
    }

    @Override
    public double getMaxRetentionTime() {
        
        return time;
    }


    @Override
    public RetentionTimeDiscrete copy() {

        return new RetentionTimeDiscrete(this);
    }

    @Override
    public boolean contains(RetentionTime retentionTime) {

        return time == retentionTime.getMinRetentionTime() && time == retentionTime.getMaxRetentionTime();
    }

    @Override
    public String toString() {

        return Double.toString(time) + 's';
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        RetentionTimeDiscrete that = (RetentionTimeDiscrete) o;

        return Double.compare(that.time, time) == 0;

    }

    @Override
    public int hashCode() {

        long temp = time != +0.0d ? Double.doubleToLongBits(time) : 0L;
        return (int) (temp ^ (temp >>> 32));
    }
}
