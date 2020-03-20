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
 * A retention time that stores a time interval.  The getTime method returns the
 * time in the middle of the interval.
 *
 * All times are stored in seconds
 *
 * @author Oliver Horlacher
 * @version 1.0
 */
public class RetentionTimeInterval implements RetentionTime {

    private final double startTime, endTime;

    /**
     * Constructor for a RetentionTimeInterval from startTime to endTime.  Both times
     * are converted to seconds
     * 
     * @param startTime the start time
     * @param endTime the end time
     * @param unit the time unit for <code>startTime</code> and <code>endTime</code>
     */
    public RetentionTimeInterval(double startTime, double endTime, TimeUnit unit) {

        if(startTime > endTime) throw new IllegalArgumentException("Start time has to be <= than the end time start time was " + startTime + " and end time was " + endTime);

        this.startTime = unit.convert(startTime, TimeUnit.SECOND);
        this.endTime = unit.convert(endTime, TimeUnit.SECOND);
    }

    /**
     * Copy constructor
     *
     * @param src the RetentionTime to copy
     */
    public RetentionTimeInterval(RetentionTimeInterval src) {
        
        startTime = src.startTime;
        endTime = src.endTime;
    }

    /**
     * Returns the time in the middle of the interval in seconds
     *
     * @return the time in the middle of the interval in seconds
     */
    @Override
    public double getTime() {

        return startTime + (endTime - startTime)/2;
    }

    @Override
    public double getMinRetentionTime() {
        
        return startTime;
    }

    @Override
    public double getMaxRetentionTime() {
        
        return endTime;
    }

    @Override
    public RetentionTimeInterval copy() {

        return new RetentionTimeInterval(this);
    }

    @Override
    public boolean contains(RetentionTime retentionTime) {

        return retentionTime.getMinRetentionTime() >= startTime && retentionTime.getMaxRetentionTime() <= endTime;
    }

    @Override
    public boolean equals(Object o) {
        
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        RetentionTimeInterval that = (RetentionTimeInterval) o;

        return Double.compare(that.endTime, endTime) == 0 && 
                Double.compare(that.startTime, startTime) == 0;

    }

    @Override
    public int hashCode() {
        
        int result;
        long temp;
        temp = startTime != +0.0d ? Double.doubleToLongBits(startTime) : 0L;
        result = (int) (temp ^ (temp >>> 32));
        temp = endTime != +0.0d ? Double.doubleToLongBits(endTime) : 0L;
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        return result;
    }

    @Override
    public String toString() {

        return String.valueOf('[') + startTime + "s - " + endTime + "s]";
    }
}
