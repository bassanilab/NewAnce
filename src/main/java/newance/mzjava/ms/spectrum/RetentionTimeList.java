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
public class RetentionTimeList extends ArrayList<RetentionTime> implements Copyable<RetentionTimeList> {

    public RetentionTimeList() {
    }

    public RetentionTimeList(RetentionTimeList retentionTimes) {

        super(retentionTimes.size());

        for(RetentionTime src : retentionTimes) {

            add(src.copy());
        }
    }

    public RetentionTimeList(RetentionTime... retentionTimes) {

        Collections.addAll(this, retentionTimes);
    }

    @Override
    public RetentionTimeList copy() {

        return new RetentionTimeList(this);
    }

    public void add(double time, TimeUnit unit) {

        add(new RetentionTimeDiscrete(time, unit));
    }

    public void add(double startTime, double endTime, TimeUnit unit) {

        add(new RetentionTimeInterval(startTime, endTime, unit));
    }

    public RetentionTime getFirst(){

        return get(0);
    }

    public RetentionTime getLast() {

        return get(size() - 1);
    }

    @Override
    public boolean contains(Object o) {

        if(!(o instanceof RetentionTime)) return false;

        RetentionTime rt = (RetentionTime)o;
        for(RetentionTime retentionTime : this) {

            if(retentionTime.contains(rt)) return true;
        }
        return false;
    }
}
