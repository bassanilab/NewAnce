/**
 * Copyright (C) 2019, SIB/LICR. All rights reserved
 *
 * SIB, Swiss Institute of Bioinformatics
 * Ludwig Institute for Cancer Research (LICR)
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer. Redistributions in binary
 * form must reproduce the above copyright notice, this list of conditions and
 * the following disclaimer in the documentation and/or other materials provided
 * with the distribution. Neither the name of the SIB/LICR nor the names of
 * its contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL SIB/LICR BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

package newance.util;

import org.junit.Assert;
import org.junit.Test;

/**
 * @author Markus Müller
 */

public class RunTime2StringTest {

    @Test
    public void getTimeDiffStringTest() {

        Assert.assertEquals("1ms",RunTime2String.getTimeDiffString(1));
        Assert.assertEquals("111ms",RunTime2String.getTimeDiffString(111));
        Assert.assertEquals("1sec, 1ms",RunTime2String.getTimeDiffString(1001));
        Assert.assertEquals("1sec, 111ms",RunTime2String.getTimeDiffString(1111));
        Assert.assertEquals("1min, 1sec, 1ms",RunTime2String.getTimeDiffString(61001));
        Assert.assertEquals("1min, 1sec, 111ms",RunTime2String.getTimeDiffString(61111));
        Assert.assertEquals("24min, 41sec, 111ms",RunTime2String.getTimeDiffString(1481111));
        Assert.assertEquals("1day, 0h, 0min, 0sec, 0ms",RunTime2String.getTimeDiffString(24*60*60*1000));
    }

}
