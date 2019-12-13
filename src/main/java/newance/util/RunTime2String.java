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

import java.text.SimpleDateFormat;

/**
 * @author Markus Müller
 */

public class RunTime2String {

    private static final SimpleDateFormat TIME_FORMATTER = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss.SSS");

    public static String getTimeDiffString(long miliSec) {
        int days = 0;
        int hours = 0;
        int mins = 0;
        int secs = 0;

        String res = "";
        secs = (int) Math.floor(miliSec/1000.0);
        res = (miliSec-secs*1000)+"ms";
        if (secs>0) {
            mins = (int) Math.floor(secs/60.0);
            res = (secs-mins*60)+"sec, "+res;
            if (mins>0) {
                hours = (int) Math.floor(mins/60.0);
                res = (mins-hours*60)+"min, "+res;
                if (hours>0) {
                    days = (int) Math.floor(hours/24.0);
                    res = (hours-days*24)+"h, "+res;
                }
                if (days>0) {
                    res = days+"day, "+res;
                }
            }
        }

        return res;

    }

    public static String getRuntimeString(Runtime runtime) {
        return ".  RAM: "+(runtime.totalMemory()-runtime.freeMemory())+" bytes"+", Time: "+TIME_FORMATTER.format(System.currentTimeMillis());

    }
}
