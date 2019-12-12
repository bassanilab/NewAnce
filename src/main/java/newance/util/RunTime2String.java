package newance.util;

import java.text.SimpleDateFormat;

/**
 * Copyright (C) 2019
 * @author Markus Müller
 * @Institutions: SIB, Swiss Institute of Bioinformatics; Ludwig Institute for Cancer Research
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
