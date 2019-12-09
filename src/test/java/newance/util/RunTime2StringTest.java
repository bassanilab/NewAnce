package newance.util;

import org.junit.Assert;
import org.junit.Test;

/**
 * Created by markusmueller on 06.11.19.
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
