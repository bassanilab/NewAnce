/*
Copyright (C) SIB - Swiss Institute of Bioinformatics, Lausanne, Switzerland
Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
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
