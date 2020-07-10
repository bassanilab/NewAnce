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

import java.util.Set;

/**
 * Created by markusmueller on 21.02.20.
 */
public class NewAnceParamsTest {

    @Test
    public void testGetSetValue() {

        String setStr = "[[a,b,c,d,e]]";

        Set<String> set = NewAnceParams.getSetValue(setStr);

        Assert.assertEquals(5,set.size());
        Assert.assertTrue(set.contains("a"));
        Assert.assertTrue(set.contains("b"));
        Assert.assertTrue(set.contains("c"));
        Assert.assertTrue(set.contains("d"));
        Assert.assertTrue(set.contains("e"));

        setStr = "[a,b,c,d,e]";

        set = NewAnceParams.getSetValue(setStr);

        Assert.assertEquals(5,set.size());
        Assert.assertTrue(set.contains("a"));
        Assert.assertTrue(set.contains("b"));
        Assert.assertTrue(set.contains("c"));
        Assert.assertTrue(set.contains("d"));
        Assert.assertTrue(set.contains("e"));


        setStr = "[a,b,c,d,e[1]]";

        set = NewAnceParams.getSetValue(setStr);

        Assert.assertEquals(5,set.size());
        Assert.assertTrue(set.contains("a"));
        Assert.assertTrue(set.contains("b"));
        Assert.assertTrue(set.contains("c"));
        Assert.assertTrue(set.contains("d"));
        Assert.assertTrue(set.contains("e[1]"));


        setStr = "[[1]a,b,c,d,e]";

        set = NewAnceParams.getSetValue(setStr);

        Assert.assertEquals(5,set.size());
        Assert.assertTrue(set.contains("[1]a"));
        Assert.assertTrue(set.contains("b"));
        Assert.assertTrue(set.contains("c"));
        Assert.assertTrue(set.contains("d"));
        Assert.assertTrue(set.contains("e"));


        setStr = "a,b,c,d,e";

        set = NewAnceParams.getSetValue(setStr);

        Assert.assertEquals(5,set.size());
        Assert.assertTrue(set.contains("a"));
        Assert.assertTrue(set.contains("b"));
        Assert.assertTrue(set.contains("c"));
        Assert.assertTrue(set.contains("d"));
        Assert.assertTrue(set.contains("e"));
    }
}
