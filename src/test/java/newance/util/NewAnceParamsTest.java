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

        Set<String> set = NewAnceParams.getSetValue("bla",setStr);

        Assert.assertEquals(5,set.size());
        Assert.assertTrue(set.contains("a"));
        Assert.assertTrue(set.contains("b"));
        Assert.assertTrue(set.contains("c"));
        Assert.assertTrue(set.contains("d"));
        Assert.assertTrue(set.contains("e"));

        setStr = "[a,b,c,d,e]";

        set = NewAnceParams.getSetValue("bla",setStr);

        Assert.assertEquals(5,set.size());
        Assert.assertTrue(set.contains("a"));
        Assert.assertTrue(set.contains("b"));
        Assert.assertTrue(set.contains("c"));
        Assert.assertTrue(set.contains("d"));
        Assert.assertTrue(set.contains("e"));


        setStr = "a,b,c,d,e";

        set = NewAnceParams.getSetValue("bla",setStr);

        Assert.assertEquals(5,set.size());
        Assert.assertTrue(set.contains("a"));
        Assert.assertTrue(set.contains("b"));
        Assert.assertTrue(set.contains("c"));
        Assert.assertTrue(set.contains("d"));
        Assert.assertTrue(set.contains("e"));
    }
}
