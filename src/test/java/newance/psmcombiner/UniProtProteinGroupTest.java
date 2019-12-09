package newance.psmcombiner;

import org.junit.Assert;
import org.junit.Test;

/**
 * Created by markusmueller on 17.06.18.
 */
public class UniProtProteinGroupTest {
    @Test
    public void test_reverse() {
        String seq = "madamImadam";
        Assert.assertEquals(seq,UniProtProteinGrouper.reverse(seq));

        seq = "0123456789";
        Assert.assertEquals("9876543210",UniProtProteinGrouper.reverse(seq));
    }
}
