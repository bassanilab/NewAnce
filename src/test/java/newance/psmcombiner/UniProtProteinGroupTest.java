package newance.psmcombiner;

import org.junit.Assert;
import org.junit.Test;

/**
 * Copyright (C) 2019
 * @author Markus Müller
 * @Institutions: SIB, Swiss Institute of Bioinformatics; Ludwig Institute for Cancer Research
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
