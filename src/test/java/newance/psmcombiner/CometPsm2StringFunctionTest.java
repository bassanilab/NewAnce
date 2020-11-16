package newance.psmcombiner;

import org.junit.Test;

/**
 * Created by markusmueller on 27.10.20.
 */
public class CometPsm2StringFunctionTest {

    @Test
    public void test_plus_minus_sign() {
        String valStr = "+0.0001";

        Double d = Double.valueOf(valStr);

        System.out.println(d);
    }
}
