package newance.mzjava.ms.spectrum;

import org.junit.Assert;
import org.junit.Test;

/**
 * @author markus
 */
public class PepLibPeakAnnotationTest {
    @Test
    public void testEquals() throws Exception {
        PepLibPeakAnnotation annotation1 = new PepLibPeakAnnotation(3, 1.0, 100.0);
        PepLibPeakAnnotation annotation2 = new PepLibPeakAnnotation(3, 1.0, 100.0);
        PepLibPeakAnnotation annotation3 = new PepLibPeakAnnotation(1, 1.0, 100.0);

        Assert.assertNotSame(annotation1, annotation2);

        Assert.assertTrue(annotation1.equals(annotation2));
        Assert.assertTrue(annotation2.equals(annotation1));

        Assert.assertFalse(annotation1.equals(annotation3));
        Assert.assertFalse(annotation3.equals(annotation1));
    }

    @Test
    public void testHashCode() throws Exception {

        PepLibPeakAnnotation annotation1 = new PepLibPeakAnnotation(3, 1.0, 100.0);
        PepLibPeakAnnotation annotation2 = new PepLibPeakAnnotation(3, 1.0, 100.0);
        PepLibPeakAnnotation annotation3 = new PepLibPeakAnnotation(1, 1.0, 100.0);

        int hash1 = annotation1.hashCode();
        int hash2 = annotation2.hashCode();
        int hash3 = annotation3.hashCode();

        Assert.assertEquals(hash1, hash2);
        Assert.assertTrue(hash1 != hash3);
    }

    @Test
    public void testCopy() throws Exception {

        PepLibPeakAnnotation annotation = new PepLibPeakAnnotation(3, 1.0, 100.0);

        PepLibPeakAnnotation copyAnnotation = annotation.copy();
        PepLibPeakAnnotation copyConstructedAnnotation = new PepLibPeakAnnotation(annotation);

        Assert.assertNotSame(copyAnnotation, copyConstructedAnnotation);
        Assert.assertTrue(copyAnnotation.equals(copyConstructedAnnotation));
        Assert.assertTrue(copyConstructedAnnotation.equals(copyAnnotation));
    }
}
