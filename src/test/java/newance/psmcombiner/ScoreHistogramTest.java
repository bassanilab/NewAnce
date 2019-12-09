package newance.psmcombiner;

import org.junit.Assert;
import org.junit.Test;

import java.util.*;

/**
 * Created by markusmueller on 15.11.19.
 */
public class ScoreHistogramTest {

    @Test
    public void sortIndexesTest() {
        List<Float> array = Arrays.asList(new Float[]{5f,2f,1f,4f});
        List<Integer> indexes = Arrays.asList(new Integer[]{2,1,3,0});
        List<Integer> sorted = ScoreHistogram.sortIndexes(array, false);

        for (int i=0;i<indexes.size();i++) Assert.assertEquals(indexes.get(i),sorted.get(i));

        indexes = Arrays.asList(new Integer[]{0,3,1,2});
        sorted = ScoreHistogram.sortIndexes(array, true);

        for (int i=0;i<indexes.size();i++) Assert.assertEquals(indexes.get(i),sorted.get(i));


        array = new ArrayList<>();
        for (int i=0;i<100;i++) array.add((float)Math.random());
        sorted = ScoreHistogram.sortIndexes(array, false);
        for (int i=1;i<100;i++) {
            Assert.assertTrue(array.get(sorted.get(i))-array.get(sorted.get(i-1))>0.0);
        }
    }

}
