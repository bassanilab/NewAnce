/*
Copyright (C) SIB - Swiss Institute of Bioinformatics, Lausanne, Switzerland
Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
*/

package newance.psmcombiner;

import org.expasy.mzjava.stats.Histogram;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.util.*;

/**
 * @author Markus Müller
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
