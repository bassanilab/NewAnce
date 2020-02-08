/**
 * Copyright (C) 2019, SIB/LICR. All rights reserved
 *
 * SIB, Swiss Institute of Bioinformatics
 * Ludwig Institute for Cancer Research (LICR)
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer. Redistributions in binary
 * form must reproduce the above copyright notice, this list of conditions and
 * the following disclaimer in the documentation and/or other materials provided
 * with the distribution. Neither the name of the SIB/LICR nor the names of
 * its contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL SIB/LICR BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
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

    @Test
    public void readPriorHistoTest() {
        String filename = getClass().getResource("prior_histo_Z2.txt").getFile();
        CometScoreHistogram histogram = CometScoreHistogram.read(new File(filename));
    }

}
