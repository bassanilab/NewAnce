/*
Copyright (C) SIB - Swiss Institute of Bioinformatics, Lausanne, Switzerland
Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
*/

package newance.psmconverter;

import newance.util.PsmPredicate;

import java.io.File;
import java.util.*;
import java.util.concurrent.CountDownLatch;

/**
 * @author Markus Müller
 */

public class MaxQuantMSMSConverter extends SinglePsmFileConverter {

    public MaxQuantMSMSConverter(File msmsFile, Map<String,List<PeptideSpectrumMatch>> psms, CountDownLatch latch) {

        super(msmsFile, psms, latch);
    }

    public void run() {

        System.out.println("Reading " + psmFile);

        final Map<String, List<PeptideSpectrumMatch>> psmMap = new HashMap<>();

        PsmPredicate psmPredicate = new PsmPredicate(params.getMinCharge(), params.getMaxCharge(), params.getMinPeptideLength(), params.getMaxPeptideLength(), params.getMaxRank(),
                params.getMaxQuantMainScore(), params.getMaxQuantMainScoreMinValue(), PsmPredicate.ScoreOrder.LARGER);

        PeptideSpectrumMatchList peptideSpectrumMatchList = new PeptideSpectrumMatchList(new SpectrumKeyFunctionImpl(), psmPredicate, psmMap);

        MaxQuantPsmReader2 psmReader = new MaxQuantPsmReader2();
        psmReader.parse(psmFile, peptideSpectrumMatchList);

        addPsms(psmMap);

        latch.countDown();
        System.out.println("Finished reading " + psmFile+". Latch count: "+latch.getCount());
    }
}
