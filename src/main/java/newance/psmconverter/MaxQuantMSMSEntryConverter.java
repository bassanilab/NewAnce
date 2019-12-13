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
package newance.psmconverter;

import newance.util.NewAnceParams;
import newance.util.PsmPredicate;

import java.io.File;
import java.util.*;
import java.util.concurrent.CountDownLatch;
import java.util.logging.Logger;

/**
 * @author Markus Müller
 */

public class MaxQuantMSMSEntryConverter  implements Runnable {

    private static final Logger LOGGER = Logger.getLogger(MaxQuantMSMSEntryConverter.class.getName());

    protected final File msmsFile;
    protected final NewAnceParams params;
    protected final Map<String,List<PeptideMatchData>> psms;
    protected final CountDownLatch latch;

    public MaxQuantMSMSEntryConverter(File msmsFile, Map<String,List<PeptideMatchData>> psms, CountDownLatch latch) {

        this.params = NewAnceParams.getInstance();
        this.msmsFile = msmsFile;
        this.psms = psms;
        this.latch = latch;
    }

    public void run() {

        System.out.println("Reading " + msmsFile);

        final Map<String, List<PeptideMatchData>> psmMap = new HashMap<>();

        PsmPredicate psmPredicate = new PsmPredicate(params.getMinCharge(), params.getMaxCharge(), params.getMinPeptideLength(), params.getMaxPeptideLength(), params.getMaxRank(),
                "score", 10.0f, PsmPredicate.ScoreOrder.LARGER);

        PSMReaderCallbackImpl callback = new PSMReaderCallbackImpl(new SpectrumKeyFunctionImpl(), psmPredicate, psmMap);

        MaxQuantPsmReader2 psmReader = new MaxQuantPsmReader2();
        psmReader.parse(msmsFile, callback);

        addPsms(psmMap);

        latch.countDown();
        System.out.println("Finished reading " + msmsFile+". Latch count: "+latch.getCount());
    }

    protected void addPsms(Map<String, List<PeptideMatchData>> psmMap) {

        for (String spectrumID : psmMap.keySet()) {

            List<PeptideMatchData> matches = psmMap.get(spectrumID);

            if (psms.containsKey(spectrumID)) {
                psms.get(spectrumID).addAll(matches);
            } else {
                psms.put(spectrumID,Collections.synchronizedList(matches));
            }
        }
    }
}
