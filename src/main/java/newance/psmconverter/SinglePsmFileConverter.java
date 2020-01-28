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

import java.io.File;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CountDownLatch;

/**
 * @author Markus Müller
 */

public abstract class SinglePsmFileConverter implements Runnable {

    protected final File psmFile;
    protected final NewAnceParams params;
    protected final Map<String,List<PeptideSpectrumMatch>> psms;
    protected final CountDownLatch latch;

    public SinglePsmFileConverter(File psmFile, Map<String,List<PeptideSpectrumMatch>> psms, CountDownLatch latch) {

        this.params = NewAnceParams.getInstance();
        this.psmFile = psmFile;
        this.psms = psms;
        this.latch = latch;
    }

    public abstract void run();

    protected void addPsms(Map<String, List<PeptideSpectrumMatch>> psmMap) {

        for (String spectrumID : psmMap.keySet()) {

            List<PeptideSpectrumMatch> matches = psmMap.get(spectrumID);

            if (psms.containsKey(spectrumID)) {
                psms.get(spectrumID).addAll(matches);
            } else {
                psms.put(spectrumID,Collections.synchronizedList(matches));
            }
        }
    }
}
