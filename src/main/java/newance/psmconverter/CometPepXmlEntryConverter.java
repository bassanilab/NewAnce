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

import newance.util.PsmPredicate;
import org.expasy.mzjava.core.ms.AbsoluteTolerance;
import org.expasy.mzjava.proteomics.io.ms.ident.PepXmlReader;
import org.expasy.mzjava.proteomics.mol.modification.Modification;
import org.expasy.mzjava.proteomics.ms.ident.ModListModMatchResolver;
import newance.util.NewAnceParams;

import java.io.File;
import java.util.*;
import java.util.concurrent.CountDownLatch;
import java.util.logging.Logger;
import java.util.stream.Collectors;

/**
 * @author Markus Müller
 */

public class CometPepXmlEntryConverter implements Runnable {

    private static final Logger LOGGER = Logger.getLogger(CometPepXmlEntryConverter.class.getName());

    protected final File pepXmlFile;
    protected final NewAnceParams params;
    protected final Map<String,List<PeptideMatchData>> psms;
    protected final CountDownLatch latch;

    public CometPepXmlEntryConverter(File pepXmlFile, Map<String,List<PeptideMatchData>> psms, CountDownLatch latch) {

        this.params = NewAnceParams.getInstance();
        this.pepXmlFile = pepXmlFile;
        this.psms = psms;
        this.latch = latch;
    }

    @Override
    public void run() {

        System.out.println("Reading " + pepXmlFile);

        final Map<String, List<PeptideMatchData>> psmMap = new HashMap<>();

        PsmPredicate psmPredicate = new PsmPredicate(params.getMinCharge(), params.getMaxCharge(), params.getMinPeptideLength(), params.getMaxPeptideLength(), params.getMaxRank(),
                "xcorr", 1.0f, PsmPredicate.ScoreOrder.LARGER);

        PSMReaderCallbackImpl callback = new PSMReaderCallbackImpl(new SpectrumKeyFunctionImpl(), psmPredicate, psmMap);

        Collection<Modification> modifications = params.getModifications();
        ModListModMatchResolver modMatchResolver = new ModListModMatchResolver(new AbsoluteTolerance(params.getModifMatchMassTol()), modifications);

        CometPepXMLReader psmReader = new CometPepXMLReader(PepXmlReader.ModMassStorage.AA_MASS_PLUS_MOD_MASS, true, modMatchResolver);
        psmReader.parse(pepXmlFile, callback);

        addPsms(psmMap);

        if (latch!=null) latch.countDown();
        System.out.println("Finished reading " + pepXmlFile+". Latch count: "+((latch!=null)?latch.getCount():-1));
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
