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

import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import newance.psmconverter.PeptideSpectrumMatch;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.BiConsumer;

/**
 * @author Markus Müller
 */

public class CometMaxQuantPsmMerger implements BiConsumer<String, List<PeptideSpectrumMatch>> {

    private final ConcurrentHashMap<String, List<PeptideSpectrumMatch>> maxQuantPsmMap;
    private final ConcurrentHashMap<String, List<PeptideSpectrumMatch>> combinedPsmMap;

    public CometMaxQuantPsmMerger(ConcurrentHashMap<String, List<PeptideSpectrumMatch>> maxQuantPsmMap,
                                  ConcurrentHashMap<String, List<PeptideSpectrumMatch>> combinedPsmMap) {
        this.maxQuantPsmMap = maxQuantPsmMap;
        this.combinedPsmMap = combinedPsmMap;
    }

    @Override
    public void accept(String specID, List<PeptideSpectrumMatch> cometPsms) {

        List<PeptideSpectrumMatch> maxQuantPsms = maxQuantPsmMap.get(specID);

        if (maxQuantPsms!=null) {
            List<PeptideSpectrumMatch> combined = combine(cometPsms,maxQuantPsms);

            if (combined != null) combinedPsmMap.put(specID, combined);
        }

    }

    private List<PeptideSpectrumMatch> combine(List<PeptideSpectrumMatch> cometPsms, List<PeptideSpectrumMatch> maxQuantPsms) {

        List<String> maxQuantPeptides = new ArrayList<>();
        maxQuantPsms.forEach(psm -> maxQuantPeptides.add(psm.getPeptide().toString()));
        List<PeptideSpectrumMatch> combined = null;

        for (PeptideSpectrumMatch cometPsm : cometPsms) {
            String cometPeptide = cometPsm.getPeptide().toString();

            for (int i=0;i<maxQuantPeptides.size();i++) {
                if (cometPeptide.equals(maxQuantPeptides.get(i))) {
                    if (combined==null) combined =Collections.synchronizedList(new ArrayList<>());
                    combined.add(merge(cometPsm, maxQuantPsms.get(i)));
                    break;
                }
            }
        }

        return combined;
    }

    private PeptideSpectrumMatch merge(PeptideSpectrumMatch cometPsm, PeptideSpectrumMatch maxQuantPsm) {

        TObjectDoubleMap<String> scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("xcorr",cometPsm.getScore("xcorr"));
        scoreMap.put("deltacn",cometPsm.getScore("deltacn"));
        scoreMap.put("spscore",cometPsm.getScore("spscore"));
        scoreMap.put("neg_log10_p",cometPsm.getScore("neg_log10_p"));
        scoreMap.put("mass_diff",cometPsm.getScore("mass_diff"));
        scoreMap.put("tot_num_ions",cometPsm.getScore("tot_num_ions"));
        scoreMap.put("matched_num_ions",cometPsm.getScore("matched_num_ions"));
        scoreMap.put("rt",maxQuantPsm.getScore("rt"));
        scoreMap.put("sn",cometPsm.getScore("sn"));
        scoreMap.put("mass",cometPsm.getScore("mass"));
        scoreMap.put("rank",cometPsm.getScore("rank"));
        scoreMap.put("Score",maxQuantPsm.getScore("Score"));
        scoreMap.put("Delta score",maxQuantPsm.getScore("Delta score"));
        scoreMap.put("Mass Error [ppm]",maxQuantPsm.getScore("Mass Error [ppm]"));
        scoreMap.put("Intensity coverage",maxQuantPsm.getScore("Intensity coverage"));
        scoreMap.put("Localization prob",maxQuantPsm.getScore("Localization prob"));

        return new PeptideSpectrumMatch(cometPsm.getSpectrumFile(), cometPsm.getPeptide(), cometPsm.getProteinAcc(),
                scoreMap,cometPsm.getCharge(),cometPsm.getRank(), cometPsm.getRetentionTime(), cometPsm.getScanNr(), cometPsm.getNeutralPrecMass(),
                cometPsm.isDecoy(), cometPsm.isVariant(), cometPsm.getVariantPositions(), cometPsm.getVariantWTAAs());

    }
}
