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
