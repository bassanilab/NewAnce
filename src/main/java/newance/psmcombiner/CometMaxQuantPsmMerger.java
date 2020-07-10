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

        cometPsm.addScores(maxQuantPsm.getScoreMap());

        return cometPsm;

    }
}
