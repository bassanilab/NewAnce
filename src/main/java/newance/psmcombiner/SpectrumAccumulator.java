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

import newance.psmconverter.PeptideSpectrumMatch;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.BiConsumer;

/**
 * @author Markus Müller
 */

public class SpectrumAccumulator implements BiConsumer<String, List<PeptideSpectrumMatch>> {

    protected final ConcurrentHashMap<String,Set<String>> proteinPeptideMap;
    protected final ConcurrentHashMap<String,Set<String>> peptideProteinMap;
    protected final ConcurrentHashMap<String,List<PeptideSpectrumMatch>> peptidePSMMap;

    protected int totSpectrumCount;

    public SpectrumAccumulator() {
        proteinPeptideMap = new ConcurrentHashMap<>();
        peptideProteinMap = new ConcurrentHashMap<>();
        peptidePSMMap = new ConcurrentHashMap<>();

        totSpectrumCount = 0;
    }

    @Override
    public void accept(String s, List<PeptideSpectrumMatch> peptideSpectrumMatchData) {

        if (peptideSpectrumMatchData.isEmpty()) return;

        for (PeptideSpectrumMatch psm : peptideSpectrumMatchData) {
            String peptideSeq = psm.getPeptide().toSymbolString();

            peptideProteinMap.putIfAbsent(peptideSeq,Collections.synchronizedSet(new HashSet<>()));
            peptidePSMMap.putIfAbsent(peptideSeq,Collections.synchronizedList(new ArrayList<>()));

            peptidePSMMap.get(peptideSeq).add(psm);

            for (String protein : psm.getProteinIDs()) {
                proteinPeptideMap.putIfAbsent(protein,Collections.synchronizedSet(new HashSet<>()));

                proteinPeptideMap.get(protein).add(peptideSeq);
                peptideProteinMap.get(peptideSeq).add(protein);
            }

            totSpectrumCount++;
        }
    }

    public ConcurrentHashMap<String, Set<String>> getProteinPeptideMap() {
        return proteinPeptideMap;
    }

    public ConcurrentHashMap<String, Set<String>> getPeptideProteinMap() {
        return peptideProteinMap;
    }

    public ConcurrentHashMap<String, List<PeptideSpectrumMatch>> getPeptidePSMMap() {
        return peptidePSMMap;
    }

    public int getTotSpectrumCount() {
        return totSpectrumCount;
    }
}
