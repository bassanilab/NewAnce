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

            for (String protein : psm.getProteinAcc()) {
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
