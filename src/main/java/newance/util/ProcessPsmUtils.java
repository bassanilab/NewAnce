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

package newance.util;

import newance.psmcombiner.Psm2PeptideStringFunction;
import newance.psmconverter.PeptideSpectrumMatch;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

/**
 * @author Markus Müller
 */

public class ProcessPsmUtils {

    public static ConcurrentHashMap<String, List<PeptideSpectrumMatch>>  removeDecoys(ConcurrentHashMap<String, List<PeptideSpectrumMatch>>  psms) {

        ConcurrentHashMap<String, List<PeptideSpectrumMatch>> noDecoyPsms = new ConcurrentHashMap<>();

        for (String specID : psms.keySet()) {
            List<PeptideSpectrumMatch> noDecoy = new ArrayList<>();
            for (PeptideSpectrumMatch psm : psms.get(specID)) {
                if (!psm.isDecoy()) noDecoy.add(psm);
            }

            if (!noDecoy.isEmpty()) noDecoyPsms.put(specID, Collections.synchronizedList(noDecoy));
        }

        return noDecoyPsms;
    }


    public static ConcurrentHashMap<String, List<PeptideSpectrumMatch>>  extractGroupPsms(
            PsmGrouper psmGrouper, ConcurrentHashMap<String, List<PeptideSpectrumMatch>>  psms, String group) {

        ConcurrentHashMap<String, List<PeptideSpectrumMatch>> filteredPsms = new ConcurrentHashMap<>();

        for (String specID : psms.keySet()) {
            List<PeptideSpectrumMatch> groupPsm = new ArrayList<>();
            for (PeptideSpectrumMatch psm : psms.get(specID)) {
                if (psmGrouper.apply(specID,psm).equals(group)) groupPsm.add(psm);
            }

            if (!groupPsm.isEmpty()) filteredPsms.put(specID, Collections.synchronizedList(groupPsm));
        }

        return filteredPsms;
    }

    public static int countPsms(ConcurrentHashMap<String, List<PeptideSpectrumMatch>>  psms) {

        int cnt = 0;

        for (String specID : psms.keySet()) {
            cnt += psms.get(specID).size();
        }

        return cnt;
    }

    public static int countUniquePeptides(ConcurrentHashMap<String, List<PeptideSpectrumMatch>>  psms) {

        final Psm2PeptideStringFunction stringFunction = new Psm2PeptideStringFunction(Psm2PeptideStringFunction.StringMode.SEQUENCE);
        final Set<String> peptides = new HashSet<>();

        for (String specID : psms.keySet()) {
            peptides.addAll(stringFunction.apply(specID,psms.get(specID)));
        }

        return peptides.size();
    }


}
