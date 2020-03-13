/*
Copyright (C) SIB - Swiss Institute of Bioinformatics, Lausanne, Switzerland
Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
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
