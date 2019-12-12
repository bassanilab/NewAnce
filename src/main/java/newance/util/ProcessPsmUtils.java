package newance.util;

import newance.psmcombiner.Psm2PeptideStringFunction;
import newance.psmconverter.PeptideMatchData;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Copyright (C) 2019
 * @author Markus Müller
 * @Institutions: SIB, Swiss Institute of Bioinformatics; Ludwig Institute for Cancer Research
 */

public class ProcessPsmUtils {

    public static ConcurrentHashMap<String, List<PeptideMatchData>>  removeDecoys(ConcurrentHashMap<String, List<PeptideMatchData>>  psms) {

        ConcurrentHashMap<String, List<PeptideMatchData>> noDecoyPsms = new ConcurrentHashMap<>();

        for (String specID : psms.keySet()) {
            List<PeptideMatchData> noDecoy = new ArrayList<>();
            for (PeptideMatchData psm : psms.get(specID)) {
                if (!psm.isDecoy()) noDecoy.add(psm);
            }

            if (!noDecoy.isEmpty()) noDecoyPsms.put(specID, Collections.synchronizedList(noDecoy));
        }

        return noDecoyPsms;
    }


    public static ConcurrentHashMap<String, List<PeptideMatchData>>  extractGroupPsms(
            PsmGrouper psmGrouper, ConcurrentHashMap<String, List<PeptideMatchData>>  psms, String group) {

        ConcurrentHashMap<String, List<PeptideMatchData>> filteredPsms = new ConcurrentHashMap<>();

        for (String specID : psms.keySet()) {
            List<PeptideMatchData> groupPsm = new ArrayList<>();
            for (PeptideMatchData psm : psms.get(specID)) {
                if (psmGrouper.apply(specID,psm).equals(group)) groupPsm.add(psm);
            }

            if (!groupPsm.isEmpty()) filteredPsms.put(specID, Collections.synchronizedList(groupPsm));
        }

        return filteredPsms;
    }

    public static int countPsms(ConcurrentHashMap<String, List<PeptideMatchData>>  psms) {

        int cnt = 0;

        for (String specID : psms.keySet()) {
            cnt += psms.get(specID).size();
        }

        return cnt;
    }

    public static int countUniquePeptides(ConcurrentHashMap<String, List<PeptideMatchData>>  psms) {

        final Psm2PeptideStringFunction stringFunction = new Psm2PeptideStringFunction(Psm2PeptideStringFunction.StringMode.SEQUENCE);
        final Set<String> peptides = new HashSet<>();

        psms.forEach(10000,stringFunction, s -> peptides.addAll(s));

        return peptides.size();
    }


}
