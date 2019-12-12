package newance.psmcombiner;

import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import newance.psmconverter.PeptideMatchData;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.BiConsumer;

/**
 * Copyright (C) 2019
 * @author Markus Müller
 * @Institutions: SIB, Swiss Institute of Bioinformatics; Ludwig Institute for Cancer Research
 */

public class CometMaxQuantPsmCombiner implements BiConsumer<String, List<PeptideMatchData>> {

    private final ConcurrentHashMap<String, List<PeptideMatchData>> maxQuantPsmMap;
    private final ConcurrentHashMap<String, List<PeptideMatchData>> combinedPsmMap;

    public CometMaxQuantPsmCombiner(ConcurrentHashMap<String, List<PeptideMatchData>> maxQuantPsmMap,
                                    ConcurrentHashMap<String, List<PeptideMatchData>> combinedPsmMap) {
        this.maxQuantPsmMap = maxQuantPsmMap;
        this.combinedPsmMap = combinedPsmMap;
    }

    @Override
    public void accept(String specID, List<PeptideMatchData> cometPsms) {

        List<PeptideMatchData> maxQuantPsms = maxQuantPsmMap.get(specID);

        if (maxQuantPsms!=null) {
            List<PeptideMatchData> combined = combine(cometPsms,maxQuantPsms);

            if (combined != null) combinedPsmMap.put(specID, combined);
        }

    }

    private List<PeptideMatchData> combine(List<PeptideMatchData> cometPsms, List<PeptideMatchData> maxQuantPsms) {

        List<String> maxQuantPeptides = new ArrayList<>();
        maxQuantPsms.forEach(psm -> maxQuantPeptides.add(psm.getPeptide().toString()));
        List<PeptideMatchData> combined = null;

        for (PeptideMatchData cometPsm : cometPsms) {
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

    private PeptideMatchData merge(PeptideMatchData cometPsm, PeptideMatchData maxQuantPsm) {
        TObjectDoubleMap<String> scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("xcorr",cometPsm.getScore("xcorr"));
        scoreMap.put("deltacn",cometPsm.getScore("deltacn"));
        scoreMap.put("spscore",cometPsm.getScore("spscore"));
        scoreMap.put("neg_log10_p",cometPsm.getScore("neg_log10_p"));
        scoreMap.put("rt",maxQuantPsm.getScore("rt"));
        scoreMap.put("rank",cometPsm.getScore("rank"));
        scoreMap.put("score",maxQuantPsm.getScore("score"));
        scoreMap.put("delta-score",maxQuantPsm.getScore("delta-score"));
        scoreMap.put("mass-error-ppm",maxQuantPsm.getScore("mass-error-ppm"));
        scoreMap.put("int-coverage",maxQuantPsm.getScore("int-coverage"));
        scoreMap.put("mod-pos-prob",maxQuantPsm.getScore("mod-pos-prob"));

        return new PeptideMatchData(cometPsm.getPeptide(), cometPsm.getProteinAcc(),scoreMap,cometPsm.getCharge(),cometPsm.isDecoy());

    }
}
