package newance.psmcombiner;

import newance.psmconverter.PeptideMatchData;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.BiConsumer;

/**
 * Created by markusmueller on 07.01.19.
 */
public class SpectrumAccumulator implements BiConsumer<String, List<PeptideMatchData>> {

    protected final ConcurrentHashMap<String,Set<String>> proteinPeptideMap;
    protected final ConcurrentHashMap<String,Set<String>> peptideProteinMap;
    protected final ConcurrentHashMap<String,List<PeptideMatchData>> peptidePSMMap;

    protected int totSpectrumCount;

    public SpectrumAccumulator() {
        proteinPeptideMap = new ConcurrentHashMap<>();
        peptideProteinMap = new ConcurrentHashMap<>();
        peptidePSMMap = new ConcurrentHashMap<>();

        totSpectrumCount = 0;
    }

    @Override
    public void accept(String s, List<PeptideMatchData> peptideMatchData) {

        if (peptideMatchData.isEmpty()) return;

        for (PeptideMatchData psm : peptideMatchData) {
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

    public ConcurrentHashMap<String, List<PeptideMatchData>> getPeptidePSMMap() {
        return peptidePSMMap;
    }

    public int getTotSpectrumCount() {
        return totSpectrumCount;
    }
}
