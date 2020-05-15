/*
Copyright (C) SIB - Swiss Institute of Bioinformatics, Lausanne, Switzerland
Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
*/

package newance.proteinmatch;


import newance.psmcombiner.SpectrumAccumulator;
import newance.psmconverter.PeptideSpectrumMatch;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

/**
 * @author Markus Müller
 */

public class OccamRazorSpectrumCounter {

    class MapValueComparator implements Comparator<Map.Entry<String,Double>> {

        public MapValueComparator() {
        }

        public int compare(Map.Entry<String,Double> entry1, Map.Entry<String,Double> entry2) {

           int val = entry2.getValue().compareTo(entry1.getValue());

           if (val==0) {
               if (entry1.getKey().startsWith("sp|") && !entry2.getKey().startsWith("sp|")) return -1;
               else if (!entry1.getKey().startsWith("sp|") && entry2.getKey().startsWith("sp|")) return 1;
               else if (entry1.getKey().startsWith("tr|") && !entry2.getKey().startsWith("tr|")) return -1;
               else if (!entry1.getKey().startsWith("tr|") && entry2.getKey().startsWith("tr|")) return +1;
           }

           return val;
        }
    }

    class ProteinComparator implements Comparator<String> {

        protected final Map<String,Double> sortedSpectrumCountsMap;

        public ProteinComparator(Map<String,Double> sortedSpectrumCountsMap) {
            this.sortedSpectrumCountsMap = sortedSpectrumCountsMap;
        }

        public int compare(String protein1, String protein2) {

            int val = sortedSpectrumCountsMap.get(protein2).compareTo(sortedSpectrumCountsMap.get(protein1));

            if (val==0) {
                if (protein1.startsWith("sp|") && !protein2.startsWith("sp|")) return -1;
                else if (!protein1.startsWith("sp|") && protein2.startsWith("sp|")) return 1;
                else if (protein1.startsWith("tr|") && !protein2.startsWith("tr|")) return -1;
                else if (!protein1.startsWith("tr|") && protein2.startsWith("tr|")) return +1;
            }

            return val;
        }
    }

    protected final ConcurrentHashMap<String,Set<String>> proteinPeptideMap;
    protected final ConcurrentHashMap<String,Set<String>> peptideProteinMap;
    protected final ConcurrentHashMap<String,List<PeptideSpectrumMatch>> peptidePSMMap;
    protected final Map<String,Set<String>> proteinPeptideMapAdj;
    protected final UniProtDB uniProtDB; // uniprot db is only required to get sequence length of protein and gene name

    protected final Map<String,Double> sortedSpectrumCountsMap;
    protected final Map<String,Double> sortedSpectrumCountsMapAdj;
    protected final Map<String,Set<String>> razorProteinMap;
    protected final Map<String,String> proteinGroupMap;
    protected final MapValueComparator mapValueComparator;

    public OccamRazorSpectrumCounter(SpectrumAccumulator spectrumAccumulator) {

        proteinPeptideMap = spectrumAccumulator.getProteinPeptideMap();
        peptideProteinMap = spectrumAccumulator.getPeptideProteinMap();
        peptidePSMMap = spectrumAccumulator.getPeptidePSMMap();

        proteinPeptideMapAdj = new HashMap<>();
        copyPeptideProteinMap();

        this.uniProtDB = null;

        mapValueComparator = new MapValueComparator();
        razorProteinMap = new HashMap<>();
        proteinGroupMap = new HashMap<>();
        sortedSpectrumCountsMap = new LinkedHashMap<>();
        sortedSpectrumCountsMapAdj = new LinkedHashMap<>();
    }

    public OccamRazorSpectrumCounter(SpectrumAccumulator spectrumAccumulator, UniProtDB uniProtDB) {

        proteinPeptideMap = spectrumAccumulator.getProteinPeptideMap();
        peptideProteinMap = spectrumAccumulator.getPeptideProteinMap();
        peptidePSMMap = spectrumAccumulator.getPeptidePSMMap();

        proteinPeptideMapAdj = new HashMap<>();
        copyPeptideProteinMap();

        this.uniProtDB = uniProtDB;

        mapValueComparator = new MapValueComparator();
        razorProteinMap = new HashMap<>();
        proteinGroupMap = new HashMap<>();
        sortedSpectrumCountsMap = new LinkedHashMap<>();
        sortedSpectrumCountsMapAdj = new LinkedHashMap<>();
    }

    protected void copyPeptideProteinMap() {

        for (String protein : proteinPeptideMap.keySet()) {
            proteinPeptideMapAdj.put(protein,new HashSet<>(proteinPeptideMap.get(protein)));
        }
    }

    protected void calcSpectrumCounts(Map<String,Set<String>> proteinPeptideMap,Map<String,Double> sortedSpectrumCountsMap) {

        Map<String,Double> spectrumCountsMap = new HashMap<>();

        for (String protein : proteinPeptideMap.keySet()) {
            double count = 0.0;
            for (String peptideSeq : proteinPeptideMap.get(protein)) {
                count += peptidePSMMap.get(peptideSeq).size();
            }
            spectrumCountsMap.put(protein,count);
        }

        List<Map.Entry<String,Double>> tmpList = new ArrayList<>(spectrumCountsMap.entrySet());
        Collections.sort(tmpList, mapValueComparator);

        sortedSpectrumCountsMap.clear();

        for (Map.Entry<String,Double> entry : tmpList) {
            sortedSpectrumCountsMap.put(entry.getKey(),entry.getValue());
        }
    }

    protected void adjustSpectrumCounts() {

        for (String protein : razorProteinMap.keySet()) {

            List<String> linkedProts = new ArrayList<>(razorProteinMap.get(protein));

            if (linkedProts.size()==1) continue;

            Collections.sort(linkedProts, new ProteinComparator(sortedSpectrumCountsMap));

            for (int i=0;i<linkedProts.size();i++) {

                Set<String> peptides = proteinPeptideMapAdj.get(linkedProts.get(i));
                for (int j=i+1;j<linkedProts.size();j++) {
                    proteinPeptideMapAdj.get(linkedProts.get(j)).removeAll(peptides);
                }
            }
        }

    }


    protected void calcSpectrumCounts() {

        calcSpectrumCounts(proteinPeptideMap,sortedSpectrumCountsMap);
        calcProteinGroups();
        adjustSpectrumCounts();
        calcSpectrumCounts(proteinPeptideMapAdj,sortedSpectrumCountsMapAdj);
    }


    protected void calcProteinGroups() {

        Set<String> remainingProteins = new HashSet<>(sortedSpectrumCountsMap.keySet());
        Set<String> proteinGrp = new HashSet<>();

        for (String protein : sortedSpectrumCountsMap.keySet()) {

            if (!remainingProteins.contains(protein)) continue;

            proteinGrp.clear();
            addToProteinGroup(protein,proteinGrp);

            razorProteinMap.put(protein,new HashSet<>(proteinGrp));

            remainingProteins.removeAll(proteinGrp);

            for (String linkedProtein : proteinGrp) {
                proteinGroupMap.put(linkedProtein,protein);
            }
        }
    }


    protected void addToProteinGroup(String protein, Set<String> proteinGrp) {

        if (proteinGrp.contains(protein)) return;

        proteinGrp.add(protein);

        // get linked proteins
        Set<String> peptides = proteinPeptideMap.get(protein);

        for (String peptide : peptides) {
            for (String linkedProt : peptideProteinMap.get(peptide)) {

                addToProteinGroup(linkedProt,proteinGrp);
            }
        }
    }

    public void write(File reportFile) throws IOException {

        if (sortedSpectrumCountsMap.isEmpty()) calcSpectrumCounts();

        FileWriter writer = new FileWriter(reportFile);
        writer.write("protein\tgeneName\tproteinGroup\tgeneNameGroup\tlength\toriginalSpectrumCount\tadjustedSpectrumCount\toriginalPeptideCount\tadjustedPeptideCount\n");

        for (String protein : sortedSpectrumCountsMap.keySet()) {

            String ac = protein;
            if (protein.startsWith("sp|") || protein.startsWith("tr|")) {
                ac = protein.split("\\|")[1];
            }

            String lenStr = "NA";
            String gn = "NA";
            String lgn = "";
            if (uniProtDB!=null && uniProtDB.containsProteinID(ac)) {
                UniProtProtein p = (UniProtProtein) uniProtDB.getProtein(ac);
                lenStr = ""+p.getSequence().length();
                gn = p.getGeneName();

                String leadingProt = proteinGroupMap.get(protein);
                if (leadingProt.startsWith("sp|") || leadingProt.startsWith("tr|")) {
                    ac = leadingProt.split("\\|")[1];

                    p = (UniProtProtein) uniProtDB.getProtein(ac);
                    lenStr = ""+p.getSequence().length();
                    lgn = p.getGeneName();
                }

            }

            writer.write(protein+"\t"+gn+"\t"+proteinGroupMap.get(protein)+"\t"+lgn+"\t"+lenStr+"\t"+
                    sortedSpectrumCountsMap.get(protein)+"\t"+sortedSpectrumCountsMapAdj.get(protein)+"\t"+
                    proteinPeptideMap.get(protein).size()+"\t"+proteinPeptideMapAdj.get(protein).size()+"\n");
            writer.flush();
        }

        writer.close();
    }
}
