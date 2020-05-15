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

import java.io.*;
import java.util.*;

/**
 * @author Markus Müller
 */

public abstract class FastaDB {

    protected class TagMatch {
        private final String ac;
        private final int pos;

        public TagMatch(String ac, int pos) {
            this.ac = ac;
            this.pos = pos;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof TagMatch)) return false;

            TagMatch tagMatch = (TagMatch) o;

            if (pos != tagMatch.pos) return false;
            return ac.equals(tagMatch.ac);
        }

        @Override
        public int hashCode() {
            int result = ac.hashCode();
            result = 31 * result + pos;
            return result;
        }

        public String getAc() {
            return ac;
        }

        public int getPos() {
            return pos;
        }
    }

    protected List<String> fastaFileNames;
    protected List<FastaProtein> proteins;
    protected Map<String,FastaProtein> acIndex;
    protected Map<Integer,Set<TagMatch>> sequenceIndex;
    protected int seqTagLength;

    public FastaDB(int seqTagLength) {

        this.fastaFileNames = new ArrayList<>();
        this.proteins = new ArrayList<>();
        this.acIndex = new HashMap<>();
        this.sequenceIndex = new HashMap<>();
        this.seqTagLength = seqTagLength;
    }

    public FastaDB(String[] fastaLines, int seqTagLength) {

        this.fastaFileNames = new ArrayList<>();
        this.proteins = new ArrayList<>();
        this.acIndex = new HashMap<>();
        this.sequenceIndex = new HashMap<>();
        this.seqTagLength = seqTagLength;

        readProteins(fastaLines);
        indexProteins();
    }

    public FastaDB(String fastaFileName, int seqTagLength) {

        this.fastaFileNames = new ArrayList<>();
        this.fastaFileNames.add(fastaFileName);
        this.proteins = new ArrayList<>();
        this.acIndex = new HashMap<>();
        this.sequenceIndex = new HashMap<>();
        this.seqTagLength = seqTagLength;

        readProteins(fastaFileName);
        indexProteins();
    }

    public FastaDB(List<String> fastaFileNames, int seqTagLength) {

        this.fastaFileNames = new ArrayList<>(fastaFileNames);
        this.proteins = new ArrayList<>();
        this.acIndex = new HashMap<>();
        this.sequenceIndex = new HashMap<>();
        this.seqTagLength = seqTagLength;

        for (String fastaFileName : this.fastaFileNames) {
            readProteins(fastaFileName);
        }

        indexProteins();
    }

    public void addFastaFiles(List<String> fastaFileNames) {

        for (String fastaFileName : this.fastaFileNames) {
            readProteins(fastaFileName);
        }

        indexProteins();
    }

    public void addFastaFile(String fastaFileName) {

        readProteins(fastaFileName);
        indexProteins();
    }

    public List<String> getFastaFileNames() {
        return fastaFileNames;
    }

    protected void readProteins(String[] fastaLines) {

        FastaProtein protein = null;
        String seq = "";
        for (String line : fastaLines) {

            if (line.isEmpty()) continue;

            if (line.startsWith(">")) {

                addFastaEntry(protein,seq);
                protein = parseHeader(line);
                seq = "";
            } else {
                seq += line;
            }
        }

        addFastaEntry(protein,seq);
    }

    protected void readProteins(String fastaFileName) {

        try {
            BufferedReader lineReader = new BufferedReader(new FileReader(new File(fastaFileName)));

            String line;
            FastaProtein protein = null;
            String seq = "";
            while ((line = lineReader.readLine()) != null) {

                if (line.isEmpty()) continue;

                if (line.startsWith(">")) {

                    addFastaEntry(protein,seq);
                    protein = parseHeader(line);
                    seq = "";
                } else {
                    seq += line;
                }
            }

            addFastaEntry(protein,seq);

            lineReader.close();

        } catch (FileNotFoundException e) {
            System.err.println(e.getMessage());
        }  catch (IOException e) {
            System.err.println(e.getMessage());
        }

    }

    protected void addFastaEntry(FastaProtein protein, String seq) {

        if (protein==null || seq.isEmpty()) return;

        protein.setSequence(seq);
        proteins.add(protein);
        add2Index(protein);
    }

    protected abstract FastaProtein parseHeader(String header);

    protected void indexProteins() {

        for (FastaProtein p : proteins) {
            acIndex.put(p.getProteinID(),p);
        }
    }


    protected List<PeptideUniProtSequenceMatch> findPeptideTagMatch(TagMatch tagMatch, char[] peptideSeq){

        FastaProtein protein = getProtein(tagMatch.getAc());
        return findPeptide(protein,tagMatch.getPos(), peptideSeq);
    }

    protected List<PeptideUniProtSequenceMatch> findPeptide(FastaProtein protein, int pos, char[] peptideSeq){

        List<PeptideUniProtSequenceMatch> ppms = new ArrayList<>();

        int matchPos = indexOf(protein.getCharArray(),peptideSeq,pos);
        if (matchPos >= 0) {
            ppms.add(new PeptideUniProtSequenceMatch(matchPos,matchPos+peptideSeq.length,peptideSeq,protein));
        }

        return ppms;
    }

    protected boolean contains(TagMatch tagMatch, char[] peptideSeq){

        FastaProtein protein = getProtein(tagMatch.getAc());
        if (protein==null) return false;

        int pos = indexOf(protein.getCharArray(),peptideSeq,tagMatch.getPos());

        return pos>=0;
    }

    public boolean containsProteinID(String proteinID) {
        return acIndex.containsKey(proteinID);
    }

    public FastaProtein getProtein(String proteinID) {

        if (acIndex.containsKey(proteinID))
            return acIndex.get(proteinID);
        else
            return null;
    }

    public boolean contains(char[] peptideSeq){

        if (peptideSeq.length<4) return false;

        Set<TagMatch> tagMatches;

        int idx = getHash(peptideSeq,0,seqTagLength-1);
        if (sequenceIndex.containsKey(idx)) {
            tagMatches = new HashSet<>(sequenceIndex.get(idx));
        } else {
            return false;
        }

        for (TagMatch tagMatch : tagMatches) {
            if (contains(tagMatch,peptideSeq))
                return true;
        }

        return false;
    }

    public Map<String,List<PeptideUniProtSequenceMatch>> findPeptide(String peptideSeq){

        if (peptideSeq.length()<seqTagLength) return new HashMap<>();

        char[] aas = peptideSeq.toCharArray();

        Set<TagMatch> tagMatches;
        int hash = getHash(aas,0,seqTagLength-1);

        if (sequenceIndex.containsKey(hash)) {
            tagMatches = new HashSet<>(sequenceIndex.get(hash));
        } else {
            return new HashMap<>();
        }

        Map<String,List<PeptideUniProtSequenceMatch>> ppms = new HashMap<>();
        for (TagMatch tagMatch : tagMatches) {

            if (contains(tagMatch,aas)) {
                if (ppms.containsKey(tagMatch.getAc())) {
                    ppms.get(tagMatch.getAc()).addAll(findPeptideTagMatch(tagMatch,aas));
                } else
                    ppms.put(tagMatch.getAc(), findPeptideTagMatch(tagMatch,aas));
            }
        }

        return ppms;

    }

    public List<FastaProtein> getProteins() {
        return proteins;
    }

    public int getSeqTagLength() {
        return seqTagLength;
    }

    protected void add2Index(FastaProtein protein) {
        char[] charArray = protein.getCharArray();
        for (int i=0;i<charArray.length-seqTagLength;i++) {

            int tag =  getHash(charArray,i,i+seqTagLength-1);

            if (!sequenceIndex.containsKey(tag)) {

                sequenceIndex.put(tag,new HashSet<>());
            }
            sequenceIndex.get(tag).add(new TagMatch(protein.getProteinID(),i));
        }
    }

    protected static int getHash(char[] seq, int start, int end) {

        end = (end>seq.length-1)?seq.length-1:end;

        int hash = 0;
        for (int i=start;i<=end;i++) {
            int v = (seq[i]=='I')?(int)'L':(int)seq[i];
            v -= (int)'A';
            hash = hash*(int)'Z'+v;
        }

        return hash;
    }

    protected int indexOf(char[] source, char[] target, int fromIndex) {

        if (fromIndex+target.length>source.length) return -1;

        for (int i = seqTagLength; i < target.length; i++) {
            if (!match(source[fromIndex+i],target[i])) return -1;
        }

        return fromIndex;
    }

    protected static boolean match(char c1, char c2) {
        if (c1==c2) return true;

        if (c1=='I' && c2=='L') return true;
        if (c1=='L' && c2=='I') return true;

        return false;
    }

    public int size() {
        return proteins.size();
    }

}
