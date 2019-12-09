package newance.proteinmatch;

import java.io.*;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.logging.Logger;

/**
 * Created by markusmueller on 17.11.17.
 */
public class UniProtDB implements Serializable {
    private static final Logger LOGGER = Logger.getLogger(UniProtDB.class.getName());

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

    protected final List<String> fastaFileNames;
    protected final List<UniProtProtein> proteins;
    protected final Map<String,UniProtProtein> acIndex;
    protected final Map<Integer,Set<TagMatch>> sequenceIndex;
    protected final int seqTagLength;
    protected final ConcurrentHashMap<String, Map<String, List<PeptideUniProtSequenceMatch>>> searchedPeptides;


    public UniProtDB(String[] fastaLines) {

        this.fastaFileNames = new ArrayList<>();
        proteins = new ArrayList<>();
        acIndex = new HashMap<>();
        sequenceIndex = new HashMap<>();
        searchedPeptides = new ConcurrentHashMap<>();
        seqTagLength = 6;

        readProteins(fastaLines);
        indexProteins();
    }

    public UniProtDB(String fastaFileName) {

        this.fastaFileNames = new ArrayList<>();
        this.fastaFileNames.add(fastaFileName);
        proteins = new ArrayList<>();
        acIndex = new HashMap<>();
        sequenceIndex = new HashMap<>();
        searchedPeptides = new ConcurrentHashMap<>();
        seqTagLength = 6;

        readProteins(fastaFileName);
        indexProteins();
    }

    public UniProtDB(List<String> fastaFileNames) {

        this.fastaFileNames = new ArrayList<>(fastaFileNames);
        proteins = new ArrayList<>();
        acIndex = new HashMap<>();
        sequenceIndex = new HashMap<>();
        searchedPeptides = new ConcurrentHashMap<>();
        seqTagLength = 6;

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

    public List<UniProtProtein> getProteins() {
        return proteins;
    }

    public boolean containsUniProtAC(String uniProtAC) {
        return acIndex.containsKey(uniProtAC);
    }

    public List<String> getFastaFileNames() {
        return fastaFileNames;
    }

    private void readProteins(String[] fastaLines) {

        Map<String,String> headerFieldMap = null;
        String seq = "";
        for (String line : fastaLines) {

            if (line.isEmpty()) continue;

            if (line.startsWith(">")) {

                addFastaEntry(headerFieldMap,seq);
                headerFieldMap = parseHeader(line);
                seq = "";
            } else {
                seq += line;
            }
        }

        addFastaEntry(headerFieldMap,seq);
    }

    private void readProteins(String fastaFileName) {

        try {
            BufferedReader lineReader = new BufferedReader(new FileReader(new File(fastaFileName)));

            String line;
            Map<String,String> headerFieldMap = null;
            String seq = "";
            while ((line = lineReader.readLine()) != null) {

                if (line.isEmpty()) continue;

                if (line.startsWith(">")) {

                    addFastaEntry(headerFieldMap,seq);
                    headerFieldMap = parseHeader(line);
                    seq = "";
                } else {
                    seq += line;
                }
            }

            addFastaEntry(headerFieldMap,seq);

            lineReader.close();

        } catch (FileNotFoundException e) {
            LOGGER.severe(e.getMessage());
        }  catch (IOException e) {
            LOGGER.severe(e.getMessage());
        }

    }


    protected void addFastaEntry(Map<String,String>  headerFieldMap, String seq) {

        if (headerFieldMap != null) {
            UniProtProtein protein = new UniProtProtein(
                    headerFieldMap.get("uniProtAC"),
                    headerFieldMap.get("standardName"),
                    headerFieldMap.get("description"),
                    headerFieldMap.get("geneName"),
                    headerFieldMap.get("dbFlag"),
                    seq);
            proteins.add(protein);
            add2Index(protein);
        }
    }

    protected static Map<String,String> parseHeader(String header) {

        if (header.startsWith(">sp|") || header.startsWith(">tr|")) {
            return parseUniProtHeader(header);
        } else {
            return parseOtherHeader(header);
        }
    }

    //>ENSP00000334393.3 \VariantSimple=(141|A)
    //>LINE_L1_L1MC4_23_154353601_154353761.1-117
    protected static Map<String,String> parseOtherHeader(String header) {

        String[] fields = header.split("\\s");
        String id = fields[0].substring(1);
        Map<String,String> headerFieldMap = new HashMap<>();

        headerFieldMap.put("dbFlag","other");
        headerFieldMap.put("uniProtAC",id);

        headerFieldMap.put("standardName",id);
        if (fields.length>1) {

            headerFieldMap.put("description",fields[1]);
        } else {

            headerFieldMap.put("description","");
        }

        headerFieldMap.put("geneName","NA");

        return headerFieldMap;
    }


    //>tr|A0A024R3B9|A0A024R3B9_HUMAN Alpha-crystallin B chain OS=Homo sapiens GN=CRYAB PE=3 SV=1
    //>sp|A0A183|LCE6A_HUMAN Late cornified envelope protein 6A OS=Homo sapiens GN=LCE6A PE=2 SV=1
    protected static Map<String,String> parseUniProtHeader(String header) {

        String[] fields = header.split("\\|");

        Map<String,String> headerFieldMap = new HashMap<>();

        headerFieldMap.put("dbFlag",fields[0].substring(1,3));
        headerFieldMap.put("uniProtAC",fields[1]);

        int idx = fields[2].indexOf(" ");
        headerFieldMap.put("standardName",fields[2].substring(0,idx));
        headerFieldMap.put("description",fields[2].substring(idx+1));

        int geneNameStart = fields[2].indexOf("GN=")+3;
        int geneNameEnd = fields[2].indexOf(" ",geneNameStart);
        if (geneNameEnd<0) geneNameEnd = fields[2].length();
        headerFieldMap.put("geneName",fields[2].substring(geneNameStart,geneNameEnd));

        return headerFieldMap;
    }


    protected void indexProteins() {

        for (UniProtProtein p : proteins) {
            acIndex.put(p.getUniProtAC(),p);
        }
    }

    public UniProtProtein getProtein(String uniProtAC) {

        if (acIndex.containsKey(uniProtAC))
            return acIndex.get(uniProtAC);
        else
            return null;
    }

    private List<PeptideUniProtSequenceMatch> findPeptideTagMatch(TagMatch tagMatch, char[] peptideSeq){

        UniProtProtein protein = getProtein(tagMatch.getAc());
        return findPeptide(protein,tagMatch.getPos(), peptideSeq);
    }

    private List<PeptideUniProtSequenceMatch> findPeptide(UniProtProtein protein, int pos, char[] peptideSeq){

        List<PeptideUniProtSequenceMatch> ppms = new ArrayList<>();

        int matchPos = indexOf(protein.getCharArray(),peptideSeq,pos);
        if (matchPos >= 0) {
            ppms.add(new PeptideUniProtSequenceMatch(matchPos,matchPos+peptideSeq.length,peptideSeq,protein));
        }

        return ppms;
    }

    protected boolean contains(TagMatch tagMatch, char[] peptideSeq){

        UniProtProtein protein = getProtein(tagMatch.getAc());
        if (protein==null) return false;

        int pos = indexOf(protein.getCharArray(),peptideSeq,tagMatch.getPos());

        return pos>=0;
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
        String peptideSeqIL = peptideSeq.replace("I","L");

        if (searchedPeptides.containsKey(peptideSeqIL))
            return searchedPeptides.get(peptideSeqIL);

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

        searchedPeptides.put(peptideSeqIL,ppms);

        return ppms;

    }

    public void clear() {
        this.fastaFileNames.clear();
        this.proteins.clear();
        this.acIndex.clear();
        this.sequenceIndex.clear();
        this.searchedPeptides.clear();
    }

    private void add2Index(UniProtProtein protein) {
        char[] sequence = protein.getCharArray();
        for (int i=0;i<sequence.length-seqTagLength;i++) {

            int tag =  getHash(sequence,i,i+seqTagLength-1);

            if (!sequenceIndex.containsKey(tag)) {

                sequenceIndex.put(tag,new HashSet<>());
            }
            sequenceIndex.get(tag).add(new TagMatch(protein.getUniProtAC(),i));
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

}
