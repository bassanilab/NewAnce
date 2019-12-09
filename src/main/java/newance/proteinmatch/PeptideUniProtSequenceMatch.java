package newance.proteinmatch;

/**
 * Created by markusmueller on 16.11.17.
 */
public class PeptideUniProtSequenceMatch {
    protected final int start;
    protected final int end;
    protected final UniProtProtein protein;
    protected final char[] peptide;


    public PeptideUniProtSequenceMatch(int start, int end, char[] peptide, UniProtProtein protein) {
        this.start = start;
        this.end = end;

        this.protein = protein;
        this.peptide = peptide;
    }

    public String toString() {
        return peptide+"\t"+protein.getUniProtAC()+"\t"+start+"\t"+end;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public UniProtProtein getProtein() {
        return protein;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        PeptideUniProtSequenceMatch that = (PeptideUniProtSequenceMatch) o;

        if (start != that.start) return false;
        if (end != that.end) return false;
        if (!protein.equals(that.protein)) return false;
        return peptide.equals(that.peptide);
    }

    @Override
    public int hashCode() {
        int result = start;
        result = 31 * result + end;
        result = 31 * result + protein.hashCode();
        result = 31 * result + peptide.hashCode();
        return result;
    }
}

