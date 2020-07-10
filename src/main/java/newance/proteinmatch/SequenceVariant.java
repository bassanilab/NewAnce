package newance.proteinmatch;

/**
 * Created by markusmueller on 16.06.20.
 */
public class SequenceVariant {

    public enum VariantType {AA_SUBSTITUTION, SUBSTITUTION, INSERTION, DELETION, STOP, UNKNOWN}

    final private String proteinID;
    final private int proteinLength;
    final private int startWT;
    final private int endWT;
    private String mutatedSequence;
    final private int length;
    final private String info;
    private VariantType type;

    private int variantPos;

    // for simple variants
    public SequenceVariant(String proteinID, int proteinLength, int position, String mutatedSequence, String info) {

        if (mutatedSequence.length()!=1)
            throw new IllegalArgumentException("SequenceVariant: mutatedSequence must have length 1");

        this.proteinID = proteinID;
        this.proteinLength = proteinLength;
        this.startWT = position;
        this.endWT = position;
        this.info = info;

        setSimpleVariantParams(proteinLength, position, mutatedSequence);

        this.length = this.mutatedSequence.length();
    }

    // for complex variants
    public SequenceVariant(String proteinID, int proteinLength, int startWT, int endWT, String mutatedSequence,
                           String info) {

        if (endWT < startWT)
            throw new IllegalArgumentException("SequenceVariant: startWT must not be larger than endWT");

        this.proteinID = proteinID;
        this.proteinLength = proteinLength;
        this.startWT = startWT;
        this.endWT = endWT;
        this.info = info;

        int len = mutatedSequence.length();

        if (startWT==endWT && len==1) {
            setSimpleVariantParams(proteinLength, startWT, mutatedSequence);
        } else if (mutatedSequence.endsWith("*")) {
            this.mutatedSequence = mutatedSequence.substring(0,len-1);
            this.type = VariantType.STOP;
        } else {

            this.mutatedSequence = mutatedSequence;

            if (len == endWT - startWT + 1)
                this.type = VariantType.SUBSTITUTION;
            else if (len > endWT - startWT + 1)
                this.type = VariantType.INSERTION;
            else if (len < endWT - startWT + 1)
                this.type = VariantType.DELETION;
        }

        this.length = this.mutatedSequence.length();
    }

    private void setSimpleVariantParams(int proteinLength, int position, String mutatedSequence) {
        if (mutatedSequence.equals("*")) {
            this.mutatedSequence = "";
            this.type = VariantType.STOP;
        }
        else if (position <= proteinLength) {
            this.mutatedSequence = mutatedSequence;
            this.type = VariantType.AA_SUBSTITUTION;
        }
        else {
            // insertion at end of protein (override stop)
            this.mutatedSequence = mutatedSequence;
            this.type = VariantType.INSERTION;
        }
    }

    public int getStartWT() {
        return startWT;
    }

    public int getEndWT() {
        return endWT;
    }

    public int getLength() {
        return length;
    }

    public VariantType getType() {
        return type;
    }

    public String getProteinID() {
        return proteinID;
    }

    public String getMutatedSequence() {
        return mutatedSequence;
    }

    public String getInfo() {
        return info;
    }

    public int getProteinLength() {
        return proteinLength;
    }

    public int getPosAfterVariant() {
        return variantPos;
    }

    public String toString() {
        return type +":" + String.format("%d,%d,%s,%s",startWT, endWT, mutatedSequence,info);
    }

    // (141|A|rs75062661_0)
    // (24|*|rs6671527_0)  : stop introduced at position 24
    // (205|R|rs4970490_0) : protein length is 204 => R is inserted at end
    public static SequenceVariant parseSimpleVariantString(String proteinID, int proteinLength, String peffString) {

        String mutationStr = peffString.trim().substring(1,peffString.indexOf(')'));

        int pos1 = mutationStr.indexOf('|');
        int pos2 = mutationStr.indexOf('|',pos1+1);

        int position = Integer.parseInt(mutationStr.substring(0,pos1));
        String mutatedSequence = mutationStr.substring(pos1+1,pos2);

        String info = "";
        //  \VariantSimple=(141|A) : info can be missing
        if (pos2>=0 && mutationStr.length()-pos2>0) info = mutationStr.substring(pos2+1);

        return new SequenceVariant(proteinID, proteinLength, position, mutatedSequence, info);
    }

    // (117|119|GL|rs10578519_3)
    public static SequenceVariant parseComplexVariantString(String proteinID, int proteinLength, String peffString) {

        String mutationStr = peffString.trim().substring(1,peffString.indexOf(')'));

        int pos1 = mutationStr.indexOf('|');
        int pos2 = mutationStr.indexOf('|',pos1+1);
        int pos3 = mutationStr.indexOf('|',pos2+1);

        int start = Integer.parseInt(mutationStr.substring(0,pos1));
        int end = Integer.parseInt(mutationStr.substring(pos1+1,pos2));
        String mutatedSequence = mutationStr.substring(pos2+1,pos3);

        String info = "";
        //  \VariantComplex=(117|119|GL) : info can be missing
        if (pos3>=0 && mutationStr.length()-pos3>0) info = mutationStr.substring(pos3+1);

        return new SequenceVariant(proteinID, proteinLength, start, end, mutatedSequence, info);
    }

    public boolean match(String variantStr, int start) {

        variantPos = -1;

        if (length==0) {
            variantPos = start;
            return true;
        }

        int len2 = variantStr.length();

        if (type==VariantType.STOP) { // variantSeq must stop
            if (len2 - start > length) return false; // variantSeq too long to match
            else if (start > 0) {
                boolean match = mutatedSequence.startsWith(variantStr.substring(start));
                if (match) variantPos = len2;
                return match; // variantSeq has to be at start
            }
            else {
                boolean match = mutatedSequence.contains(variantStr);
                if (match) variantPos = len2;
                return match; // variantSeq can be anywhere
            }
        } else { // variantSeq can go on
            if (start == 0) return containsWithoutStop(variantStr); // variantSeq can start anywhere
            else {
                int len1 = mutatedSequence.length();
                variantPos = (len2-start<len1)?len2:start+len1;
                if (mutatedSequence.startsWith(variantStr.substring(start, variantPos))) {
                    return true;
                } else {
                    variantPos = -1;
                    return false;
                }
            }
        }
    }

    private boolean containsWithoutStop(String variantStr) {

        int len1 = mutatedSequence.length();
        int len2 = variantStr.length();
        for (int i=0; i < mutatedSequence.length(); i++) {
            int end1 = (i+len2>len1)?len1:i+len2;
            variantPos = (len1-i>len2)?len2:len1-i;
            if (mutatedSequence.substring(i,end1).equals(variantStr.substring(0,variantPos))) {
                return true;
            } else {
                variantPos = -1;
            }
        }

        return false;
    }

/*
    // rr5840_0:Y17C -> single aa change
    // rs78783575_0:G163GLAGRPRRGGRGARARP*  -> frame shift
    // rs4970490_0:*205R  -> insertion at end
    // rs144636354_0:K4KRK -> insertion
    // rs6671527_0:R24* -> deletion of R at 24
    // rs10578519_2:GLR75GL -> deletion
    // rs72362780_0:TSQQPSPESTP27TS -> deletion
    // rr189_0:W50CVWWCGSLEAVRESKGGNYQLCTRKKI
    public static SequenceVariant parseMaxQuantString(String mqString) {

        int pos1 = mqString.indexOf(':');
        String info = mqString.substring(0,pos1);
        int pos2 = mutationStr.indexOf('|',pos1+1);

        int position = Integer.parseInt(mutationStr.substring(0,pos1));
        String mutatedSequence = mutationStr.substring(pos1+1,pos1+2);

        String info = "";
        // old versions of fasta files have the format \VariantSimple=(141|A)
        if (pos2>=0 && mutationStr.length()-pos2>0) info = mutationStr.substring(pos2+1);

        return new SequenceVariant(position, mutatedSequence, info);
    }
*/
}
