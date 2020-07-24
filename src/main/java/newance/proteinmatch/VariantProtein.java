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


import newance.psmconverter.PeptideSpectrumMatch;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Markus Müller
 */

public class VariantProtein extends FastaProtein {
    protected final VariantInfo variantInfo;
    protected int[][] variantPos;
    protected final List<SequenceVariant> variants;
    protected int peptideStart;
    protected int peptideEnd;


    public VariantProtein(VariantInfo variantInfo, String sequence) {
        super(variantInfo.getProteinID(), sequence);
        this.variantInfo = variantInfo;
        this.variantPos = new int[this.sequence.length][];
        this.variants = new ArrayList<>();
        this.peptideStart = -1;

        setVariantPos();
    }

    public VariantProtein(VariantInfo variantInfo) {
        super(variantInfo.getProteinID());
        this.variantInfo = variantInfo;
        this.variants = new ArrayList<>();
        this.peptideStart = -1;
    }

    public String toString() {
        return variantInfo.toString();
    }

    public VariantInfo getVariantInfo() {
        return variantInfo;
    }

    public boolean hasVariants() {
        return variantInfo.hasVariants();
    }

    public int getPeptideStart() {
        return peptideStart;
    }

    public int getPeptideEnd() {
        return peptideEnd;
    }

    public void setSequence(String sequence) {
        super.setSequence(sequence);
        this.variantPos = new int[this.sequence.length+1][]; // insertion can be at end of sequence

        setVariantPos();
    }

    private void setVariantPos() {

        for (int i=0; i<variantInfo.size(); i++) {
            SequenceVariant variant = variantInfo.get(i);
            int j=variant.getStartWT()-1; // variantInfo starts position counting from 1

            if (j >= variantPos.length) continue;

            if (variantPos[j] == null) {
                variantPos[j] = new int[5]; // max 5 variant per pos
                variantPos[j][0] = i;
                for (int k=1;k<5;k++) variantPos[j][k] = -1; // -1:  no variant; >=0 variant index
            } else {
                int k;
                for (k=0;k<5;k++) if (variantPos[j][k] == -1) break;
                if (k<5) variantPos[j][k] = i;  //
            }
        }

    }

    public List<SequenceVariant>  matchWithVariant(String variantSeq) {

        if (variantSeq.isEmpty()) return null;

        variants.clear();

        String s = new String(sequence);
        if (s.contains(variantSeq)) return variants;

        char startCh = variantSeq.charAt(0);
        char[] variantChars = variantSeq.toCharArray();

        for (int i=0;i<variantPos.length;i++) {
            if (variantPos[i] == null && i < sequence.length && sequence[i] != startCh) continue;

            if (matchWithVariant(variantChars, variantSeq, i, 0)) {
                peptideStart = i;
                return variants;
            }

        }

        return variants;
    }

    private boolean matchWithVariant(char[] variantChars, String variantSeq, int proteinPos, int peptidePos) {

        if (proteinPos >= sequence.length) return false;

        SequenceVariant variant = null;

        int vi = peptidePos;
        int si = proteinPos;

        while (vi < variantChars.length && si < sequence.length) {
            if (variantPos[si] != null) {
                for (int k=0;k<5;k++) {
                    if (variantPos[si][k] == -1) break;

                    variant = variantInfo.get(variantPos[si][k]);

                    if (variant.match(variantSeq, vi)) {

                        if (matchWithVariant(variantChars, variantSeq, variant.getEndWT(),
                                variant.getPosAfterVariant())) {
                            variants.add(variant);
                            return true; // only look for first variant that matches
                        }
                    }
                }
            }

            // if no variant match
            if (sequence[si] != variantChars[vi]) return false;

            vi++;
            si++;
        }

        peptideEnd = si;
        return vi == variantChars.length;
    }

    public String getWTSequence() {

        return String.copyValueOf(sequence, peptideStart, peptideEnd-peptideStart);
    }


}
