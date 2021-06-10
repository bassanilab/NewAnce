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

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

/**
 * @author Markus Müller
 */

public class VariantInfo {

    final protected String proteinID;
    final protected List<SequenceVariant> variants;

    private VariantInfo(String proteinID, List<SequenceVariant> variants) {

        this.proteinID = proteinID;
        this.variants = variants;
    }

    public String getProteinID() {
        return proteinID;
    }

    public List<SequenceVariant> getVariants() {
        return variants;
    }

    public int size() {return variants.size();}

    public SequenceVariant get(int idx) {
        return variants.get(idx);
    }

    public String toString() {
        String varStr = proteinID + " ";

        if (!variants.isEmpty()) {
            for (SequenceVariant variant : variants) varStr += variant.toString()+" ";
            return varStr.substring(0,varStr.length()-1);
        } else
            return varStr;

    }

    public boolean hasVariants() {
        return !variants.isEmpty();
    }

    //>ENSP00000343864.2 \Length=676 \VariantSimple=(579|Q|rs13302979_0) (410|Q|rs34903276_0) (406|E|rs13303368_0) (398|G|rs13302983_0) (663|A|rs7417106_0)
    public static VariantInfo parseFastaHeader(String header) {

        int pos1 = header.indexOf("\\VariantSimple=");
        int pos2 = header.indexOf("\\VariantComplex=");
        int pos4 = header.indexOf("\\Length=");

        int pos3 = header.indexOf(" ");
        pos3 = (pos3 > 0)?pos3:header.length();

        String proteinID = header.substring(1,pos3);

        int length = -1;
        if (pos4 >= 0) {
            pos3 = header.indexOf('\\', pos4 + 1);
            pos3 = (pos3 < 0) ? header.length() : pos3;
            length = Integer.parseInt(header.substring(pos4 + "\\Length=".length(), pos3).trim());
        }

        List<SequenceVariant> variants = new ArrayList<>();

        if (pos1>=0) {
            pos3 = header.indexOf('\\', pos1 + 1);
            pos3 = (pos3 < 0) ? header.length() : pos3;
            String varSimpleStr = header.substring(pos1 + "\\VariantSimple=".length(), pos3);

            parseSimpleVariants(proteinID, length, varSimpleStr, variants);
        }

        if (pos2>=0) {
            pos3 = header.indexOf('\\', pos2 + 1);
            pos3 = (pos3 < 0) ? header.length() : pos3;
            String varComplexStr = header.substring(pos2 + "\\VariantComplex=".length(), pos3);

            parseComplexVariants(proteinID, length, varComplexStr, variants);
        }

        return new VariantInfo(proteinID, variants);
    }

    private static void parseSimpleVariants(String proteinID, int proteinLength,
                                            String varSimpleStr,List<SequenceVariant> variants) {

        int pos = 0;
        int pos2 = varSimpleStr.indexOf(')')+1;

        while (pos>-1) {
            if (pos2-pos > 5)
                variants.add(SequenceVariant.parseSimpleVariantString(proteinID, proteinLength,
                        varSimpleStr.substring(pos,pos2)));
            pos = varSimpleStr.indexOf('(',pos+1);;
            pos2 = varSimpleStr.indexOf(')',pos+1)+1;
        }
    }

    private static void parseComplexVariants(String proteinID, int proteinLength,
                                             String varComplexStr,List<SequenceVariant> variants) {

        int pos = 0;
        int pos2 = varComplexStr.indexOf(')')+1;

        while (pos>-1) {
            if (pos2-pos > 5)
                variants.add(SequenceVariant.parseComplexVariantString(proteinID, proteinLength,
                        varComplexStr.substring(pos,pos2)));
            pos = varComplexStr.indexOf('(',pos+1);;
            pos2 = varComplexStr.indexOf(')',pos+1)+1;
        }
    }
}
