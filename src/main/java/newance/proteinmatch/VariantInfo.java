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

/**
 * @author Markus Müller
 */

public class VariantInfo {

    public class VariantSimple {

        final private int position;
        final private char mutatedAA;
        final private String info;

        // \VariantSimple=(141|A|rs75062661_0)
        public VariantSimple(String mutation) {

            mutation = mutation.substring(1,mutation.indexOf(')'));
            String[] fields = mutation.split("\\|");

            this.position = Integer.parseInt(fields[0]);
            this.mutatedAA = fields[1].charAt(0);
            // old versions of fasta files have the format \VariantSimple=(141|A)
            if (fields.length==3)
                this.info = fields[2];
            else
                this.info = "";
        }

        public int getPosition() {
            return position;
        }

        public char getMutatedAA() {
            return mutatedAA;
        }

        public String getInfo() {
            return info;
        }

        public String toString() {
            return String.format("(%d|%s|%s)",position,mutatedAA,info);
        }
    }

    final protected String proteinID;
    final protected List<VariantSimple> variants;


    //>ENSP00000343864.2 \Length=676 \VariantSimple=(579|Q|rs13302979_0) (410|Q|rs34903276_0) (406|E|rs13303368_0) (398|G|rs13302983_0) (663|A|rs7417106_0)
    public VariantInfo(String header) {

        String[] fields = header.split("\\\\");
        this.proteinID = fields[0].substring(1).trim();
        this.variants = new ArrayList<>();
        for (int i=1;i < fields.length;i++) {

            if (fields[i].startsWith("VariantSimple=")) {
                String[] variantStrs = fields[i].split("=")[1].split(" ");
                for (String variant : variantStrs) variants.add(new VariantSimple(variant));
                break;
            }
        }
    }

    public String getProteinID() {
        return proteinID;
    }

    public List<VariantSimple> getVariants() {
        return variants;
    }

    public String toString() {
        String varStr = proteinID;

        if (!variants.isEmpty()) {
            varStr += " \\VariantSimple=";
            for (VariantSimple variant : variants) varStr += variant.toString()+" ";
            return varStr.substring(0,varStr.length()-1);
        } else
            return varStr;

    }

    public boolean hasVariants() {
        return !variants.isEmpty();
    }

    public String getVariantInfo(int start, char mutatedAA, int peptideVarPos) {

        for (VariantSimple variant : variants) {
            // start+peptideVarPos+1 because PEFF format start counting from 1
            if (variant.getPosition() == start+peptideVarPos+1 && variant.getMutatedAA() == mutatedAA) return variant.getInfo();
        }

        return "";
    }
}
