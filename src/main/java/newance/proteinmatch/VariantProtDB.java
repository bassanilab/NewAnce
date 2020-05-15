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

import java.util.List;

/**
 * @author Markus Müller
 */

public class VariantProtDB extends FastaDB {

    public VariantProtDB() {

        super(6);
    }

    public VariantProtDB(String[] fastaLines) {

        super( fastaLines, 6);
    }

    public VariantProtDB(String fastaFileName) {

        super(fastaFileName, 6);
    }

    public VariantProtDB(List<String> fastaFileNames) {

        super(fastaFileNames, 6);
    }

    protected void addFastaEntry(FastaProtein protein, String seq) {

        if (protein==null || seq.isEmpty()) return;

        VariantProtein variantProtein = (VariantProtein)protein;
        if (variantProtein.hasVariants()) {
            protein.setSequence(seq);
            proteins.add(protein);
            add2Index(protein);
        }
    }


    @Override
    protected FastaProtein parseHeader(String header) {

        VariantInfo variantInfo = new VariantInfo(header);

        return new VariantProtein(variantInfo);
    }
}
