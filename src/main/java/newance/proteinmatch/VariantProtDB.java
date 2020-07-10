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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Markus Müller
 */

public class VariantProtDB extends FastaDB {

    private Map<String, List<VariantProtein>> variantProteinMap;
    private Map<String, List<SequenceVariant>> variantInfoMap;

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

        if (variantProteinMap == null) variantProteinMap = new HashMap<>();
        if (variantInfoMap == null) variantInfoMap = new HashMap<>();

        VariantProtein variantProtein = (VariantProtein)protein;
        protein.setSequence(seq);
        proteins.add(protein);
//        add2Index(protein);

        for (SequenceVariant sequenceVariant : variantProtein.getVariantInfo().getVariants()) {
            variantProteinMap.putIfAbsent(sequenceVariant.getProteinID(), new ArrayList<>());
            variantProteinMap.get(sequenceVariant.getProteinID()).add(variantProtein);
            variantInfoMap.putIfAbsent(sequenceVariant.getInfo(), new ArrayList<>());
            variantInfoMap.get(sequenceVariant.getInfo()).add(sequenceVariant);
        }
    }


    @Override
    protected FastaProtein parseHeader(String header) {

        VariantInfo variantInfo = VariantInfo.parseFastaHeader(header);

        return new VariantProtein(variantInfo);
    }

    public List<SequenceVariant> getSequenceVariant(String variantID) {
        return variantInfoMap.get(variantID);
    }

    public boolean containsVariant(String variantID) {
        return variantInfoMap.containsKey(variantID);
    }
}
