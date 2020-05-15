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


/**
 * @author Markus Müller
 */

public class VariantProtein extends FastaProtein {
    final protected VariantInfo variantInfo;

    public VariantProtein(VariantInfo variantInfo, String sequence) {
        super(variantInfo.getProteinID(), sequence);
        this.variantInfo = variantInfo;
    }

    public VariantProtein(VariantInfo variantInfo) {
        super(variantInfo.getProteinID());
        this.variantInfo = variantInfo;
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
}
