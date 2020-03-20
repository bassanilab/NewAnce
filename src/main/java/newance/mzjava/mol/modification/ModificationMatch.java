/**
 * Copyright (c) 2010, SIB. All rights reserved.
 *
 * SIB (Swiss Institute of Bioinformatics) - http://www.isb-sib.ch Host -
 * http://mzjava.expasy.org
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer. Redistributions in binary
 * form must reproduce the above copyright notice, this list of conditions and
 * the following disclaimer in the documentation and/or other materials provided
 * with the distribution. Neither the name of the SIB/GENEBIO nor the names of
 * its contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL SIB/GENEBIO BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package newance.mzjava.mol.modification;


import newance.mzjava.mol.AminoAcid;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

/**
 * Class to hold information on modifications that were found for a PSM. Each modification
 * match has the measured mass shift and a list of 0 or more candidate Modifications.
 *
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class ModificationMatch {

    private final double massShift;
    private final List<Modification> modificationCandidates = new ArrayList<>();
    private final AminoAcid residue;
    private final int position;
    private final ModAttachment modAttachment;

    /**
     * Create a ModificationMatch given a mass shift
     *
     * @param massShift the modifications mass (negative if loss)
     */
    public ModificationMatch(double massShift, AminoAcid residue, Integer position, ModAttachment modAttachment) {

        this.massShift = massShift;
        this.residue = residue;
        this.position = position;
        this.modAttachment = modAttachment;
    }

    /**
     * Create a ModificationMatch given a modification. The massShift is set from
     * the modifications mass.
     *
     * @param modification the modification
     */
    public ModificationMatch(Modification modification, AminoAcid residue, int position, ModAttachment modAttachment) {

        this.massShift = modification.getMolecularMass();
        this.modificationCandidates.add(modification);
        this.residue = residue;
        this.position = position;
        this.modAttachment = modAttachment;
    }

    /**
     * Add a modification to the list of potential modifications
     *
     * @param modification the modification to add
     */
    public void addPotentialModification(Modification modification) {

        modificationCandidates.add(modification);
    }

    /**
     * Returns the mass shift
     *
     * @return the mass shift
     */
    public double getMassShift() {

        return massShift;
    }


    /**
     * Checks whether all modification candidates are within tolerance of the massShift.
     *
     * @param tolerance the tolerance
     * @return true if all modification candidates are within tolerance, false otherwise
     */
    public boolean isWithinTolerance(Tolerance tolerance) {

        for(Modification mod : modificationCandidates) {

            if(!tolerance.withinTolerance(massShift, mod.getMolecularMass()))
                return false;
        }
        return true;
    }

    public AminoAcid getResidue() {
        return residue;
    }

    public int getPosition() {
        return position;
    }

    public ModAttachment getModAttachment() {
        return modAttachment;
    }

    @Override
    public String toString() {

        if (modificationCandidates.isEmpty()) {

            return residue + " : " + Double.toString(massShift);
        } else {

            return residue + " : " + Double.toString(massShift) + ":" + modificationCandidates;
        }
    }

    public String toString(NumberFormat numberFormat) {

        if (modificationCandidates.isEmpty()) {

            return numberFormat.format(massShift);
        } else {

            return numberFormat.format(massShift) + ":" + modificationCandidates;
        }
    }

    @Override
    public boolean equals(Object o) {

        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        ModificationMatch that = (ModificationMatch) o;

        return Double.compare(that.massShift, massShift) == 0 &&
                modificationCandidates.equals(that.modificationCandidates);
    }

    @Override
    public int hashCode() {

        int result;
        long temp;
        temp = massShift != +0.0d ? Double.doubleToLongBits(massShift) : 0L;
        result = (int) (temp ^ (temp >>> 32));
        result = 31 * result + modificationCandidates.hashCode();
        return result;
    }

    /**
     * Returns the number of modification candidates
     *
     * @return the number of modification candidates
     */
    public int getCandidateCount(){

        return modificationCandidates.size();
    }

    public Modification getModificationCandidate(int index) {

        return modificationCandidates.get(index);
    }
}
