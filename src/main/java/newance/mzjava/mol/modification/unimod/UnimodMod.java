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
package newance.mzjava.mol.modification.unimod;

import gnu.trove.map.TObjectIntMap;
import newance.mzjava.mol.Atom;
import newance.mzjava.mol.Composition;
import newance.mzjava.mol.modification.Modification;
import newance.mzjava.mol.modification.NeutralLoss;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * A extension of the Modification class to store UniMod specific fields.
 * <p>
 * Adds fields for:
 * <ul>
 *     <li>approved</li>
 *     <li>recordId</li>
 *     <li>interim name</li>
 *     <li>full name</li>
 *     <li>sites</li>
 *     <li>neutral losses</li>
 *     <li>precursor neutral losses</li>
 * </ul>
 *
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class UnimodMod extends Modification {

    private final boolean approved;
    private final int recordId;
    private final String interimName;
    private final String fullName;

    private final List<String> sites;
    private final List<NeutralLoss> neutralLosses;
    private final List<NeutralLoss> precursorNeutralLosses;

    protected UnimodMod(int recordId, boolean approved, String psiMsName, String interimName, String fullName,
                        TObjectIntMap<Atom> atomCounterMap, List<String> sites, List<NeutralLoss> neutralLosses, List<NeutralLoss> precursorNeutralLosses) {

        super(psiMsName, new Composition(atomCounterMap, 0));

        this.approved = approved;
        this.recordId = recordId;
        this.interimName = interimName;
        this.fullName = fullName;

        this.sites = Collections.unmodifiableList(new ArrayList<>(sites));
        this.neutralLosses = Collections.unmodifiableList(new ArrayList<>(neutralLosses));
        this.precursorNeutralLosses = Collections.unmodifiableList(new ArrayList<>(precursorNeutralLosses));
    }

    /**
     * Return the list of modification sites
     *
     * @return  Return the list of modification sites
     */
    public List<String> getSites() {

        return sites;
    }

    /**
     * Return an unmodifiable list containing the neutral losses.
     *
     * @return an unmodifiable list containing the neutral losses
     */
    public List<NeutralLoss> getNeutralLosses() {

        return neutralLosses;
    }

    /**
     * Returns an unmodifiable list containing the precursor neutral losses.
     *
     * @return an unmodifiable list containing the precursor neutral losses
     */
    public List<NeutralLoss> getPrecursorNeutralLosses() {

        return precursorNeutralLosses;
    }

    /**
     * Returns the psi ms ms name
     *
     * @return the psi ms ms name
     */
    public String getPsiMsMsName() {

        return getLabel();
    }

    /**
     * Returns the approved flag.
     *
     * @return the approved flag
     */
    public boolean isApproved() {

        return approved;
    }

    /**
     * Returns the UniMod record id.
     *
     * @return the UniMod record id
     */
    public int getRecordId() {

        return recordId;
    }

    /**
     * Returns the UniMod interim name.
     *
     * @return the UniMod interim name
     */
    public String getInterimName() {

        return interimName;
    }

    /**
     * Returns the UniMod full name.
     *
     * @return the UniMod full name
     */
    public String getFullName() {

        return fullName;
    }

    /**
     * Returns the Composition associated with this UnimodMod.
     *
     * @return the Composition associated with this UnimodMod
     */
    public Composition getComposition() {

        return (Composition)getMass();
    }

    public static final class Builder {

        private boolean approved;
        private int recordId;
        private String interimName;
        private String psiMsName;
        private String fullName;
        private TObjectIntMap<Atom> atomCountMap;
        private List<String > sites;
        private List<NeutralLoss> neutralLosses;
        private List<NeutralLoss> precursorNeutralLosses;

        /**
         * Set the approved flag
         *
         * @param approved the approved flag
         * @return the builder
         */
        public Builder setApproved(boolean approved) {

            this.approved = approved;
            return this;
        }

        /**
         * Set the record id
         *
         * @param recordId the record id
         * @return the builder
         */
        public Builder setRecordId(int recordId) {

            this.recordId = recordId;
            return this;
        }

        /**
         * Set the interim name.
         *
         * @param interimName the interim name
         * @return the builder
         */
        public Builder setInterimName(String interimName) {

            this.interimName = interimName;
            return this;
        }

        /**
         * Set the psi ms name
         *
         * @param psiMsName the psi ms name
         * @return the builder
         */
        public Builder setPsiMsName(String psiMsName) {

            this.psiMsName = psiMsName;
            return this;
        }

        /**
         * Set the full name.
         *
         * @param fullName the full name
         * @return the builder
         */
        public Builder setFullName(String fullName) {

            this.fullName = fullName;
            return this;
        }

        /**
         * Set the atom counter map that contains the composition
         *
         * @param atomCountMap the atom counter map
         * @return the builder
         */
        public Builder setComposition(TObjectIntMap<Atom> atomCountMap) {

            this.atomCountMap = atomCountMap;
            return this;
        }

        /**
         * Set the sites list
         *
         * @param sites the list of modification sites
         * @return the builder
         */
        public Builder setSites(List<String> sites) {

            this.sites = sites;
            return this;
        }

        /**
         * Set the neutral loss list
         *
         * @param neutralLosses the neutral loss list
         * @return the builder
         */
        public Builder setNeutralLosses(List<NeutralLoss> neutralLosses) {

            this.neutralLosses = neutralLosses;
            return this;
        }

        /**
         * Set the precursor neutral loss list
         *
         * @param precursorNeutralLosses the presursor neutral loss list
         * @return the builder
         */
        public Builder setPrecursorNeutralLosses(List<NeutralLoss> precursorNeutralLosses) {

            this.precursorNeutralLosses = precursorNeutralLosses;
            return this;
        }

        /**
         * Build the UnimodMod.
         *
         * @return the new UnimodMod
         */
        public UnimodMod build() {

            return new UnimodMod(recordId, approved, psiMsName, interimName, fullName, atomCountMap,
                    sites != null ? sites : Collections.<String>emptyList(),
                    neutralLosses != null ? neutralLosses : Collections.<NeutralLoss>emptyList(),
                    precursorNeutralLosses!= null ? precursorNeutralLosses : Collections.<NeutralLoss>emptyList());
        }

    }
}
