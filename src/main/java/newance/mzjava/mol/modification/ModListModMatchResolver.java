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

import com.google.common.base.Optional;
import com.google.common.collect.Lists;
import newance.mzjava.mol.Weighable;

import java.util.*;

/**
 * Resolves ModificationMatches against the list of Modifications held by this ModListModMatchResolver.
 *
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class ModListModMatchResolver implements ModificationMatchResolver {

    private final Tolerance tolerance;

    private final Key key = new Key();
    private final List<Modification> modifications;

    public ModListModMatchResolver(Tolerance tolerance, Collection<Modification> mods) {

        this(tolerance, new ArrayList<>(mods));
    }

    public ModListModMatchResolver(Tolerance tolerance, Modification... mods) {

        this(tolerance, Lists.newArrayList(mods));
    }

    private ModListModMatchResolver(Tolerance tolerance, List<Modification> mods) {

        this.tolerance = tolerance;
        modifications = mods;
        Collections.sort(modifications, new Comparator<Modification>() {
            @Override
            public int compare(Modification o1, Modification o2) {

                return Double.compare(o1.getMolecularMass(), o2.getMolecularMass());
            }
        });

        for (int i = 0, modificationsSize = modifications.size() - 1; i < modificationsSize; i++) {

            Modification mod = modifications.get(i);
            Modification nextMod = modifications.get(i + 1);

            if(tolerance.withinTolerance(nextMod.getMolecularMass(), tolerance.getMax(mod.getMolecularMass())))
                throw new IllegalStateException(mod + " and " + nextMod + " masses overlap");
        }
    }

    @Override
    public Optional<Modification> resolve(ModificationMatch modMatch) {

        double mass = modMatch.getMassShift();
        key.setMass(tolerance.getMin(mass));
        int index = Collections.binarySearch(modifications, key, new Comparator<Weighable>(){

            @Override
            public int compare(Weighable w1, Weighable w2) {

                return Double.compare(w1.getMolecularMass(), w2.getMolecularMass());
            }
        });

        if (index < 0) index = -1 * (index + 1);


        if(index > -1 && index < modifications.size()) {

            Modification mod = modifications.get(index);
            return tolerance.withinTolerance(mod.getMolecularMass(), mass) ? Optional.of(mod) : Optional.<Modification>absent();
        } else {

            return Optional.absent();
        }
    }

    private static class Key implements Weighable {

        double mass;

        @Override
        public double getMolecularMass() {

            return mass;
        }

        public void setMass(double mass) {

            this.mass = mass;
        }
    }
}
