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

import com.google.common.base.Optional;
import newance.mzjava.mol.modification.Modification;
import newance.mzjava.mol.modification.ModificationResolver;

import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import static com.google.common.base.Preconditions.checkNotNull;

/**
 * A Modification  that converts Strings to modifications.  The modifications are retrieved from
 * Unimod.
 *
 * This function can be customised by overriding the Modification that is unimod stores for a name by
 * using the (String name, Modification modification).  Or by specifying a translation from a name to
 * a unimod name.
 *
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class UnimodModificationResolver implements ModificationResolver {

    private static final Logger LOGGER = Logger.getLogger(UnimodModificationResolver.class.getName());

    private final Map<String, Modification> overrideMap = new HashMap<String, Modification>();
    private final Map<String, String> translateMap = new HashMap<String, String>();

    public UnimodModificationResolver() {

        translateMap.put("ICAT_light", "ICAT-C");
        translateMap.put("ICAT_heavy", "ICAT-C:13C(9)");
        translateMap.put("AB_old_ICATd0", "ICAT-D");
        translateMap.put("AB_old_ICATd8", "ICAT-D:2H(8)");
        translateMap.put("Deamidation", "Deamidated");
        translateMap.put("Pyro-cmC", "Pyro-carbamidomethyl");
        translateMap.put("Pyro-glu", "Gln->pyro-Glu ");
        translateMap.put("Pyro_glu", "Glu->pyro-Glu");
        translateMap.put("Amide", "Amidated");
        translateMap.put("USM_C_517.202881", "Bodipy");
        translateMap.put("USM_n_230.170762", "TMT6plex");
        translateMap.put("USM_K_357.257892", "TMT6plex");

        overrideMap.put("Propionamide:13C(3)", Modification.parseModification("H5C[13]3NO"));
    }

    @Override
    public Optional<Modification> resolve(String input) {

        if (translateMap.containsKey(input)) input = translateMap.get(input);

        if (overrideMap.containsKey(input)) return Optional.of(overrideMap.get(input));

        Optional<UnimodMod> modificationOp = UnimodManager.getModification(input);
        if (modificationOp.isPresent()) {

            return Optional.<Modification>of(modificationOp.get());
        } else {

            try {

                return Optional.<Modification>of(Modification.parseModification(input));
            } catch (IllegalArgumentException e) {

                LOGGER.log(Level.INFO, "Could not resolve " + input, e);
                return Optional.absent();
            }
        }
    }

    /**
     * Override the modification that the name is mapped to in Unimod.
     *
     * @param name         the modification name
     * @param modification the modification to return
     */
    public void putOverrideUnimod(String name, Modification modification) {

        checkNotNull(name);
        checkNotNull(modification);

        overrideMap.put(name, modification);
    }

    public void clearOverrides(){

        overrideMap.clear();
    }

    /**
     * Add a translation for the name.
     *
     * @param name       the name to translate
     * @param unimodName the unimod name that the name translates to
     */
    public void putTranslate(String name, String unimodName) {

        checkNotNull(name);
        checkNotNull(unimodName);

        translateMap.put(name, unimodName);
    }

    public void clearTranslations(){

        translateMap.clear();
    }
}
