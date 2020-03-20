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
import newance.mzjava.mol.Weighable;
import newance.mzjava.mol.modification.Modification;
import newance.mzjava.mol.modification.Tolerance;
import newance.mzjava.mol.modification.unimod.jaxb.UnimodT;
import org.xml.sax.*;
import org.xml.sax.helpers.XMLReaderFactory;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBElement;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;
import javax.xml.transform.sax.SAXSource;
import java.io.InputStream;
import java.util.*;


/**
 * This manager handles modifications registering and look ups from many points
 * of view: the names, the masses and the modified sites (AA + terms).
 *
 * @author nikitin
 * @version 1.0
 */
public class UnimodManager {

    private static UnimodManager INSTANCE;

    private final Map<String, UnimodMod> psiNameToModMap = new HashMap<String, UnimodMod>();
    private final Map<String, UnimodMod> interimNameToModMap = new HashMap<String, UnimodMod>();
    private final List<UnimodMod> modificationList = new ArrayList<UnimodMod>();


    private UnimodManager(){

        parseUnimodXml(getClass().getResourceAsStream("unimod.xml"));
    }


    public static UnimodManager getInstance() {

        if(INSTANCE == null) {

            INSTANCE = new UnimodManager();
        }
        return INSTANCE;
    }


    /**
     *  parse the unimod xml file
     *
     * @param unimodInputStream the path to the unimod.xml file
     * @return
     */
    private void parseUnimodXml(InputStream unimodInputStream){
        try {

            Class<UnimodT> docClass = UnimodT.class;

            String packageName = docClass.getPackage().getName();
            XMLReader reader = XMLReaderFactory.createXMLReader();
            reader.setFeature("http://apache.org/xml/features/allow-java-encodings", true);

            InputSource is = new InputSource(unimodInputStream);
            SAXSource ss = new SAXSource(reader, is);

            JAXBContext jc = JAXBContext.newInstance(packageName);
            Unmarshaller u = jc.createUnmarshaller();
            UnimodT doc = (UnimodT)((JAXBElement)u.unmarshal(ss)).getValue();

            for(UnimodMod mod : doc.getModifications().getMod()) {

                psiNameToModMap.put(mod.getLabel(), mod);
                interimNameToModMap.put(mod.getInterimName(),  mod);
                modificationList.add(mod);
            }

        } catch (JAXBException e) {

            throw new IllegalStateException("Could not read unimod.xml", e);
        } catch (SAXNotSupportedException e) {

            throw new IllegalStateException("Could not read unimod.xml", e);
        } catch (SAXNotRecognizedException e) {

            throw new IllegalStateException("Could not read unimod.xml", e);
        } catch (SAXException e) {

            throw new IllegalStateException("Could not read unimod.xml", e);
        }

        // sort the masses
        Collections.sort(modificationList, new Comparator<Modification>() {
            @Override
            public int compare(Modification o1, Modification o2) {

                return Double.compare(o1.getMolecularMass(), o2.getMolecularMass());
            }
        });
    }


    /**
     * Get the modification for the name.
     *
     * The name can be either a unimod PSI-MS or Interim name.
     *
     * @param name the modification name to get.
     * @return Optional containing the modification if there is a modification with <code>name</code>, Optional.absent otherwise.
     */
    private Optional<UnimodMod> getModificationFromName(String name) {

        UnimodMod mod = psiNameToModMap.get(name);
        if(mod != null) return Optional.of(mod);

        mod = interimNameToModMap.get(name);
        if(mod != null) return Optional.of(mod);

        return Optional.absent();
    }

    public static List<Modification> getModifications(final double mass, final Tolerance tolerance) {

        return getInstance().doGetModifications(mass, tolerance);
    }

    public static List<Modification> getModifications(final double mass, final Tolerance tolerance, final Set<String> aminoAcids) {

        List<Modification> modifs = getInstance().doGetModifications(mass, tolerance);

        List<Modification> matchedModifs = new ArrayList<>();
        for (Modification modif : modifs) {
            UnimodMod uniMod = (UnimodMod) modif;


            for (String aa : aminoAcids) {
                if (uniMod.getSites().contains(aa)) {
                    matchedModifs.add(modif);
                    break;
                }
            }
        }

        return matchedModifs;
    }

    /**
     * Get all modifications at the given mass.
     *
     * @param mass      the modification mass.
     * @param tolerance the precision of the mass.
     * @return the modifications with the given mass.
     */
    private List<Modification> doGetModifications(final double mass, final Tolerance tolerance) {

        List<Modification> modifications = new ArrayList<Modification>();

        Weighable key = new Weighable(){
            @Override
            public double getMolecularMass() {

                return tolerance.getMin(mass);
            }
        };

        int index = Collections.binarySearch(modificationList, key, new Comparator<Weighable>() {
            @Override
            public int compare(Weighable o1, Weighable o2) {

                return Double.compare(o1.getMolecularMass(), o2.getMolecularMass());
            }
        });

        if (index < 0) index = -1 * (index + 1);

        for(; index < modificationList.size(); index++) {

            Modification mod = modificationList.get(index);

            Tolerance.Location location = tolerance.check(mass, mod.getMolecularMass());
            if(location == Tolerance.Location.WITHIN) {
                modifications.add(mod);
            }
            else if(location == Tolerance.Location.LARGER) break;
            else throw new IllegalStateException("Have bug in code");
        }

        return modifications;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();

        for (String key : psiNameToModMap.keySet()) {
            sb.append(key);
            sb.append(" => ");
            sb.append(psiNameToModMap.get(key));
            sb.append("\n");
        }

        return sb.toString();
    }

    public static Optional<UnimodMod> getModification(String modName) {

        return UnimodManager.getInstance().getModificationFromName(modName);
    }

    public static int size() {

        return getInstance().modificationList.size();
    }

    /**
     * Return an unmodifiable list containing all the modifications that are managed by this
     *
     * @return an unmodifiable list containing all the modifications that are managed by this
     */
    public List<UnimodMod> getModificationList() {

        return Collections.unmodifiableList(modificationList);
    }

}
