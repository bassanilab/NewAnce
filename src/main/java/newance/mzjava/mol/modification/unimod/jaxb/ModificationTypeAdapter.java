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
package newance.mzjava.mol.modification.unimod.jaxb;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import newance.mzjava.mol.Atom;
import newance.mzjava.mol.AtomicSymbol;
import newance.mzjava.mol.PeriodicTable;
import newance.mzjava.mol.modification.NeutralLoss;
import newance.mzjava.mol.modification.unimod.UnimodMod;

import javax.xml.bind.annotation.adapters.XmlAdapter;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;
import java.io.IOException;
import java.io.InputStream;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * XmlAdapter to convert ModT objects into UnimodMod objects.  Because the current unimod.xml name does not have
 * the Interim Name, despite there being a field for it in the style sheet, the unimod_tables.xml file is read to
 * get this data.
 * The unimod_tables.xml file distributed in the mzjava jar has been trimmed to only contain the ex_code_name,
 * record_id, full_name and code_name attributes of the modifications table.
 * <pre>
 * {@code
 * <modifications>
 *      <modifications_row ex_code_name="Acetyl" record_id="1" full_name="Acetylation" code_name="Acetyl"/>
 *      ...
 * </modifications>
 * }
 * </pre>
 *
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class ModificationTypeAdapter extends XmlAdapter<ModT, UnimodMod> {

    public static final Pattern isotopePattern = Pattern.compile("(\\d+)?([A-Z][a-z]?)");

    private final UnimodMod.Builder builder = new UnimodMod.Builder();

    private final TIntObjectMap<String> interimNameMap = new TIntObjectHashMap<>(980); //code_name in xml file

    public ModificationTypeAdapter() {

        try {

            parseUnimodTables();
        } catch (XMLStreamException | IOException e) {

            throw new IllegalStateException(e);
        }
    }

    @Override
    public UnimodMod unmarshal(ModT v) throws Exception {

        TObjectIntMap<Atom> atomCounterMap = extractComposition(v.getDelta());

        Multimap<TObjectIntMap<Atom>, String> neutralLossMap = HashMultimap.create();
        Multimap<TObjectIntMap<Atom>, String> precursorLossMap = HashMultimap.create();
        List<String> sites = new ArrayList<>();
        for(SpecificityT specificity : v.getSpecificity()){

            String site = specificity.getSite();
            sites.add(site);
            for(NeutralLossT neutralLoss : specificity.getNeutralLoss()){

                neutralLossMap.put(extractComposition(neutralLoss), site);
            }
            for(PepNeutralLossT neutralLoss : specificity.getPepNeutralLoss()){

                precursorLossMap.put(extractComposition(neutralLoss), site);
            }
        }

        int recordId = v.getRecordId().intValue();
        builder.setRecordId(recordId).
                setApproved(v.isApproved()).
                setPsiMsName(v.getTitle()).
                setInterimName(interimNameMap.get(recordId)).
                setFullName(v.getFullName()).
                setComposition(atomCounterMap).
                setSites(sites).
                setNeutralLosses(makeLossList(neutralLossMap)).
                setPrecursorNeutralLosses(makeLossList(precursorLossMap)).
                build();

        return builder.build();
    }

    private List<NeutralLoss> makeLossList(Multimap<TObjectIntMap<Atom>, String> lossMap) {

        if(lossMap.isEmpty()) return Collections.emptyList();

        Set<TObjectIntMap<Atom>> keySet = lossMap.keySet();
        List<NeutralLoss> neutralLosses = new ArrayList<>(keySet.size());

        for(TObjectIntMap<Atom> composition : keySet) {

            Set<String> sites = new HashSet<>(lossMap.get(composition));
            neutralLosses.add(new NeutralLoss(composition, 0, sites));
        }

        if (neutralLosses.size() > 1) {
            Collections.sort(neutralLosses, new Comparator<NeutralLoss>() {
                @Override
                public int compare(NeutralLoss o1, NeutralLoss o2) {

                    return Double.compare(o1.getMolecularMass(), o2.getMolecularMass());
                }
            });
        }

        return neutralLosses;
    }

    private TObjectIntMap<Atom> extractComposition(CompositionT compositionT) {

        TObjectIntMap<Atom> atomCounterMap = new TObjectIntHashMap<>();

        for(ElemRefT elementRef : compositionT.getElement()) {

            Matcher matcher = isotopePattern.matcher(elementRef.getSymbol());

            Atom atom;
            if(matcher.matches()) {

                String isotopeCount = matcher.group(1);
                String atomicSymbol = matcher.group(2);
                if (isotopeCount == null) {

                    atom = PeriodicTable.getInstance().getAtom(AtomicSymbol.valueOf(atomicSymbol));
                } else {

                    atom = PeriodicTable.getInstance().getAtom(AtomicSymbol.valueOf(atomicSymbol), Integer.valueOf(isotopeCount));
                }
            } else {

                throw new IllegalStateException("Could not parseAndAdd symbol " + elementRef.getSymbol());
            }

            int count = elementRef.getNumber();
            if (count != 0) {

                atomCounterMap.adjustOrPutValue(atom, count, count);
            }
        }
        return atomCounterMap;
    }

    @Override
    public ModT marshal(UnimodMod v) throws Exception {

        throw new UnsupportedOperationException("Can only read unimod, not write");
    }

    private void parseUnimodTables() throws XMLStreamException, IOException {

        InputStream in = getClass().getResourceAsStream("unimod_tables.xml");

        XMLInputFactory factory = XMLInputFactory.newInstance();
        XMLEventReader reader = factory.createXMLEventReader(in);
        try {
            while (reader.hasNext()) {

                XMLEvent event = reader.peek();
                if (event.isStartElement()) {

                    StartElement element = event.asStartElement();
                    if ("modifications_row".equals(element.getName().getLocalPart())) {
                        readModificationRow(element);
                    }
                }
                reader.nextEvent();
            }
        } catch(XMLStreamException e){

            in.close();
            throw e;
        }finally {

            reader.close();
            in.close();
        }
    }

    private void readModificationRow(StartElement element) {

        int recordId = -1;
        String codeName = "";
        for (Iterator i = element.getAttributes(); i.hasNext(); ) {

            Attribute attr = (Attribute) i.next();
            String localName = attr.getName().getLocalPart();

            if("code_name".equals(localName)) codeName = attr.getValue();
            else if("record_id".equals(localName)) recordId = Integer.parseInt(attr.getValue());
        }

        interimNameMap.put(recordId, codeName);
    }
}
