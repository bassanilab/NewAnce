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
package newance.mzjava.mol;

import com.google.common.base.Optional;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import gnu.trove.set.TCharSet;
import gnu.trove.set.hash.TCharHashSet;
import newance.mzjava.mol.modification.*;
import newance.mzjava.mol.modification.unimod.UnimodModificationResolver;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Oliver Horlacher
 * @version sqrt -1
 */
class PeptideParser {

    private final TCharSet aaChars = new TCharHashSet();

    private final ModificationResolver defaultModResolver = new CompositeModResolver(
            new UnimodModificationResolver(),
            new NumericModificationResolver()
            );

    private static PeptideParser peptideParser;

    static PeptideParser getInstance(){

        if(peptideParser == null) peptideParser = new PeptideParser();

        return peptideParser;
    }

    private PeptideParser() {

        for(AminoAcid aa : AminoAcid.values()) {

            aaChars.add(aa.getSymbol().charAt(0));
        }
    }

    Peptide parsePeptide(String s) {

        return parsePeptide(s, defaultModResolver);
    }

    Peptide parsePeptide(String s, ModificationResolver modResolver ) {

        List<AminoAcid> residues = new ArrayList<AminoAcid>(s.length());
        Multimap<Integer, Modification> modMap = ArrayListMultimap.create();
        Multimap<ModAttachment, Modification> termModMap = ArrayListMultimap.create();

        parsePeptide(s, residues, modMap, termModMap, modResolver);

        return new Peptide(residues, modMap, termModMap);
    }

    /**
     * Extract the sequence and modification information from <code>s</code>
     *
     * @param s               the string to parse
     * @param sequence        list to hold the AA sequence
     * @param sideChainModMap multi map to hold the side chain modifications
     * @param termModMap      multi map to hold the n-term and c-term modifications
     * @throws PeptideParseException if s is not a valid string
     */
    public void parsePeptide(String s, List<AminoAcid> sequence, Multimap<Integer, Modification> sideChainModMap, Multimap<ModAttachment, Modification> termModMap) {

        parsePeptide(s, sequence, sideChainModMap, termModMap, defaultModResolver);
    }

    /**
     * Extract the sequence and modification information from <code>s</code>
     *
     * @param s               the string to parse
     * @param sequence        list to hold the AA sequence
     * @param sideChainModMap multi map to hold the side chain modifications
     * @param termModMap      multi map to hold the n-term and c-term modifications
     * @param modResolver     function to convert the modification strings found in s to Modification
     * @throws PeptideParseException if s is not a valid string
     */
    void parsePeptide(String s, List<AminoAcid> sequence, Multimap<Integer, Modification> sideChainModMap,
                      Multimap<ModAttachment, Modification> termModMap, ModificationResolver modResolver){

        String input = s;
        s = stripTerminalMods(s, termModMap, modResolver);

        int bracketCount = 0;
        int bracketStartIndex = -1;
        int lastModIndex = -1;
        for(int i = 0; i < s.length(); i++) {

            char c = s.charAt(i);

            if(bracketCount == 0 && aaChars.contains(c)) {

                sequence.add(AminoAcid.valueOf(c));
            } else if(c == '(') {

                if(bracketCount == 0) bracketStartIndex = i;
                bracketCount += 1;
            } else if(c == ')') {

                bracketCount -= 1;
                if(bracketCount == 0) {

                    if(lastModIndex == sequence.size() - 1)  //Prevent A(mod1)(mod2) from being legal
                        throw new PeptideParseException(input + " has illegal modification format");
                    parseMods(s.substring(bracketStartIndex + 1, i), sideChainModMap, sequence.size() - 1, modResolver);
                    lastModIndex = sequence.size() - 1;
                }
            } else if(bracketCount == 0) {

                throw new PeptideParseException(input + " has illegal character " + c);
            }

            if(bracketCount < 0) throw new PeptideParseException(input + " has to many closing brackets" + i);
        }

        if(bracketCount != 0)
            throw new PeptideParseException(input + " has to many opening brackets");

        if(sequence.isEmpty())
            throw new PeptideParseException("The peptide '" + input + "' contains no amino acid residues");
    }

    private String stripTerminalMods(String s, Multimap<ModAttachment, Modification> termModMap, ModificationResolver modResolver) {

        String input = s;
        int length = s.length();
        if(length == 0) return s;

        if(s.charAt(0) == '(') {

            boolean foundNTerm = false;
            for(int i = 1; i < length; i++) {

                if(s.charAt(i) == '_') {

                    try {

                        parseMods(s.substring(1, i - 1), ModAttachment.N_TERM, termModMap, modResolver);
                    } catch (Exception e) {

                        throw new PeptideParseException("Illegal n-term mod format on " + input, e);
                    }
                    foundNTerm = true;
                    s = s.substring(i + 1);
                    break;
                }
            }

            if(!foundNTerm)
                throw new PeptideParseException("Invalid n-term mod format on " + input);
        }

        length = s.length();
        if(length > 2 && s.charAt(length - 1) == ')'){

            for(int i = length - 2; i >= 0; i--) {

                if(s.charAt(i) == '_') {

                    try {

                        parseMods(s.substring(i + 2, length - 1), ModAttachment.C_TERM, termModMap, modResolver);
                    } catch (IllegalStateException e) {

                        throw new PeptideParseException("Illegal c-term mod format on " + s, e);
                    } catch (StringIndexOutOfBoundsException e) {

                        throw new PeptideParseException("Illegal c-term mod format on " + s, e);
                    }
                    s = s.substring(0, i);
                    break;
                }
            }
        }

        return s;
    }

    private void parseMods(String substring, Multimap<Integer, Modification> modMap, Integer index, ModificationResolver modResolver) {

        for(String sMod : substring.split(",\\s?")) {

            Optional<Modification> optional = modResolver.resolve(sMod);
            if (optional.isPresent()) {

                modMap.put(index, optional.get());
            } else {

                throw new PeptideParseException("Cannot resolve modification " + sMod + " on residue " + index);
            }
        }
    }

    private void parseMods(String substring, ModAttachment modAttachment, Multimap<ModAttachment, Modification> termModMap, ModificationResolver modResolver) {

        for(String sMod : substring.split(",\\s?")) {

            Optional<Modification> optional = modResolver.resolve(sMod);
            if (optional.isPresent()) {

                termModMap.put(modAttachment, optional.get());
            } else {

                throw new PeptideParseException("Cannot resolve modification " + sMod);
            }
        }
    }
}
