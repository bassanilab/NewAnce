/**
 * Copyright (C) 2019, SIB/LICR. All rights reserved
 *
 * SIB, Swiss Institute of Bioinformatics
 * Ludwig Institute for Cancer Research (LICR)
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer. Redistributions in binary
 * form must reproduce the above copyright notice, this list of conditions and
 * the following disclaimer in the documentation and/or other materials provided
 * with the distribution. Neither the name of the SIB/LICR nor the names of
 * its contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL SIB/LICR BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

package newance.psmcombiner;

import newance.psmconverter.PeptideMatchData;
import newance.util.PsmGrouper;

import java.util.HashSet;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author Markus Müller
 */

public class RegExpProteinGrouper extends PsmGrouper {
    private final Pattern regExp;
    private final Set<String> exclude;
    private final String masterGroup;
    private final String otherGroup;

    public RegExpProteinGrouper(Pattern regExp, String masterGroup, String otherGroup) {
        this.regExp = regExp;
        this.masterGroup = masterGroup;
        this.otherGroup = otherGroup;
        this.exclude = null;
    }

    public RegExpProteinGrouper(Pattern regExp, Set<String> exclude, String masterGroup, String otherGroup) {
        this.regExp = regExp;
        this.masterGroup = masterGroup;
        this.otherGroup = otherGroup;
        this.exclude = exclude;
    }

    @Override
    public String apply(String specID, PeptideMatchData psm) {

        Set<String> proteins = psm.getProteinAcc();

        for (String protein : proteins) {

            if (exclude!=null) {
                for (String excludedProt : exclude) {
                    if (protein.contains(excludedProt))
                        return otherGroup;
                }
            }

            Matcher m = regExp.matcher(protein);
            if (m.find( ))
                return masterGroup;
        }

        return otherGroup;
    }

    @Override
    public String getMasterGroup() {
        return masterGroup;
    }

    @Override
    public Set<String> getGroups() {

        Set<String> groups = new HashSet<>();
        groups.add(masterGroup);
        groups.add(otherGroup);

        return groups;
    }

}
