/*
Copyright (C) SIB - Swiss Institute of Bioinformatics, Lausanne, Switzerland
Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
*/

package newance.util;

import java.util.function.Function;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author Markus Müller
 */

public class SpectrumFilter implements Function<String, Boolean> {

    private Pattern regExp;
    private final boolean match;

    public SpectrumFilter(Pattern regExp) {
        this.regExp = regExp;
        this.match = true;
    }

    public SpectrumFilter(String regExp) {
        this.regExp = Pattern.compile(regExp);
        this.match = true;
    }

    public SpectrumFilter(Pattern regExp, boolean match) {
        this.regExp = regExp;
        this.match = match;
    }

    public SpectrumFilter(String regExp, boolean match) {
        this.regExp = Pattern.compile(regExp);
        this.match = match;
    }

    @Override
    public Boolean apply(String spectrumID)  {

        try {
            if (regExp == null) return true;

            Matcher m = regExp.matcher(spectrumID);
            if (m.find()) return match;

            return !match;
        } catch (Exception e) {
            return false;
        }
    }

}
