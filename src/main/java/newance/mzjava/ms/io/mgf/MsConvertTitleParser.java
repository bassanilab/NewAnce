/*
Copyright (C) SIB - Swiss Institute of Bioinformatics, Lausanne, Switzerland
Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
*/

package newance.mzjava.ms.io.mgf;

import newance.mzjava.ms.spectrum.MsnSpectrum;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created by markusmueller on 20.03.20.
 */
public class MsConvertTitleParser implements TitleParser {

    private final Pattern pattern = Pattern.compile("^([^\\.]+)\\.(\\d+)\\.(\\d+)\\.(\\d+)");

    @Override
    public boolean parseTitle(String title, MsnSpectrum spectrum) {

        Matcher matcher = pattern.matcher(title);

        if (matcher.find()) {

            spectrum.addScanNumber(Integer.parseInt(matcher.group(2)));
            return true;
        } else {

            return false;
        }

    }

}

