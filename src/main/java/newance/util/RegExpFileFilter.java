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

import java.io.File;
import java.io.FileFilter;
import java.util.regex.Pattern;

/**
 * @author Markus Müller
 */

public class RegExpFileFilter implements FileFilter {

    private final Pattern regex;

    public RegExpFileFilter(Pattern regex){
        this.regex = regex;
    }

    @Override
    public boolean accept(File pathname) {
        return  (regex.matcher(pathname.getAbsolutePath()).find());
    }
}

