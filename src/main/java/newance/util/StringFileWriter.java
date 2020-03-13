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

import newance.psmcombiner.Psm2StringFunction;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;
import java.util.function.Consumer;

/**
 * @author Markus Müller
 */

public class StringFileWriter implements Consumer<String> {

    private final String fileName;
    private BufferedWriter writer;
    private final Psm2StringFunction stringFunction;
    private final Set<String> writtenStrings;
    private final boolean uniqueStrings;

    public StringFileWriter(String filename, Psm2StringFunction stringFunction) {
        this.fileName = filename;
        this.stringFunction = stringFunction;
        this.writtenStrings = new HashSet<>();
        this.uniqueStrings = false;

        try {
            writer = new BufferedWriter(new FileWriter(new File(fileName)));
            String header = stringFunction.getHeader();
            if (!header.isEmpty()) writer.write(header+"\n");
        }
        catch (IOException e) {
            writer = null;
            System.out.println("Cannot writeHistograms to file "+filename+"!");
        }
    }


    public StringFileWriter(String filename, Psm2StringFunction stringFunction, boolean uniqueStrings) {
        this.fileName = filename;
        this.stringFunction = stringFunction;
        this.writtenStrings = new HashSet<>();
        this.uniqueStrings = uniqueStrings;

        try {
            writer = new BufferedWriter(new FileWriter(new File(fileName)));
            String header = stringFunction.getHeader();
            if (!header.isEmpty()) writer.write(header+"\n");
        }
        catch (IOException e) {
            writer = null;
            System.out.println("Cannot writeHistograms to file "+filename+"!");
        }
    }

    @Override
    public void accept(String s) {

        if (uniqueStrings) {
            if (writtenStrings.contains(s)) return;
        }

        try {
            if (writer!=null) writer.write(s+"\n");
        }
        catch (IOException e) {
            System.out.println("Cannot writeHistograms to file " + fileName);
        }
    }

    public void flush() {
        try {
            if (writer!=null) writer.flush();
        }
        catch (IOException e) {
            System.out.println("Cannot flush file " + fileName);
        }
    }

    public void close() {
        try {
            if (writer!=null) writer.close();
        }
        catch (IOException e) {
            System.out.println("Cannot close file " + fileName);
        }
    }
}
