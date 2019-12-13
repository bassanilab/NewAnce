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
