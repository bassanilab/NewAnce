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
package newance.scripts;

import java.io.*;
import java.util.logging.Logger;

/**
 * @author Markus Müller
 */

public class ExtractVarSequFromPEFF {

    private static final Logger LOGGER = Logger.getLogger(ExtractVarSequFromPEFF.class.getName());

    protected final File peffFastaFileIn;
    protected final File peffFastaFileOut;

    public ExtractVarSequFromPEFF(String peffFastaFileInName, String peffFastaFileOutName) {

        peffFastaFileIn = new File(peffFastaFileInName);
        peffFastaFileOut = new File(peffFastaFileOutName);

        transform();
    }

    protected void transform() {

        try {
            BufferedReader lineReader = new BufferedReader(new FileReader(peffFastaFileIn));
            BufferedWriter writer = new BufferedWriter(new FileWriter(peffFastaFileOut));

            String line;
            boolean doCopy = false;
            while ((line = lineReader.readLine()) != null) {

                if (line.isEmpty()) continue;

                if (line.startsWith(">")) {
                     doCopy = line.contains("VariantSimple");
                }

                if (doCopy) {
                    writer.write(line+"\n");
                    writer.flush();
                }
            }
            lineReader.close();
            writer.close();
        } catch (FileNotFoundException e) {
            LOGGER.severe(e.getMessage());
        } catch (IOException e) {
            LOGGER.severe(e.getMessage());
        }

    }


    public static void main(String[] args) {

        ExtractVarSequFromPEFF extractor = new ExtractVarSequFromPEFF(args[0],args[1]);
        extractor.transform();
    }
}