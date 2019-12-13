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
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * @author Markus Müller
 */

public class MQFasta2PEFF {

    private class SAAV {
        private final int position;
        private final char newAA;
        private final char oldAA;

        public SAAV(int position, char newAA, char oldAA) {
            this.position = position;
            this.newAA = newAA;
            this.oldAA = oldAA;
        }

        //\VariantSimple=(21|L)(68|R)
        public String toPEFFString() {

            String peffVar = "(" + position + "|" + newAA + ")";
            return peffVar;
        }
    }

    private static final Logger LOGGER = Logger.getLogger(MQFasta2PEFF.class.getName());

    protected final File mqFastaFile;
    protected final File peffFastaFile;

    public MQFasta2PEFF(String mqFastaFileName, String peffFastaFileName) {

        mqFastaFile = new File(mqFastaFileName);
        peffFastaFile = new File(peffFastaFileName);

        transform();
    }

    protected void transform() {

        try {
            BufferedReader lineReader = new BufferedReader(new FileReader(mqFastaFile));
            BufferedWriter writer = new BufferedWriter(new FileWriter(peffFastaFile));

            String line;
            String seq = "";
            while ((line = lineReader.readLine()) != null) {

                if (line.isEmpty()) continue;

                if (line.startsWith(">")) {
                    writer.flush();
                    List<String> headerFields = parseHeader(line);
                    writer.write(headerFields.get(0));
                    if (headerFields.size() > 1) {
                        writer.write(" \\VariantSimple=");
                        for (int i = 1; i < headerFields.size(); i++) {
                            writer.write(headerFields.get(i));
                        }
                    }
                    writer.write("\n");
                } else {
                    writer.write(line + "\n");
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

    //>ENSP00000390635.1      rs139896904_0:E171K
    //>ENSP00000428660.1      rs7418389_0:F72L;rs12755088_0:M102T
    protected List<String> parseHeader(String header) {

        String[] fields = header.split("\\t");

        List<String> headerFields = new ArrayList<>();
        headerFields.add(fields[0]);

        if (fields.length > 1) {
            String[] variants = header.split(";");
            for (int i = 0; i < variants.length; i++) {
                String variant = parseMQVariant(variants[i]);
                if (!variant.isEmpty())
                    headerFields.add(variant);
            }

        }

        return headerFields;
    }

    protected String parseMQVariant(String variant) {

        String[] fields = variant.split(":");

        int position = 0;
        try {
            position = Integer.parseInt(fields[1].substring(1, fields[1].length() - 1));
        } catch(java.lang.NumberFormatException e) {
           return "";
        }

        SAAV saav = new SAAV(position,fields[1].charAt(fields[1].length() - 1),fields[1].charAt(0));

        return saav.toPEFFString();
    }

    public static void main(String[] args) {

        MQFasta2PEFF mqFasta2PEFF = new MQFasta2PEFF(args[0],args[1]);
        mqFasta2PEFF.transform();
    }
}