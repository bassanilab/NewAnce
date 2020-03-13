/*
Copyright (C) SIB - Swiss Institute of Bioinformatics, Lausanne, Switzerland
Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
*/


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