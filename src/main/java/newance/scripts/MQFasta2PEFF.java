package newance.scripts;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * Created by markusmueller on 19.02.18.
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