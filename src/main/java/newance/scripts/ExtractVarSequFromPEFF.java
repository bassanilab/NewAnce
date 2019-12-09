package newance.scripts;

import java.io.*;
import java.util.logging.Logger;

/**
 * Created by markusmueller on 19.02.18.
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