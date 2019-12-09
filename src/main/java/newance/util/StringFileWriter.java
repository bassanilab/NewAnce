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
 * Created by markusmueller on 03.12.19.
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
