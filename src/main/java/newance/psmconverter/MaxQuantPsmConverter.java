package newance.psmconverter;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.regex.Pattern;

import static com.google.common.base.Preconditions.checkState;

/**
 * Copyright (C) 2019
 * @author Markus Müller
 * @Institutions: SIB, Swiss Institute of Bioinformatics; Ludwig Institute for Cancer Research
 */

public class MaxQuantPsmConverter extends PsmConverter {


    public MaxQuantPsmConverter(String psmRootDirName, Pattern regex) {

        super(psmRootDirName, regex);    }


    public void run() throws IOException{

        long start = System.currentTimeMillis();

        System.out.println("MaxQuant psm input path " + psmRootDirName);

        final List<File> psmFileList = new ArrayList<>();
        Files.walk(Paths.get(psmRootDirName))
                .filter(Files::isRegularFile)
                .forEach((f)->{
                    if( f.toString().endsWith("msms.txt")) psmFileList.add(f.toFile());
                });

        checkState(!psmFileList.isEmpty());

        int nrTasks = psmFileList.size();
        ExecutorService exe = Executors.newFixedThreadPool(nrTasks);

        CountDownLatch latch = new CountDownLatch(nrTasks);
        for (int i = 0; i < nrTasks; i++) {
            exe.submit(new MaxQuantMSMSEntryConverter(psmFileList.get(i), psms, latch));
        }

        try {
            latch.await();

            System.out.println("Number of MaxQuant spectra converted: "+psms.size());
            System.out.println("MaxQuant msms.txt conversion ran in " + (System.currentTimeMillis() - start) / 1000d + "s");

        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        exe.shutdown();
    }
}
