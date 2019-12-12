package newance.psmconverter;

import newance.util.RegExpFileFilter;

import java.io.File;
import java.io.IOException;
import java.util.*;
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

public class CometPsmConverter extends PsmConverter {

    public CometPsmConverter(String psmRootDirName, Pattern regex) {

        super(psmRootDirName, regex);
    }

    public void run() throws IOException{

        long start = System.currentTimeMillis();

        System.out.println("Comet psm input path " + psmRootDirName);
        List<File> psmFileList = Arrays.asList(new File(psmRootDirName).listFiles(new RegExpFileFilter(regex)));
        checkState(!psmFileList.isEmpty());

        int nrTasks = psmFileList.size();
        ExecutorService exe = Executors.newFixedThreadPool(nrTasks);

        CountDownLatch latch = new CountDownLatch(nrTasks);
        for (int i = 0; i < nrTasks; i++) {
            exe.submit(new CometPepXmlEntryConverter(psmFileList.get(i), psms, latch));
        }

        try {
            latch.await();

            System.out.println("Number of Comet spectra converted: "+psms.size());
            System.out.println("Comet PepXML conversion ran in " + (System.currentTimeMillis() - start) / 1000d + "s");

        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        exe.shutdown();

    }
}
