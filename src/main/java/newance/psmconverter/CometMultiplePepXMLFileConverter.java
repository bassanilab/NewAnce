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

package newance.psmconverter;

import newance.psmcombiner.GroupedFDRCalculator;
import newance.util.NewAnceParams;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.regex.Pattern;

import static com.google.common.base.Preconditions.checkState;

/**
 * @author Markus Müller
 */

public class CometMultiplePepXMLFileConverter extends MultiplePsmFileConverter {

    private final GroupedFDRCalculator groupedFDRCalculator;
    private final boolean reportHistosOnly;

    public CometMultiplePepXMLFileConverter(String psmRootDirName, Pattern regex) {

        super(psmRootDirName, regex);

        this.groupedFDRCalculator = null;
        this.reportHistosOnly = false;
    }

    public CometMultiplePepXMLFileConverter(String psmRootDirName, Pattern regex, GroupedFDRCalculator groupedFDRCalculator, boolean reportHistosOnly) {

        super(psmRootDirName, regex);

        this.groupedFDRCalculator = groupedFDRCalculator;
        this.reportHistosOnly = reportHistosOnly;
    }

    public void run() throws IOException{

        long start = System.currentTimeMillis();

        System.out.println("Comet psm input path " + psmRootDirName);

        final List<File> psmFileList = new ArrayList<>();
        Files.walk(Paths.get(psmRootDirName))
                .filter(f -> regex.matcher(f.getFileName().toString()).find())
                .forEach(f -> psmFileList.add(f.toFile()));
        checkState(!psmFileList.isEmpty());

        int nrTasks = psmFileList.size();
        int maxNrThreads = NewAnceParams.getInstance().getNrThreads();

        int nrIter = nrTasks/maxNrThreads;
        int threadCnt = 0;
        for (int i=0;i<=nrIter;i++) {

            final ConcurrentHashMap<String,List<PeptideSpectrumMatch>> psmBuffer = new ConcurrentHashMap<>();
            int nrThreads = (nrTasks-threadCnt>maxNrThreads)?maxNrThreads:nrTasks-threadCnt;
            if (nrThreads<=0) break;

            ExecutorService exe = Executors.newFixedThreadPool(nrThreads);

            CountDownLatch latch = new CountDownLatch(nrThreads);
            for (int j =0;j<nrThreads;j++) {
                exe.submit(new CometPepXmlConverter(psmFileList.get(threadCnt), psmBuffer, latch));
                threadCnt++;
            }

            try {
                latch.await();

                if (groupedFDRCalculator!=null) groupedFDRCalculator.addAll(psmBuffer);
                if (!reportHistosOnly) psms.putAll(psmBuffer);

            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            exe.shutdown();
        }

        System.out.println("Number of Comet spectra converted: "+psms.size());
        System.out.println("Comet PepXML conversion ran in " + (System.currentTimeMillis() - start) / 1000d + "s");

    }

}
