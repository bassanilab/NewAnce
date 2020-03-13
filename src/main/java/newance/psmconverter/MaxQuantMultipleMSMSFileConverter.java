/*
Copyright (C) SIB - Swiss Institute of Bioinformatics, Lausanne, Switzerland
Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
*/

package newance.psmconverter;

import newance.util.NewAnceParams;

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
 * @author Markus Müller
 */

public class MaxQuantMultipleMSMSFileConverter extends MultiplePsmFileConverter {


    public MaxQuantMultipleMSMSFileConverter(String psmRootDirName, Pattern regex) {

        super(psmRootDirName, regex);    }


    public void run() throws IOException{

        long start = System.currentTimeMillis();

        System.out.println("MaxQuant psm input path " + psmRootDirName);

        final List<File> psmFileList = new ArrayList<>();
        Files.walk(Paths.get(psmRootDirName))
                .filter(Files::isRegularFile)
                .filter(f -> regex.matcher(f.getFileName().toString()).find())
                .forEach(f -> psmFileList.add(f.toFile()));

        checkState(!psmFileList.isEmpty());

        int nrTasks = psmFileList.size();
        int maxNrThreads = NewAnceParams.getInstance().getNrThreads();

        int nrIter = nrTasks/maxNrThreads;
        int threadCnt = 0;
        for (int i=0;i<=nrIter;i++) {

            int nrThreads = (nrTasks-threadCnt>maxNrThreads)?maxNrThreads:nrTasks-threadCnt;
            if (nrThreads<=0) break;

            ExecutorService exe = Executors.newFixedThreadPool(nrThreads);

            CountDownLatch latch = new CountDownLatch(nrThreads);
            for (int j =0;j<nrThreads;j++) {
                exe.submit(new MaxQuantMSMSConverter(psmFileList.get(threadCnt), psms, latch));
                threadCnt++;
            }

            try {
                latch.await();

            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            exe.shutdown();
        }

        System.out.println("Number of MaxQuant spectra converted: " + psms.size());
        System.out.println("MaxQuant msms.txt conversion ran in " + (System.currentTimeMillis() - start) / 1000d + "s");

    }
}
