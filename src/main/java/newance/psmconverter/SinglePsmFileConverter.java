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
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CountDownLatch;

/**
 * @author Markus Müller
 */

public abstract class SinglePsmFileConverter implements Runnable {

    protected final File psmFile;
    protected final NewAnceParams params;
    protected final Map<String,List<PeptideSpectrumMatch>> psms;
    protected final CountDownLatch latch;

    public SinglePsmFileConverter(File psmFile, Map<String,List<PeptideSpectrumMatch>> psms, CountDownLatch latch) {

        this.params = NewAnceParams.getInstance();
        this.psmFile = psmFile;
        this.psms = psms;
        this.latch = latch;
    }

    public abstract void run();

    protected void addPsms(Map<String, List<PeptideSpectrumMatch>> psmMap) {

        for (String spectrumID : psmMap.keySet()) {

            List<PeptideSpectrumMatch> matches = psmMap.get(spectrumID);

            if (psms.containsKey(spectrumID)) {
                psms.get(spectrumID).addAll(matches);
            } else {
                psms.put(spectrumID,Collections.synchronizedList(matches));
            }
        }
    }
}
