/*
Copyright (C) SIB - Swiss Institute of Bioinformatics, Lausanne, Switzerland
Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
*/

package newance.psmcombiner;

import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import newance.mzjava.mol.Peptide;
import newance.psmconverter.PeptideSpectrumMatch;
import newance.util.NewAnceParams;
import org.junit.Assert;
import org.junit.Test;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

/**
 * @author Markus Müller
 */

public class CometMaxQuantPsmMergerTest {

    @Test
    public void test_merge() {

        ConcurrentHashMap<String, List<PeptideSpectrumMatch>> cometPsms = getCometPsms();
        ConcurrentHashMap<String, List<PeptideSpectrumMatch>> mqPsms = getMQPsms();
        ConcurrentHashMap<String, List<PeptideSpectrumMatch>> combined = new ConcurrentHashMap<>();

        CometMaxQuantPsmMerger combiner = new CometMaxQuantPsmMerger(mqPsms,combined);
        cometPsms.forEach(combiner);

        Assert.assertTrue(combined.size()==1);
        Assert.assertTrue(combined.get("spec1").size()==2);
        Assert.assertEquals("PEPTIDE",combined.get("spec1").get(0).getPeptide().toString());
        Assert.assertEquals("PEPTIDER",combined.get("spec1").get(1).getPeptide().toString());

    }

    public static ConcurrentHashMap<String, List<PeptideSpectrumMatch>> getCometPsms()
    {

        ConcurrentHashMap<String, List<PeptideSpectrumMatch>> psms = new ConcurrentHashMap<>();

        List<String> prots = new ArrayList<>();
        prots.add("protein1");
        prots.add("protein2");

        TObjectDoubleMap<String> scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("xcorr",1.0);
        scoreMap.put("deltacn",0.1);
        scoreMap.put("spscore",100.0);

        PeptideSpectrumMatch peptideSpectrumMatch1 = new PeptideSpectrumMatch("spectrumFile", Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, false, false, null);

        scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("xcorr",2.0);
        scoreMap.put("deltacn",0.2);
        scoreMap.put("spscore",200.0);

        PeptideSpectrumMatch peptideSpectrumMatch2 = new PeptideSpectrumMatch("spectrumFile",Peptide.parse("PEPTIDER"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, false, false, null);

        psms.put("spec1",Collections.synchronizedList(new ArrayList<>()));
        psms.get("spec1").add(peptideSpectrumMatch1);
        psms.get("spec1").add(peptideSpectrumMatch2);

        scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("xcorr",3.0);
        scoreMap.put("deltacn",0.3);
        scoreMap.put("spscore",300.0);

        peptideSpectrumMatch1 = new PeptideSpectrumMatch("spectrumFile",Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
            100, 101, 1001.1, false, false, null);
        psms.put("spec2",Collections.synchronizedList(new ArrayList<>()));
        psms.get("spec2").add(peptideSpectrumMatch1);

        peptideSpectrumMatch1 = new PeptideSpectrumMatch("spectrumFile",Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, false, false, null);
        psms.put("spec5",Collections.synchronizedList(new ArrayList<>()));
        psms.get("spec5").add(peptideSpectrumMatch1);

        return psms;
    }

    public static ConcurrentHashMap<String, List<PeptideSpectrumMatch>> getMQPsms()
    {

        ConcurrentHashMap<String, List<PeptideSpectrumMatch>> psms = new ConcurrentHashMap<>();

        List<String> prots = new ArrayList<>();
        prots.add("protein1");
        prots.add("protein2");

        NewAnceParams params = NewAnceParams.getInstance();


        TObjectDoubleMap<String> scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("Score",1.0);

        PeptideSpectrumMatch peptideSpectrumMatch1 = new PeptideSpectrumMatch("spectrumFile",Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, false, false, null);

        scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("Score",2.0);

        PeptideSpectrumMatch peptideSpectrumMatch2 = new PeptideSpectrumMatch("spectrumFile",Peptide.parse("PEPTIDER"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, false, false, null);

        psms.put("spec1",Collections.synchronizedList(new ArrayList<>()));
        psms.get("spec1").add(peptideSpectrumMatch1);
        psms.get("spec1").add(peptideSpectrumMatch2);

        scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("Score",3.0);

        peptideSpectrumMatch1 = new PeptideSpectrumMatch("spectrumFile",Peptide.parse("PEPTIDER"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, false, false, null);
        psms.put("spec2",Collections.synchronizedList(new ArrayList<>()));
        psms.get("spec2").add(peptideSpectrumMatch1);

        scoreMap = new TObjectDoubleHashMap<>();
        scoreMap.put("Score",4.0);

        peptideSpectrumMatch1 = new PeptideSpectrumMatch("spectrumFile",Peptide.parse("PEPTIDE"), prots, scoreMap, 1, 1,
                100, 101, 1001.1, false, false, null);
        psms.put("spec3",Collections.synchronizedList(new ArrayList<>()));
        psms.get("spec3").add(peptideSpectrumMatch1);

        return psms;
    }
}
