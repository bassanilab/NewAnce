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

import newance.mzjava.mol.modification.AbsoluteTolerance;
import newance.mzjava.mol.modification.ModListModMatchResolver;
import newance.mzjava.mol.modification.Modification;
import newance.psmcombiner.CometPsm2PsqlStringFunction;
import newance.psmcombiner.CometPsm2StringFunction;
import newance.psmcombiner.GroupedFDRCalculator;
import newance.psmcombiner.Psm2StringFunction;
import newance.util.NewAnceParams;
import newance.util.PsmPredicate;
import newance.util.StringFileWriter;
import org.junit.Ignore;
import org.junit.Test;

import javax.xml.bind.JAXBException;
import javax.xml.stream.XMLStreamException;
import java.io.File;
import java.io.IOException;
import java.nio.file.*;
import java.util.*;
import java.util.regex.Pattern;

/**
 * Created by markusmueller on 27.01.20.
 */
public class PepXMLParserTest {

    @Test
    @Ignore
    public void testPeptXMLParser() throws IOException, XMLStreamException, JAXBException {

        String psmFile = "/Users/markusmueller/Documents/data/NewAnce/0D5P/lncRNA/Comet/mgf/20170913_QEh1_LC1_CHC_SA_HLAIp_OD5P_DAC_1_R1.pep.xml";

        final Map<String, List<PeptideSpectrumMatch>> psmMap = new HashMap<>();

        NewAnceParams params = NewAnceParams.getInstance();

        PsmPredicate psmPredicate = new PsmPredicate(params);

        PeptideSpectrumMatchList peptideSpectrumMatchList = new PeptideSpectrumMatchList(new SpectrumKeyFunctionImpl(), psmPredicate, psmMap);

        Collection<Modification> modifications = params.getModifications();
        ModListModMatchResolver modMatchResolver = new ModListModMatchResolver(new AbsoluteTolerance(params.getModifMatchMassTol()), modifications);

        CometPEFFPepXmlReader psmReader = new CometPEFFPepXmlReader( null, true, modMatchResolver);
        psmReader.parse(new File(psmFile), peptideSpectrumMatchList);

        final Psm2StringFunction stringFunction = new CometPsm2StringFunction(null,null);

        for (String specID : psmMap.keySet()) {
            for (PeptideSpectrumMatch psm : psmMap.get(specID)) {
                if (psm.isVariant()) {
                    System.out.println(stringFunction.apply(specID, psmMap.get(specID)));
                }
            }

        }
    }

    @Test
    @Ignore
    public void testMaxQuantParser() {
        String msmsFile = "/Users/markusmueller/Documents/data/0D5P/lncRNA/complete/MaxQuant/msms.txt";

        final Map<String, List<PeptideSpectrumMatch>> psmMap = new HashMap<>();

        NewAnceParams params = NewAnceParams.getInstance();

        PsmPredicate psmPredicate = new PsmPredicate(params);

        PeptideSpectrumMatchList peptideSpectrumMatchList = new PeptideSpectrumMatchList(new SpectrumKeyFunctionImpl(), psmPredicate, psmMap);

        MaxQuantPsmReader psmReader = new MaxQuantPsmReader();
        psmReader.parse(new File(msmsFile), peptideSpectrumMatchList);

    }

    @Test
    @Ignore
    public void testListFiles() throws IOException {
        File dir = new File("/Users/markusmueller/Documents/data/0D5P/lncRNA/complete");

        final Pattern pattern = Pattern.compile(".*txt$");
        Files.walk(Paths.get(dir.getAbsolutePath()))
                .filter(f -> pattern.matcher(f.getFileName().toString()).find())
                .forEach(System.out::println);

        dir = new File("/Users/markusmueller/Documents/data/0D5P/lncRNA/complete/MaxQuant");

        Path path = FileSystems.getDefault().getPath(dir.getAbsolutePath()+"/msms.txt");
        System.out.println(path.getFileName());

        Files.walk(Paths.get(dir.getAbsolutePath()))
                .filter(Files::isRegularFile)
                .filter(f -> f.getFileName().equals("msms.txt"))
                .forEach(f -> System.out.println(f.getFileName()));
    }
}
