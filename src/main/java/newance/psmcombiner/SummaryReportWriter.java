package newance.psmcombiner;

import newance.psmconverter.PeptideMatchData;
import newance.util.NewAnceParams;
import newance.util.ProcessPsmUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Created by markusmueller on 03.12.19.
 */
public class SummaryReportWriter {
    private final String fileName;
    private BufferedWriter reportWriter;


    public SummaryReportWriter(String fileName, boolean includeMaxQuant) {

        this.fileName = fileName;
        try {
            reportWriter = new BufferedWriter(new FileWriter(new File(fileName)));
            writeHeader(includeMaxQuant);

        } catch (IOException e) {
            System.out.println("Unable to writeHistograms to file "+ fileName);
            reportWriter = null;
        }
    }

    private void writeHeader(boolean includeMaxQuant) throws IOException {
        reportWriter.write(NewAnceParams.getInstance().toString()+"\n\n");
        if (includeMaxQuant) reportWriter.write("Group\tNrCometPSMs\tNrMaxQuantPSMs\tNrCombinedPSMs\tNrCometPepts\tNrMaxQuantPepts\tNrCombinedPepts\n");
        else reportWriter.write("Group\tNrCometPSMs\tNrCometPepts\n");
    }

    public void write(String group,
                      ConcurrentHashMap<String, List<PeptideMatchData>> combinedPsms,
                      ConcurrentHashMap<String, List<PeptideMatchData>>  cometPsms,
                      ConcurrentHashMap<String, List<PeptideMatchData>>  maxQuantPsms) {

        if (reportWriter==null) return;

        ConcurrentHashMap<String, List<PeptideMatchData>>  noDecoycombinedPsms = ProcessPsmUtils.removeDecoys(combinedPsms);
        ConcurrentHashMap<String, List<PeptideMatchData>>  noDecoyCometPsms = ProcessPsmUtils.removeDecoys(cometPsms);
        ConcurrentHashMap<String, List<PeptideMatchData>>  noDecoyMaxQuantPsms = ProcessPsmUtils.removeDecoys(maxQuantPsms);

        try {
            reportWriter.write(group+"\t"+ProcessPsmUtils.countPsms(noDecoyCometPsms)+"\t"+ ProcessPsmUtils.countPsms(noDecoyMaxQuantPsms)+"\t"+
                    ProcessPsmUtils.countPsms(noDecoycombinedPsms)+"\t"+ProcessPsmUtils.countUniquePeptides(noDecoyCometPsms)+"\t"+
                    ProcessPsmUtils.countUniquePeptides(noDecoyMaxQuantPsms)+"\t"+ ProcessPsmUtils.countUniquePeptides(noDecoycombinedPsms)+"\n");

        } catch (IOException e) {
        }
    }


    public void write(String group, ConcurrentHashMap<String, List<PeptideMatchData>> cometPsms) {

        if (reportWriter==null) return;

        ConcurrentHashMap<String, List<PeptideMatchData>>  noDecoyCometPsms = ProcessPsmUtils.removeDecoys(cometPsms);

        try {
            reportWriter.write(group+"\t"+ProcessPsmUtils.countPsms(noDecoyCometPsms)+"\t"+ProcessPsmUtils.countUniquePeptides(noDecoyCometPsms)+"\n");

        } catch (IOException e) {
        }
    }

    public void close() {
        try {
            if (reportWriter!=null) reportWriter.close();
            reportWriter = null;
        } catch (IOException e) {
            System.out.println("Unable to close file "+ fileName);
        }
    }

}
