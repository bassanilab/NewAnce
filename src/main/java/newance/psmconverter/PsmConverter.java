package newance.psmconverter;

import newance.proteinmatch.UniProtDB;

import java.io.IOException;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.regex.Pattern;

/**
 * @author Markus Muller
 */
public abstract class PsmConverter {

    protected final String psmRootDirName;
    protected final Pattern regex;
    protected final ConcurrentHashMap<String,List<PeptideMatchData>> psms;

    public PsmConverter(String psmRootDirName, Pattern regex) {

        this.psmRootDirName = psmRootDirName;
        this.regex = regex;
        this.psms = new ConcurrentHashMap<>();
    }


    public abstract void run() throws IOException;

    public void addDBProteins(UniProtDB uniProtDB) {

        if (uniProtDB==null) return;

        psms.forEach(100000,new AddUniProtIds(uniProtDB));
    }

    public void reportDBProteins() {

        psms.forEach((id,psm) -> System.out.println(id+", "+psm.get(0).getPeptide().toString()+", "+psm.get(0).getProteinAcc().toString()));
    }

    public ConcurrentHashMap<String, List<PeptideMatchData>> getPsms() {
        return psms;
    }
}
