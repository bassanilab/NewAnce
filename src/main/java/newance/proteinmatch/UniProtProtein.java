package newance.proteinmatch;

import java.io.Serializable;

/**
 * Copyright (C) 2019
 * @author Markus Müller
 * @Institutions: SIB, Swiss Institute of Bioinformatics; Ludwig Institute for Cancer Research
 */

public class UniProtProtein implements Serializable{
    final protected String uniProtAC;
    final protected String uniProtName;
    final protected String description;
    final protected char[] sequence;
    final protected String geneName;
    final protected String dbFlag;


    public UniProtProtein(String uniProtAC, String uniProtName, String description, String geneName, String dbFlag, String sequence) {
        this.uniProtAC = uniProtAC;
        this.uniProtName = uniProtName;
        this.description = description;
        this.geneName = geneName;
        this.dbFlag = dbFlag;
        this.sequence = sequence.toCharArray();
    }

    public String getUniProtAC() {
        return uniProtAC;
    }

    public String getUniProtName() {
        return uniProtName;
    }

    public String getDescription() {
        return description;
    }

    public String getSequence() {
        return new String(sequence);
    }

    public char[] getCharArray() { return sequence; }

    public String getGeneName() {
        return geneName;
    }

    public String getDbFlag() {
        return dbFlag;
    }

    public String getFastaUniProtName() {
        return dbFlag+"|"+uniProtAC+"|"+uniProtName;
    }
}
