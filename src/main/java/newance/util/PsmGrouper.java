package newance.util;

import newance.psmconverter.PeptideMatchData;

import java.util.Set;
import java.util.function.BiFunction;

/**
 * Copyright (C) 2019
 * @author Markus Müller
 * @Institutions: SIB, Swiss Institute of Bioinformatics; Ludwig Institute for Cancer Research
 */

public abstract class PsmGrouper implements BiFunction<String, PeptideMatchData, String> {

    public abstract String getMasterGroup();
    public abstract Set<String> getGroups();
}
