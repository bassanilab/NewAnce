package newance.util;

import newance.psmconverter.PeptideMatchData;

import java.util.Set;
import java.util.function.BiFunction;

/**
 * Created by markusmueller on 08.05.18.
 */
public abstract class PsmGrouper implements BiFunction<String, PeptideMatchData, String> {

    public abstract String getMasterGroup();
    public abstract Set<String> getGroups();
}
