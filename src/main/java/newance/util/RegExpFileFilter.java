package newance.util;

import java.io.File;
import java.io.FileFilter;
import java.util.regex.Pattern;

/**
 * Copyright (C) 2019
 * @author Markus Müller
 * @Institutions: SIB, Swiss Institute of Bioinformatics; Ludwig Institute for Cancer Research
 */

public class RegExpFileFilter implements FileFilter {

    private final Pattern regex;

    public RegExpFileFilter(Pattern regex){
        this.regex = regex;
    }

    @Override
    public boolean accept(File pathname) {
        return  (regex.matcher(pathname.getAbsolutePath()).find());
    }
}

