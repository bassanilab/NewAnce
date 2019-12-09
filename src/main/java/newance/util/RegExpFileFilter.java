package newance.util;

import java.io.File;
import java.io.FileFilter;
import java.util.regex.Pattern;

/**
 * @author Markus Muller
 * @version 0.0
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

