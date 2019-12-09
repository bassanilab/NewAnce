package newance.util;

import java.util.function.Function;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created by markusmueller on 06.04.17.
 */
public class SpectrumFilter implements Function<String, Boolean> {

    private Pattern regExp;
    private final boolean match;

    public SpectrumFilter(Pattern regExp) {
        this.regExp = regExp;
        this.match = true;
    }

    public SpectrumFilter(String regExp) {
        this.regExp = Pattern.compile(regExp);
        this.match = true;
    }

    public SpectrumFilter(Pattern regExp, boolean match) {
        this.regExp = regExp;
        this.match = match;
    }

    public SpectrumFilter(String regExp, boolean match) {
        this.regExp = Pattern.compile(regExp);
        this.match = match;
    }

    @Override
    public Boolean apply(String spectrumID)  {

        try {
            if (regExp == null) return true;

            Matcher m = regExp.matcher(spectrumID);
            if (m.find()) return match;

            return !match;
        } catch (Exception e) {
            return false;
        }
    }

}
