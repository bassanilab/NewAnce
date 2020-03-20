package newance.mzjava.ms.io.mgf;

import newance.mzjava.ms.spectrum.MsnSpectrum;

/**
 * This interface has to be implemented by parsers that can extract information from
 * the MGF {@code TITLE} line.
 *
 * @version 1.0
 */
public interface TitleParser {

    boolean parseTitle(String title, MsnSpectrum spectrum);
}
