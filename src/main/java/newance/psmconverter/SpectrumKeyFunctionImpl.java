package newance.psmconverter;

import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.expasy.mzjava.proteomics.ms.ident.SpectrumIdentifier;
import newance.util.SpectrumKeyFunction;

import java.util.regex.Pattern;

/**
 * Copyright (C) 2019
 * @author Markus Müller
 * @Institutions: SIB, Swiss Institute of Bioinformatics; Ludwig Institute for Cancer Research
 */

public class SpectrumKeyFunctionImpl implements SpectrumKeyFunction<MsnSpectrum> {

    private Pattern regEx = Pattern.compile("\\.[0]+");

    @Override
    public String apply(SpectrumIdentifier identifier) {

        String comment = identifier.getSpectrum();

        int firstSpace = comment.indexOf(' ');

        if (firstSpace>=1) comment = comment.substring(0,firstSpace);

        comment = regEx.matcher(comment).replaceAll(".");

        return comment;
    }

}
