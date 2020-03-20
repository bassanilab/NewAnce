package newance.mzjava.ms.io.mgf;

import newance.mzjava.ms.peaklist.PeakList;

import java.text.NumberFormat;

import static com.google.common.base.Preconditions.checkNotNull;

/**
 * @author fnikitin
 */
public class PeakListPrecisionFormat {

    private PeakListPrecisionFormat() {

    }

    public static NumberFormat getMzFormat(PeakList.Precision precision) {

        switch (precision) {

            case DOUBLE:

                return NumberFormatFactory.DOUBLE_PRECISION;
            case FLOAT:

                return NumberFormatFactory.FLOAT_PRECISION;
            case DOUBLE_FLOAT:

                return NumberFormatFactory.DOUBLE_PRECISION;

            case DOUBLE_CONSTANT:

                return NumberFormatFactory.DOUBLE_PRECISION;

            case FLOAT_CONSTANT:

                return NumberFormatFactory.FLOAT_PRECISION;
            default:

                throw new IllegalArgumentException(precision + ": not a valid peak list precision!");
        }
    }

    public static NumberFormat getIntensityFormat(PeakList.Precision precision) {

        checkNotNull(precision);

        switch (precision) {

            case DOUBLE:

                return NumberFormatFactory.DOUBLE_PRECISION;
            case FLOAT:

                return NumberFormatFactory.FLOAT_PRECISION;
            case DOUBLE_FLOAT:

                return NumberFormatFactory.FLOAT_PRECISION;
            case DOUBLE_CONSTANT:

                return NumberFormatFactory.INT_PRECISION;
            case FLOAT_CONSTANT:

                return NumberFormatFactory.INT_PRECISION;
            default:

                throw new IllegalArgumentException(precision + ": not a valid peak list precision!");
        }
    }

}
