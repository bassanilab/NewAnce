package newance.mzjava.ms.peaklist;

/**
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class PeakListUtils {

    private PeakListUtils() {}

    public static double calcLength(double[] intensityList, int size) {

        double dist = 0;

        for (int i = 0; i < size; i++) {

            double x = intensityList[i];

            dist += x * x;
        }

        return Math.sqrt(dist);
    }

    public static double calcLength(float[] intensityList, int size) {

        double dist = 0;

        for (int i = 0; i < size; i++) {

            double x = intensityList[i];

            dist += x * x;
        }

        return Math.sqrt(dist);
    }
}
