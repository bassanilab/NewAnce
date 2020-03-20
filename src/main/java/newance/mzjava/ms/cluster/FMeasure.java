package newance.mzjava.ms.cluster;

import com.google.common.base.Preconditions;

import java.util.Collection;
import java.util.Set;

/**
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class FMeasure {

    private FMeasure() {}

    private static double calcClusterFMeasure(int tp, int expectedSize, int actualSize) {

        if (expectedSize == 0 || actualSize == 0) {

            return 0;
        } else {

            double precision = (double)tp / actualSize;
            double recall = (double)tp / expectedSize;
            return 2 * precision * recall / (precision + recall);
        }
    }

    protected static <E> int intersect(Set<E> expectedCluster, Set<E> actualCluster) {

        Set<E> smallerCluster;
        Set<E> largerCluster;

        if (expectedCluster.size() < actualCluster.size()) {

            smallerCluster = expectedCluster;
            largerCluster = actualCluster;
        } else {

            smallerCluster = actualCluster;
            largerCluster = expectedCluster;
        }

        int overlap = 0;
        for (E event : smallerCluster) {

            if (largerCluster.contains(event)) {

                overlap = overlap + 1;
            }
        }

        return overlap;
    }

    protected static int calcTotalMembers(Collection<? extends Set> clusters) {

        int totalEvent = 0;
        for (Set cluster : clusters) {

            totalEvent += cluster.size();
        }

        return totalEvent;
    }

    public static <E> double calcFMeasure(Collection<Set<E>> expectedClusters, Collection<Set<E>> actualClusters) {

        Preconditions.checkArgument(!actualClusters.isEmpty());

        double fMeasure = 0;
        double totalMembers = calcTotalMembers(expectedClusters);

        for (Set<E> expectedCluster : expectedClusters) {

            double maxClusterF = 0.0;
            final int expectedSize = expectedCluster.size();

            for (Set<E> actualCluster : actualClusters) {

                int tp = intersect(expectedCluster, actualCluster);
                double clusterF = tp > 0 ? calcClusterFMeasure(tp, expectedSize, actualCluster.size()) : 0;
                maxClusterF = Math.max(maxClusterF, clusterF);
            }

            fMeasure += maxClusterF * expectedSize / totalMembers;
        }

        return fMeasure;
    }
}
