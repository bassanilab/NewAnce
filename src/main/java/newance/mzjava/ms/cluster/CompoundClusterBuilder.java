package newance.mzjava.ms.cluster;

import java.util.Collection;
import java.util.Set;

/**
 * Cluster builder that uses the output of builder1 as the starting clusters of builder2.
 *
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class CompoundClusterBuilder<V> implements ClusterBuilder<V> {

    private final ClusterBuilder<V> builder1;
    private final ClusterBuilder<V> builder2;

    /**
     * Construct a new CompoundClusterBuilder
     *
     * @param builder1 the first builder
     * @param builder2 the second builder
     */
    public CompoundClusterBuilder(ClusterBuilder<V> builder1, ClusterBuilder<V> builder2) {

        this.builder1 = builder1;
        this.builder2 = builder2;
    }

    @Override
    public Collection<Set<V>> cluster(SimilarityGraph<V> graph) {

        return builder2.cluster(graph, builder1.cluster(graph));
    }

    @Override
    public Collection<Set<V>> cluster(SimilarityGraph<V> graph, Collection<Set<V>> startingClusters) {

        return builder2.cluster(graph, builder1.cluster(graph, startingClusters));
    }
}
