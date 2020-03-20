package newance.mzjava.ms.cluster;

import com.google.common.base.Optional;
import com.google.common.base.Preconditions;
import com.google.common.collect.Sets;
import newance.mzjava.ms.library.Procedure;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomAdaptor;
import org.apache.commons.math3.stat.clustering.Cluster;
import org.apache.commons.math3.stat.clustering.Clusterable;
import org.apache.commons.math3.stat.clustering.KMeansPlusPlusClusterer;

import java.util.*;

/**
 * KMeans cluster builder that uses the algorithm provided by commons math KMeansPlusPlusClusterer.
 *
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class KMeansPlusPlusClusterBuilder<V> implements ClusterBuilder<V> {

    private final int numClusters;
    private final Optional<Long> seed;

    public KMeansPlusPlusClusterBuilder(int numClusters) {

        this.numClusters = numClusters;
        this.seed = Optional.absent();
    }

    public KMeansPlusPlusClusterBuilder(int numClusters, Optional<Long> seed) {

        Preconditions.checkNotNull(seed);

        this.numClusters = numClusters;
        this.seed = seed;
    }

    @Override
    public Collection<Set<V>> cluster(final SimilarityGraph<V> graph) {

        if (graph.getVertexCount() == 1)
            return Collections.singletonList((Set<V>) Sets.newHashSet(graph.getVertices()));

        final List<SimPoint<V>> points = new ArrayList<>(graph.getVertexCount());
        graph.forEachVertex(new Procedure<V>() {
            @Override
            public void execute(V v) {

                points.add(new SimPoint<>(v, graph));
            }
        });

        KMeansPlusPlusClusterer<SimPoint<V>> kMeans = new KMeansPlusPlusClusterer<>(new RandomAdaptor(
                new MersenneTwister(seed.or(System.currentTimeMillis() + System.identityHashCode(this))))
        );
        List<Cluster<SimPoint<V>>> clusters = kMeans.cluster(points, numClusters, 1000);

        List<Set<V>> results = new ArrayList<>(clusters.size());
        for (Cluster<SimPoint<V>> cluster : clusters) {

            final List<SimPoint<V>> pointsInCluster = cluster.getPoints();
            final Set<V> set = new HashSet<>(pointsInCluster.size());
            for (SimPoint<V> point : pointsInCluster) {

                set.add(point.vertex);
            }

            if (!set.isEmpty()) {

                results.add(set);
            }
        }
        return results;
    }

    @Override
    public Collection<Set<V>> cluster(SimilarityGraph<V> graph, Collection<Set<V>> startingClusters) {

        throw new UnsupportedOperationException("SimpleKMeansClusterBuilder cant use starting clusters");
    }

    private static class SimPoint<V> implements Clusterable<SimPoint<V>> {

        private final V vertex;
        private final SimilarityGraph<V> similarityGraph;

        private SimPoint(V vertex, SimilarityGraph<V> similarityGraph) {

            this.vertex = vertex;
            this.similarityGraph = similarityGraph;
        }

        @Override
        public double distanceFrom(SimPoint<V> p) {

            if (vertex.equals(p.vertex)) {

                return 0;
            } else {

                final Optional<SimEdge<V>> edgeOptional = similarityGraph.findEdge(vertex, p.vertex);
                return edgeOptional.isPresent() ? 1 - edgeOptional.get().getScore() : 1;
            }
        }

        @Override
        public SimPoint<V> centroidOf(Collection<SimPoint<V>> p) {

            if (p.size() == 1) return p.iterator().next();

            SimPoint<V> centroid = null;
            double maxScore = 0;
            for (SimPoint<V> candidate : p) {

                double score = 0;

                for (SimPoint<V> neighbour : p) {

                    if (candidate != neighbour) {
                        final Optional<SimEdge<V>> edgeOptional = similarityGraph.findEdge(candidate.vertex, neighbour.vertex);
                        if (edgeOptional.isPresent()) {

                            score += edgeOptional.get().getScore();
                        }
                    }
                }

                if (centroid == null || score > maxScore) {

                    maxScore = score;
                    centroid = candidate;
                }
            }

            if (centroid == null)
                throw new IllegalStateException("Could not find a centroid");

            return centroid;
        }
    }
}