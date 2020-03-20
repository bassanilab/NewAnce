package newance.mzjava.ms.cluster;

import com.google.common.base.Optional;
import com.google.common.collect.Sets;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.custom_hash.TObjectIntCustomHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.procedure.TObjectProcedure;
import gnu.trove.strategy.IdentityHashingStrategy;
import newance.mzjava.ms.library.Procedure;

import java.util.*;

/**
 * A similarity graph that is backed by a matrix
 *
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class DenseSimilarityGraph<V> implements SimilarityGraph<V> {

    private final TObjectIntMap<V> vertexIndexMap;
    private final SimEdge<V>[][] adjacencyList;
    private final int edgeCount;
    private static final String MESSAGE_DOES_NOT_EXIST = " does not exist in this graph";

    protected DenseSimilarityGraph(TObjectIntMap<V> vertexIndexMap, SimEdge<V>[][] adjacencyList, int edgeCount) {

        this.vertexIndexMap = vertexIndexMap;
        this.adjacencyList = adjacencyList;
        this.edgeCount = edgeCount;
    }

    @Override
    public int getVertexCount() {

        return adjacencyList.length;
    }

    @Override
    public Iterable<V> getVertices() {

        return Collections.unmodifiableSet(vertexIndexMap.keySet());
    }

    @Override
    public void forEachVertex(final Procedure<V> procedure) {

        vertexIndexMap.forEachKey(new TObjectProcedure<V>() {
            @Override
            public boolean execute(V object) {

                procedure.execute(object);
                return true;
            }
        });
    }

    @Override
    public Iterable<SimEdge<V>> getEdges() {

        int size = adjacencyList.length;
        final ArrayList<SimEdge<V>> edgeList = new ArrayList<>(size * (size - 1));

        forEachEdge(new Procedure<SimEdge<V>>() {
            @Override
            public void execute(SimEdge<V> edge) {

                edgeList.add(edge);
            }
        });

        edgeList.trimToSize();
        return edgeList;
    }

    @Override
    public void forEachEdge(Procedure<SimEdge<V>> procedure) {

        for (int i = 0; i < adjacencyList.length; i++) {

            SimEdge<V>[] simEdges = adjacencyList[i];
            for (int j = i; j < simEdges.length; j++) {

                SimEdge<V> simEdge = simEdges[j];
                if (simEdge != null) {

                    procedure.execute(simEdge);
                }
            }
        }
    }

    @Override
    public Iterable<V> getNeighbors(V vertex) {

        final List<V> neighbours = new ArrayList<>(adjacencyList.length);

        forEachNeighbour(vertex, new Procedure<V>() {
            @Override
            public void execute(V v) {

                neighbours.add(v);
            }
        });

        return neighbours;
    }

    @Override
    public void forEachNeighbour(V vertex, Procedure<V> procedure) {

        int index = vertexIndexMap.get(vertex);
        if (index == -1)
            throw new IllegalStateException(vertex + MESSAGE_DOES_NOT_EXIST);

        SimEdge<V>[] simEdges = adjacencyList[index];
        for (SimEdge<V> simEdge : simEdges) {

            if (simEdge != null)
                procedure.execute(simEdge.getOther(vertex));
        }
    }

    @Override
    public Iterable<SimEdge<V>> getEdges(V vertex) {

        final List<SimEdge<V>> outEdge = new ArrayList<>(adjacencyList.length);

        forEachEdge(vertex, new Procedure<SimEdge<V>>() {
            @Override
            public void execute(SimEdge<V> simEdge) {

                outEdge.add(simEdge);
            }
        });

        return outEdge;
    }

    @Override
    public void forEachEdge(V vertex, Procedure<SimEdge<V>> procedure) {

        int index = vertexIndexMap.get(vertex);
        if (index == -1)
            throw new IllegalStateException(vertex + MESSAGE_DOES_NOT_EXIST);

        SimEdge<V>[] simEdges = adjacencyList[index];
        for (SimEdge<V> simEdge : simEdges) {

            if (simEdge != null)
                procedure.execute(simEdge);
        }
    }

    @Override
    public Optional<SimEdge<V>> findEdge(V vertex1, V vertex2) {

        int index1 = vertexIndexMap.get(vertex1);
        if(index1 == -1)
            throw new IllegalStateException(vertex1 + MESSAGE_DOES_NOT_EXIST);
        int index2 = vertexIndexMap.get(vertex2);
        if(index2 == -1)
            throw new IllegalStateException(vertex1 + MESSAGE_DOES_NOT_EXIST);

        return Optional.fromNullable(adjacencyList[index1][index2]);
    }

    @Override
    public int degree(V vertex) {

        final int[] count = new int[]{0};
        forEachNeighbour(vertex, new Procedure<V>() {
            @Override
            public void execute(V v) {
                count[0]++;
            }
        });

        return count[0];
    }

    @Override
    public boolean containsVertex(V vertex) {

        return vertexIndexMap.containsKey(vertex);
    }

    @Override
    public int getEdgeCount() {

        return edgeCount;
    }

    public static <V> Builder<V> builder() {

        return new Builder<>(false);
    }

    public static class Builder<V> extends AbstractBuilder<V, DenseSimilarityGraph<V>> {

        public Builder() {

            super(false);
        }

        public Builder(boolean useIdentityHash) {

            super(useIdentityHash);
        }

        protected DenseSimilarityGraph<V> doBuild(TObjectIntMap<V> vertexIndexMap, int edgeCount, SimEdge<V>[][] adjacencyList) {

            return new DenseSimilarityGraph<>(vertexIndexMap, adjacencyList, edgeCount);
        }
    }

    protected static abstract class AbstractBuilder<V, G extends SimilarityGraph<V>> implements SimilarityGraphBuilder<V, G> {

        private final Set<V> vertices;
        private final List<SimEdge<V>> edgeList = new ArrayList<>();
        private final boolean useIdentityHash;

        public AbstractBuilder(boolean useIdentityHash) {

            this.useIdentityHash = useIdentityHash;

            if (useIdentityHash) {

                vertices = Sets.newIdentityHashSet();
            } else {

                vertices = new HashSet<>();
            }
        }

        @Override
        public SimilarityGraphBuilder<V, G> add(V vertex) {

            vertices.add(vertex);

            return this;
        }

        @Override
        public SimEdge<V> add(V vertex1, V vertex2, double score) {

            add(vertex1);
            add(vertex2);

            SimEdge<V> edge = new SimEdge<>(vertex1, vertex2, score);

            edgeList.add(edge);

            return edge;
        }

        public G build() {

            TObjectIntMap<V> vertexIndexMap;
            if (useIdentityHash) {
                vertexIndexMap = new TObjectIntCustomHashMap<>(new IdentityHashingStrategy<>(), vertices.size(), 0.1f, -1);
            } else {
                vertexIndexMap = new TObjectIntHashMap<>(vertices.size(), 0.1f, -1);
            }
            int edgeCount = edgeList.size();
            int id = 0;
            for (V vertex : vertices) {

                vertexIndexMap.put(vertex, id++);
            }

            int size = vertices.size();
            //noinspection unchecked
            SimEdge<V>[][] adjacencyList = new SimEdge[size][size];
            for (SimEdge<V> edge : edgeList) {

                int index1 = vertexIndexMap.get(edge.getVertex1());
                int index2 = vertexIndexMap.get(edge.getVertex2());

                adjacencyList[index1][index2] = edge;
                adjacencyList[index2][index1] = edge;
            }

            vertices.clear();
            edgeList.clear();

            return doBuild(vertexIndexMap, edgeCount, adjacencyList);
        }

        protected abstract G doBuild(TObjectIntMap<V> vertexIndexMap, int edgeCount, SimEdge<V>[][] adjacencyList);

        @Override
        public boolean isReusable() {

            return true;
        }

        public int edgeCount() {

            return edgeList.size();
        }

        @Override
        public SimilarityGraphBuilder<V, G> add(SimEdge<V> simEdge) {

            add(simEdge.getVertex1());
            add(simEdge.getVertex2());
            edgeList.add(simEdge);

            return this;
        }
    }
}