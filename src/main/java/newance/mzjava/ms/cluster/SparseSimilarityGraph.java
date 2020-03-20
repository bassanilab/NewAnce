package newance.mzjava.ms.cluster;

import com.google.common.base.Optional;
import com.google.common.base.Preconditions;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import newance.mzjava.ms.library.Procedure;

import java.util.*;

/**
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class SparseSimilarityGraph<V> implements SimilarityGraph<V> {

    private final List<SimEdge<V>> edges;
    private final Map<V, TObjectIntMap<V>> adjacencyMap;
    private static final String MESSAGE_VERTEX_IS_NOT_A_NODE_IN_THIS_GRAPH = "vertex is not a node in this graph";

    private SparseSimilarityGraph(List<SimEdge<V>> edges, Map<V, TObjectIntMap<V>> adjacencyMap) {

        this.edges = edges;
        this.adjacencyMap = adjacencyMap;
    }

    @Override
    public int getVertexCount() {

        return adjacencyMap.size();
    }

    @Override
    public Collection<V> getVertices() {

        return Collections.unmodifiableSet(adjacencyMap.keySet());
    }

    @Override
    public void forEachVertex(Procedure<V> procedure) {

        for(V vertex : adjacencyMap.keySet()) {

            procedure.execute(vertex);
        }
    }

    @Override
    public Iterable<SimEdge<V>> getEdges() {

        return Collections.unmodifiableList(edges);
    }

    @Override
    public void forEachEdge(Procedure<SimEdge<V>> procedure) {

        for(SimEdge<V> edge : edges) {

            procedure.execute(edge);
        }
    }

    @Override
    public Collection<V> getNeighbors(V vertex) {

        TObjectIntMap<V> neighbourMap = adjacencyMap.get(vertex);
        if (neighbourMap == null)
            throw new IllegalStateException(MESSAGE_VERTEX_IS_NOT_A_NODE_IN_THIS_GRAPH);

        return Collections.unmodifiableSet(neighbourMap.keySet());
    }

    @Override
    public void forEachNeighbour(V vertex, Procedure<V> procedure) {

        TObjectIntMap<V> neighbourMap = adjacencyMap.get(vertex);
        if (neighbourMap == null)
            throw new IllegalStateException(MESSAGE_VERTEX_IS_NOT_A_NODE_IN_THIS_GRAPH);

        for(V current : neighbourMap.keySet()) {

            procedure.execute(current);
        }
    }

    @Override
    public Iterable<SimEdge<V>> getEdges(V vertex) {

        TObjectIntMap<V> neighbourMap = adjacencyMap.get(vertex);
        if (neighbourMap == null)
            throw new IllegalStateException(MESSAGE_VERTEX_IS_NOT_A_NODE_IN_THIS_GRAPH);

        Set<SimEdge<V>> outEdges = new HashSet<>(neighbourMap.size());
        for(int i : neighbourMap.values()) {

            outEdges.add(edges.get(i));
        }

        return Collections.unmodifiableSet(outEdges);
    }

    @Override
    public void forEachEdge(V vertex, Procedure<SimEdge<V>> procedure) {

        TObjectIntMap<V> neighbourMap = adjacencyMap.get(vertex);
        if (neighbourMap == null)
            throw new IllegalStateException(MESSAGE_VERTEX_IS_NOT_A_NODE_IN_THIS_GRAPH);

        for(int i : neighbourMap.values()) {

            procedure.execute(edges.get(i));
        }
    }

    @Override
    public Optional<SimEdge<V>> findEdge(V vertex1, V vertex2) {

        TObjectIntMap<V> neighbourMap = adjacencyMap.get(vertex1);
        if (neighbourMap == null)
            throw new IllegalStateException("vertex1 is not anode in this graph");
        if(!adjacencyMap.containsKey(vertex2))
            throw new IllegalStateException("vertex2 is not anode in this graph");

        int edgeIndex = neighbourMap.get(vertex2);

        if(edgeIndex == neighbourMap.getNoEntryValue())
            return Optional.absent();

        return Optional.of(edges.get(edgeIndex));
    }

    @Override
    public int degree(V vertex) {

        TObjectIntMap<V> neighbourMap = adjacencyMap.get(vertex);
        if (neighbourMap == null)
            throw new IllegalStateException(MESSAGE_VERTEX_IS_NOT_A_NODE_IN_THIS_GRAPH);

        return neighbourMap.size();
    }

    @Override
    public boolean containsVertex(V vertex) {

        return adjacencyMap.containsKey(vertex);
    }

    @Override
    public int getEdgeCount() {

        return edges.size();
    }

    public static <V> Builder<V> builder(int edgeCapacity) {

        return new Builder<>(edgeCapacity);
    }

    public static class Builder<V> implements SimilarityGraphBuilder<V, SparseSimilarityGraph<V>> {

        private final List<SimEdge<V>> edges = new ArrayList<>();
        private final Map<V, TObjectIntMap<V>> adjacencyMap = new HashMap<>();
        private final int edgeCapacity;

        private boolean built = false;

        public Builder(int edgeCapacity) {

            this.edgeCapacity = edgeCapacity;
        }

        public Builder<V> add(V vertex) {

            Preconditions.checkNotNull(vertex);
            if (!adjacencyMap.containsKey(vertex)) {

                adjacencyMap.put(vertex, new TObjectIntHashMap<V>(edgeCapacity, 0.5f, -1));
            }

            return this;
        }

        public SimEdge<V> add(V vertex1, V vertex2, double score){

            SimEdge<V> edge = new SimEdge<>(vertex1, vertex2, score);
            add(edge);
            return edge;
        }

        public Builder<V> add(SimEdge<V> edge){

            Preconditions.checkArgument(edge.getScore() >= 0 && edge.getScore() <= 1);

            add(edge.getVertex1());
            add(edge.getVertex2());
            edges.add(edge);
            return this;
        }

        public SparseSimilarityGraph<V> build(){

            if(built)
                throw new IllegalStateException("This builder cannot be reused");

            built = true;

            for (int i = 0; i < edges.size(); i++) {

                SimEdge<V> edge = edges.get(i);
                this.adjacencyMap.get(edge.getVertex1()).put(edge.getVertex2(), i);
                this.adjacencyMap.get(edge.getVertex2()).put(edge.getVertex1(), i);
            }
            return new SparseSimilarityGraph<>(edges, adjacencyMap);
        }

        @Override
        public boolean isReusable() {

            return false;
        }

        public int edgeCount() {

            return edges.size();
        }
    }
}
