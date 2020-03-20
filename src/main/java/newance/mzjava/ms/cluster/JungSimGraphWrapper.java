package newance.mzjava.ms.cluster;

import com.google.common.base.Optional;
import com.google.common.collect.Lists;
import edu.uci.ics.jung.graph.AbstractTypedGraph;
import edu.uci.ics.jung.graph.UndirectedGraph;
import edu.uci.ics.jung.graph.util.EdgeType;
import edu.uci.ics.jung.graph.util.Pair;

import java.util.Collection;

/**
 * Wraps a SimilarityGraph in a jung.UndirectedGraph so that the JUNG algorithms can be used.
 *
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class JungSimGraphWrapper<V> extends AbstractTypedGraph<V, SimEdge<V>> implements UndirectedGraph<V, SimEdge<V>> {

    private static final String MESSAGE_GRAPH_IS_IMMUTABLE = "SimilarityGraphs are immutable";

    private final SimilarityGraph<V> simGraph;

    public JungSimGraphWrapper(SimilarityGraph<V> simGraph) {

        super(EdgeType.UNDIRECTED);
        this.simGraph = simGraph;
    }

    @Override
    public boolean addEdge(SimEdge<V> edge, Pair<? extends V> endpoints, EdgeType edgeType) {

        throw new UnsupportedOperationException(MESSAGE_GRAPH_IS_IMMUTABLE);
    }

    @Override
    public Collection<SimEdge<V>> getInEdges(V vertex) {

        return this.getIncidentEdges(vertex);
    }

    @Override
    public Collection<SimEdge<V>> getOutEdges(V vertex) {

        return this.getIncidentEdges(vertex);
    }

    @Override
    public Collection<V> getPredecessors(V vertex) {

        return this.getNeighbors(vertex);
    }

    @Override
    public Collection<V> getSuccessors(V vertex) {

        return this.getNeighbors(vertex);
    }

    @Override
    public SimEdge<V> findEdge(V v1, V v2) {

        if (!containsVertex(v1) || !containsVertex(v2)){
            //The jung API wants null here, see UndirectedSparseGraph
            return null;
        }

        return simGraph.findEdge(v1, v2).orNull();
    }

    @Override
    public Pair<V> getEndpoints(SimEdge<V> edge) {

        return new Pair<V>(edge.getVertex1(), edge.getVertex2());
    }

    @Override
    public V getSource(SimEdge<V> directedEdge) {

        //The jung API wants null here, see UndirectedSparseGraph
        return null;
    }

    @Override
    public V getDest(SimEdge<V> directedEdge) {

        //The jung API wants null here, see UndirectedSparseGraph
        return null;
    }

    @Override
    public boolean isSource(V vertex, SimEdge<V> edge) {

        return false;
    }

    @Override
    public boolean isDest(V vertex, SimEdge<V> edge) {

        return false;
    }

    @Override
    public Collection<SimEdge<V>> getEdges() {

        return newArrayList(simGraph.getEdges());
    }

    @Override
    public Collection<V> getVertices() {

        return newArrayList(simGraph.getVertices());
    }

    @Override
    public boolean containsVertex(V vertex) {

        return simGraph.containsVertex(vertex);
    }

    @Override
    public boolean containsEdge(SimEdge<V> edge) {

        Optional<SimEdge<V>> edgeOpt = simGraph.findEdge(edge.getVertex1(), edge.getVertex2());
        return edgeOpt.isPresent() && edgeOpt.get().equals(edge);
    }

    @Override
    public int getEdgeCount() {

        return simGraph.getEdgeCount();
    }

    @Override
    public int getVertexCount() {

        return simGraph.getVertexCount();
    }

    @Override
    public Collection<V> getNeighbors(V vertex) {

        if (!containsVertex(vertex)) {
            //The jung API wants null here, see UndirectedSparseGraph
            return null;
        }
        return newArrayList(simGraph.getNeighbors(vertex));
    }

    @Override
    public Collection<SimEdge<V>> getIncidentEdges(V vertex) {

        if (!containsVertex(vertex)) {
            //The jung API wants null here, see UndirectedSparseGraph
            return null;
        }
        return newArrayList(simGraph.getEdges(vertex));
    }

    private <N> Collection<N> newArrayList(Iterable<N> iterable) {

        if(iterable instanceof Collection)
            return (Collection<N>) iterable;
        else
            return Lists.newArrayList(iterable);
    }

    @Override
    public boolean addVertex(V vertex) {

        throw new UnsupportedOperationException(MESSAGE_GRAPH_IS_IMMUTABLE);
    }

    @Override
    public boolean removeVertex(V vertex) {

        throw new UnsupportedOperationException(MESSAGE_GRAPH_IS_IMMUTABLE);
    }

    @Override
    public boolean removeEdge(SimEdge<V> edge) {

        throw new UnsupportedOperationException(MESSAGE_GRAPH_IS_IMMUTABLE);
    }
}
