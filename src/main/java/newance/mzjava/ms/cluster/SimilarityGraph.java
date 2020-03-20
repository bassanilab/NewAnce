/**
 * Copyright (c) 2010, SIB. All rights reserved.
 *
 * SIB (Swiss Institute of Bioinformatics) - http://www.isb-sib.ch Host -
 * https://sourceforge.net/projects/javaprotlib/
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer. Redistributions in binary
 * form must reproduce the above copyright notice, this list of conditions and
 * the following disclaimer in the documentation and/or other materials provided
 * with the distribution. Neither the name of the SIB/GENEBIO nor the names of
 * its contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL SIB/GENEBIO BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package newance.mzjava.ms.cluster;

import com.google.common.base.Optional;
import newance.mzjava.ms.library.Procedure;

/**
 * Interface for a graph that represents the similarity between vertexes.
 *
 * @param <V> the type of vertexes in this graph
 *
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public interface SimilarityGraph<V> {

    /**
     * Return the number of vertexes in this graph.
     *
     * @return the number of vertexes in this graph
     */
    int getVertexCount();

    /**
     * Return an Iterable that iterates over the vertices in this graph.
     *
     * @return an Iterable that iterates over the vertices in this graph
     */
    Iterable<V> getVertices();

    /**
     * Applies the procedure to each vertex in this graph
     *
     * @param procedure a <code>Procedure</code> value
     */
    void forEachVertex(Procedure<V> procedure);

    /**
     * Return an Iterable that iterates over the SimEdge's in this graph.
     *
     * @return an Iterable that iterates over the SimEdge's in this graph
     */
    Iterable<SimEdge<V>> getEdges();

    /**
     * Applies the procedure to each edge in this graph
     *
     * @param procedure a <code>Procedure</code> value
     */
    void forEachEdge(Procedure<SimEdge<V>> procedure);

    /**
     * Return an Iterable that iterates over the neighbours of <code>vertex</code>.
     *
     * @param vertex the vertex whose neighbours are to be iterated
     * @return an Iterable that iterates over the neighbours of <code>vertex</code>
     */
    Iterable<V> getNeighbors(V vertex);

    /**
     * Applies the procedure to each vertex that is a neighbour of <code>vertex</code>
     *
     * @param vertex the vertex whose neighbours are to have the procedure applied
     * @param procedure a <code>Procedure</code> value
     */
    void forEachNeighbour(V vertex, Procedure<V> procedure);

    /**
     * Return an Iterable that iterates over all the edges connected to <code>vertex</code>.
     *
     * @param vertex the vertex whose edges are to be iterated
     * @return an Iterable that iterates over all the edges connected to <code>vertex</code>
     */
    Iterable<SimEdge<V>> getEdges(V vertex);

    /**
     * Applies the procedure to each edge that is connected to <code>vertex</code>
     *
     * @param vertex the vertex whose edges are to have the procedure applied
     * @param procedure a <code>Procedure</code> value
     */
    void forEachEdge(V vertex, Procedure<SimEdge<V>> procedure);

    /**
     * Find the edge between <code>vertex1</code> and <code>vertex2</code>. If there is no edge between
     * <code>vertex1</code> and <code>vertex2</code> an Optional that contains no SimEdge is returned.
     *
     * @param vertex1 the start vertex
     * @param vertex2 the end vertex
     * @return an Optional that constrains the SimEdge between <code>vertex1</code> and <code>vertex2</code>. If no
     * such edge exists an Optional that contains no reference is returned
     */
    Optional<SimEdge<V>> findEdge(V vertex1, V vertex2);

    /**
     * Return the number of edges incident to <code>vertex</code>.
     * @param vertex the vertex to check
     * @return the number of edges incident to <code>vertex</code>
     */
    int degree(V vertex);

    /**
     * Returns true if this graph contains <code>vertex</code>, false otherwise.
     *
     * @param vertex the vertex to check
     * @return true if this graph contains <code>vertex</code>, false otherwise
     */
    boolean containsVertex(V vertex);

    /**
     * Returns the number of edges in this graph
     *
     * @return the number of edges in this graph
     */
    int getEdgeCount();
}
