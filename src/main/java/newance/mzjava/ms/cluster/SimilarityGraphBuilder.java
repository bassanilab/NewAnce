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

/**
 * Common interface for similarity graph builders.
 *
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public interface SimilarityGraphBuilder<V, G extends SimilarityGraph<V>> {

    /**
     * Add the <code>vertex</code> to the graph
     *
     * @param vertex the vertex to add
     * @return this builder
     */
    SimilarityGraphBuilder<V, G> add(V vertex);

    /**
     * Create and add an edge from <code>vertex1</code> to <code>vertex2</code> with similarity
     * equal to <code>score</code> to the graph.
     * <p/>
     * If a vertex is not present in the graph the vertex is added to the graph.
     *
     * @param vertex1 the parent vertex
     * @param vertex2 the child vertex
     * @param score   the similarity of the edge
     * @return the created edge
     */
    SimEdge<V> add(V vertex1, V vertex2, double score);

    /**
     * Return the number of edges currently in the graph
     *
     * @return the number of edges currently in the graph
     */
    int edgeCount();

    /**
     * Add the <code>simEdge</code> to the graph.
     *
     * @param simEdge the edge to add
     * @return this builder
     */
    SimilarityGraphBuilder<V, G> add(SimEdge<V> simEdge);

    G build();

    /**
     * Returns true if build can be called more than once, false otherwise.
     *
     * @return true if build can be called more than once, false otherwise
     */
    boolean isReusable();
}
