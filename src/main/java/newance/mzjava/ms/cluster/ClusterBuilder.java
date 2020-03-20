/**
 * Copyright (c) 2010, SIB. All rights reserved.
 *
 * SIB (Swiss Institute of Bioinformatics) - http://www.isb-sib.ch Host -
 * http://mzjava.expasy.org
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

import java.util.Collection;
import java.util.Set;

/**
 * Interface for any object that can cluster the vertexes of a SimilarityGraph
 *
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public interface ClusterBuilder<V> {

    /**
     * Cluster the vertices given a SimilarityGraph
     *
     * @param graph the graph
     * @return collection containing sets of vertices that belong to the same cluster
     */
    Collection<Set<V>> cluster(SimilarityGraph<V> graph);

    /**
     * Cluster the vertices in <code>graph</code> using the <code>startingClusters</code> as a starting point.
     * <p>
     * Cluster  builders that where the starting clusters have no impact on clustering, such as the
     * threshold cluster builder, can delegate to the cluster(SimilarityGraph<V> graph) method.
     *
     * @param graph the graph
     * @param startingClusters the clusters to start with
     * @return collection containing sets of vertices that belong to the same cluster
     */
    Collection<Set<V>> cluster(SimilarityGraph<V> graph, Collection<Set<V>> startingClusters);
}
