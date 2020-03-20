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

import newance.mzjava.ms.library.Procedure;
import org.apache.commons.collections15.Buffer;
import org.apache.commons.collections15.buffer.UnboundedFifoBuffer;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

/**
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class GlobalThresholdClusterBuilder<V> implements ClusterBuilder<V> {

    private final double scoreThreshold;

    public GlobalThresholdClusterBuilder(double scoreThreshold) {

        this.scoreThreshold = scoreThreshold;
    }

    @Override
    public Collection<Set<V>> cluster(SimilarityGraph<V> graph) {

        Set<Set<V>> clusters = new HashSet<>();

        final Set<V> unvisitedVertices = new HashSet<>(graph.getVertexCount());
        graph.forEachVertex(new Procedure<V>() {
            @Override
            public void execute(V v) {
                unvisitedVertices.add(v);
            }
        });

        while (!unvisitedVertices.isEmpty()) {
            final Set<V> cluster = new HashSet<>();
            V root = unvisitedVertices.iterator().next();
            unvisitedVertices.remove(root);
            cluster.add(root);

            final Buffer<V> queue = new UnboundedFifoBuffer<>();
            queue.add(root);

            while (!queue.isEmpty()) {

                final V currentVertex = queue.remove();
                graph.forEachEdge(currentVertex, new Procedure<SimEdge<V>>() {

                    @Override
                    public void execute(SimEdge<V> simEdge) {

                        if(simEdge.getScore() < scoreThreshold)
                            return;

                        V neighbor = simEdge.getOther(currentVertex);
                        if (unvisitedVertices.contains(neighbor)) {
                            queue.add(neighbor);
                            unvisitedVertices.remove(neighbor);
                            cluster.add(neighbor);
                        }
                    }
                });
            }
            clusters.add(cluster);
        }

        return clusters;
    }

    @Override
    public Collection<Set<V>> cluster(SimilarityGraph<V> graph, Collection<Set<V>> startingClusters) {

        return cluster(graph);
    }
}
