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

import edu.uci.ics.jung.algorithms.cluster.WeakComponentClusterer;
import edu.uci.ics.jung.algorithms.shortestpath.MinimumSpanningForest;
import edu.uci.ics.jung.graph.DelegateForest;
import edu.uci.ics.jung.graph.Forest;
import newance.mzjava.ms.library.Procedure;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * @author Markus Muller
 * @version sqrt -1
 */
public class MSTClusterBuilder<V> implements ClusterBuilder<V> {

    private final double scoreThreshold;
    private final DoubleFunction<SimEdge<V>> scoreMapper;

    public MSTClusterBuilder(double scoreThreshold) {

        this.scoreThreshold = scoreThreshold;
        this.scoreMapper = new DoubleFunction<SimEdge<V>>() {
            @Override
            public double apply(SimEdge<V> simEdge) {
                return 1.0 - simEdge.getScore();
            }
        };
    }

    public MSTClusterBuilder(double scoreThreshold,DoubleFunction<SimEdge<V>> scoreMapper) {

        this.scoreThreshold = scoreThreshold;
        this.scoreMapper = scoreMapper;
    }

    @Override
    public Collection<Set<V>> cluster(SimilarityGraph<V> graph) {

        final Map<SimEdge<V>, Double> scoreMap = new HashMap<SimEdge<V>, Double>();
        graph.forEachEdge(new Procedure<SimEdge<V>>() {
            @Override
            public void execute(SimEdge<V> simEdge) {
                scoreMap.put(simEdge, scoreMapper.apply(simEdge));
            }
        });

        JungSimGraphWrapper<V> jGraphWrapper = new JungSimGraphWrapper<V>(graph);

        MinimumSpanningForest<V,SimEdge<V>> mstforest =
                new MinimumSpanningForest<V,SimEdge<V>>(jGraphWrapper,new DelegateForest<V, SimEdge<V>>(),null,scoreMap);

        Forest<V,SimEdge<V>> forest = mstforest.getForest();

        Forest<V,SimEdge<V>> cutForest = cutForest(forest);

        return getClusters(cutForest);
    }

    @Override
    public Collection<Set<V>> cluster(SimilarityGraph<V> graph, Collection<Set<V>> startingClusters) {
        return null;
    }

    private Forest<V,SimEdge<V>> cutForest(Forest<V,SimEdge<V>> forest) {

        Forest<V, SimEdge<V>> cutForest = new DelegateForest<V, SimEdge<V>>();
        for (V node : forest.getVertices()) {

            cutForest.addVertex(node);
        }
        for (SimEdge<V> edge : forest.getEdges()) {

            double score = edge.getScore();
            if (score >= scoreThreshold) {

                cutForest.addEdge(edge, edge.getVertex1(), edge.getVertex2());
            }
        }
        return cutForest;
    }

    private Collection<Set<V>> getClusters(Forest<V,SimEdge<V>> forest) {
        WeakComponentClusterer<V, SimEdge<V>> weakComponentClusterer = new WeakComponentClusterer<V, SimEdge<V>>();
        Set<Set<V>> clusters = weakComponentClusterer.transform(forest);

        return clusters;
    }
}
