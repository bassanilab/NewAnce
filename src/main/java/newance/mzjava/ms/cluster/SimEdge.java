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

import static com.google.common.base.Preconditions.checkNotNull;

/**
 * Class to hols an undirected edge that is weighted by a score between two vertexes.
 * <p>
 * Because SimEdge is undirected edge A1-A2:0.9 and A2-A1:0.9 are equal and have the same hash code
 * <p>
 * This class is immutable
 *
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class SimEdge<V> {

    private final V vertex1;
    private final V vertex2;
    private final double score;

    /**
     * Constructor
     *
     * @param vertex1 the first vertex
     * @param vertex2 the second vertex
     * @param score the score of the vertex
     */
    public SimEdge(V vertex1, V vertex2, double score) {

        checkNotNull(vertex1);
        checkNotNull(vertex2);

        if (vertex1.equals(vertex2))
            throw new IllegalStateException("SimEdges cannot be a loop. The loop was on " + vertex1 + " - " + vertex2);

        this.vertex1 = vertex1;
        this.vertex2 = vertex2;
        this.score = score;
    }

    public V getVertex1() {

        return vertex1;
    }

    public V getVertex2() {

        return vertex2;
    }

    public double getScore() {

        return score;
    }

    @Override
    public String toString() {

        return "SimEdge{" +
                vertex1 +
                " - " + vertex2 +
                ", score=" + score +
                '}';
    }

    public V getOther(V vertex) {

        if(vertex.equals(vertex1))
            return vertex2;
        else if(vertex.equals(vertex2)) {
            return vertex1;
        } else {
            throw new IllegalStateException(vertex + " is not a vertex in this edge");
        }
    }

    @Override
    public boolean equals(Object o) {

        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        SimEdge<V> simEdge = (SimEdge<V>) o;

        return Double.compare(simEdge.score, score) == 0 &&
                contains(simEdge.vertex1) &&
                contains(simEdge.vertex2);
    }

    public boolean contains(V vertex){

        return vertex1.equals(vertex) || vertex2.equals(vertex);
    }

    @Override
    public int hashCode() {

        int result;
        long temp;
        result = vertex1.hashCode() + vertex2.hashCode();
        temp = Double.doubleToLongBits(score);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        return result;
    }
}
