package newance.mzjava.ms.cluster;

import newance.mzjava.mol.modification.Tolerance;
import newance.mzjava.ms.peaklist.PeakAnnotation;
import newance.mzjava.ms.peaklist.PeakList;
import newance.mzjava.ms.spectrasim.SimFunc;

import java.util.*;

/**
 * Class to evaluate whether a set of spectra cluster together. If the set contains alien spectra,
 * the SpectrumClusterEvaluator.evaluate method tries to find the largest subset of similar spectra.
 * If no such subset can be found, it returns the spectrum with the largest similarity to the other
 * set members.
 *
 * @author Markus Muller
 * @version 0.0
 */
public class SpectrumClusterEvaluator<A extends PeakAnnotation, S extends PeakList<A>> {

    private final Tolerance precTolerance;
    private final SimFunc<A, A> simFunc;
    private final ClusterBuilder<S> clusterBuilder;
    private final Comparator<S> spectrumComparator;

    protected double minSimScore;
    protected double maxSimScore;
    protected double meanSimScore;
    protected double sdSimScore;
    protected int nrMembersOrg;
    protected int nrMembersCleaned;


    public SpectrumClusterEvaluator(Tolerance tolerance, SimFunc<A, A> simFunc, ClusterBuilder<S> clusterBuilder) {

        this.precTolerance = tolerance;
        this.simFunc = simFunc;
        this.clusterBuilder = clusterBuilder;
        this.spectrumComparator = new Comparator<S>() {
            @Override
            public int compare(S spectrum1, S spectrum2) {

                return (spectrum1.getPrecursor().compareTo(spectrum2.getPrecursor()));
            }
        };
    }

    public Set<S> evaluate(Collection<S> spectra) {

        if (spectra.isEmpty())
            return new HashSet<>();

        List<S> sortedSpectra = new ArrayList<>(spectra);
        Collections.sort(sortedSpectra, spectrumComparator);
        nrMembersOrg = spectra.size();

        SimilarityGraph<S> simGraph = PeakListSimGraphFactory.build(sortedSpectra, simFunc, precTolerance, 0.0);

        Map<S, Double> totClusterSims = calcSimScoreStats(simGraph);

        Collection<Set<S>> clusters = clusterBuilder.cluster(simGraph);

        Set<S> maxCluster = getBestSubset(clusters, totClusterSims);

        // recalculate cluster stats
        sortedSpectra = new ArrayList<>(maxCluster);
        Collections.sort(sortedSpectra, spectrumComparator);
        simGraph = PeakListSimGraphFactory.build(sortedSpectra, simFunc, precTolerance, 0.0);
        calcSimScoreStats(simGraph);

        nrMembersCleaned = maxCluster.size();

        return maxCluster;
    }

    private Set<S> getBestSubset(Collection<Set<S>> clusters, Map<S, Double> totClusterSims) {
        // get largest cluster
        double maxClusterScore = 0.0;
        Set<S> maxCluster = null;
        for (Set<S> cluster : clusters) {
            double clusterScore = 0.0;
            for (S spectrum : cluster) {
                clusterScore += totClusterSims.get(spectrum);
            }
            if (clusterScore >= maxClusterScore) {
                maxClusterScore = clusterScore;
                maxCluster = cluster;
            }
        }

        return maxCluster;
    }

    private Map<S, Double> calcSimScoreStats(SimilarityGraph<S> simGraph) {

        meanSimScore = 0.0;
        sdSimScore = 0.0;
        minSimScore = 0.0;
        maxSimScore = 0.0;
        int n = 0;

        Map<S, Double> totClusterSims = new HashMap<>();

        for (S vertex : simGraph.getVertices()) {
            totClusterSims.put(vertex, 0.0);
        }

        if (simGraph.getEdgeCount() > 0) {
            minSimScore = Double.MAX_VALUE;
            maxSimScore = Double.MIN_VALUE;
        }

        for (SimEdge<S> edge : simGraph.getEdges()) {
            n++;

            double score = edge.getScore();
            double delta = score - meanSimScore;
            meanSimScore += delta / n;
            sdSimScore += delta * (score - meanSimScore);

            if (score < minSimScore) minSimScore = score;
            if (score > maxSimScore) maxSimScore = score;

            totClusterSims.put(edge.getVertex1(), totClusterSims.get(edge.getVertex1()) + score);
            totClusterSims.put(edge.getVertex2(), totClusterSims.get(edge.getVertex2()) + score);
        }

        if (n > 1) sdSimScore /= n - 1;

        sdSimScore = Math.sqrt(sdSimScore);

        return totClusterSims;
    }

    public double getMinSimScore() {

        return minSimScore;
    }

    public double getMaxSimScore() {

        return maxSimScore;
    }

    public double getMeanSimScore() {

        return meanSimScore;
    }

    public double getSdSimScore() {

        return sdSimScore;
    }

    public int getNrMembersOrg() {

        return nrMembersOrg;
    }

    public int getNrMembersCleaned() {

        return nrMembersCleaned;
    }
}
