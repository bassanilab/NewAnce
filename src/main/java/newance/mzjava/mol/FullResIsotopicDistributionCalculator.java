package newance.mzjava.mol;

import newance.mzjava.ms.peaklist.DoublePeakList;
import newance.mzjava.ms.peaklist.PeakAnnotation;
import newance.mzjava.ms.peaklist.PeakList;

import java.security.InvalidParameterException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Set;

/**
 * Calculation of the isotopic distribution of a Composition. This method is exact and calculates all isotopic peaks at arbitrary
 * resolution. Peaks that are within a mass tolerance can be grouped and their probabilities summed up.
 *
 * @author Markus Muller
 * @version 0.0
 */
public class FullResIsotopicDistributionCalculator {

    public static final class Peak {

        public Peak(double mz, double prob){
            this.mz = mz;
            this.prob = prob;
        }

        public double mz;
        public double prob;
    }

    class PeakComparator implements Comparator<Peak> {

        @Override
        public int compare(Peak p1, Peak p2) {
            if (Math.abs(p1.mz-p2.mz)< massTol) return 0;
            if (p1.mz<p2.mz) return -1;

            return 1;
        }
    }

    protected double[][] isotopeMasses;
    protected double[][] isotopeProbs;
    protected double[][] logIsotopeProbs;
    protected int[][] nrNeutrons;

    protected int maxNrIsotopes;
    protected double minProb;
    protected double massTol;

    protected final double minIsoProb;  // avoid having probs==0.0
    protected final PeakComparator peakComparator;
    protected Peak[] isotopicPeakBuffer;
    protected int peakListIdx;
    protected int peakArrayIncrease;

    protected final int maxNrAtoms;
    protected final double[] logInts;
    protected double logFactorialN;
    protected final double[] logFactorials;

    public FullResIsotopicDistributionCalculator() {
        this.maxNrIsotopes = 0;
        this.minProb = 0.0;
        this.massTol = 0.0;
        this.minIsoProb = 0.0000001;
        this.peakComparator = new PeakComparator();

        this.maxNrAtoms = 10000;
        this.peakArrayIncrease = 100;

        this.logInts = new double[maxNrAtoms];
        for (int i=1;i<maxNrAtoms;i++) logInts[i] = Math.log(i);
        this.logFactorials = new double[maxNrAtoms];
        for (int i=0;i<maxNrAtoms;i++) logFactorials[i] = logfactorial(i);

        fetchIsotopeData();
    }

    public FullResIsotopicDistributionCalculator(int maxNrAtoms) {
        this.maxNrIsotopes = 0;
        this.minProb = 0.0;
        this.massTol = 0.0;
        this.minIsoProb = 0.0000001;
        this.peakComparator = new PeakComparator();

        this.maxNrAtoms = maxNrAtoms;
        this.peakArrayIncrease = 100;

        this.logInts = new double[maxNrAtoms];
        for (int i=1;i<maxNrAtoms;i++) logInts[i] = Math.log(i);
        this.logFactorials = new double[maxNrAtoms];
        for (int i=0;i<maxNrAtoms;i++) logFactorials[i] = logfactorial(i);

        fetchIsotopeData();
    }

    protected void fetchIsotopeData() {

        PeriodicTable periodicTable = PeriodicTable.getInstance();

        AtomicSymbol[] atoms = AtomicSymbol.values();
        isotopeMasses = new double[atoms.length][];
        isotopeProbs = new double[atoms.length][];
        logIsotopeProbs = new double[atoms.length][];
        nrNeutrons = new int[atoms.length][];

        for (AtomicSymbol atomicSymbol : AtomicSymbol.values()){

            List<Atom> isotopes = periodicTable.getAtoms(atomicSymbol);
            Atom monoIsotope = periodicTable.getAtom(atomicSymbol);

            int isoCnt = 0;
            for (Atom isotope : isotopes) if (isotope.getAbundance() >= minIsoProb) isoCnt++;

            int i = atomicSymbol.ordinal();
            isotopeMasses[i] = new double[isoCnt];
            isotopeProbs[i] = new double[isoCnt];
            logIsotopeProbs[i] = new double[isoCnt];
            nrNeutrons[i] = new int[isoCnt];

            int j = 0;
            for (Atom isotope : isotopes) {

                if (isotope.getAbundance()<minIsoProb) continue;

                isotopeMasses[i][j] = isotope.getMass();
                isotopeProbs[i][j] = isotope.getAbundance();
                logIsotopeProbs[i][j] = Math.log(isotope.getAbundance());
                nrNeutrons[i][j] = isotope.getMassNumber()-monoIsotope.getMassNumber();
                j++;
            }
        }
    }

    /**
     * Calculate isotopic distribution.
     *
     * @param composition   : chemical composition of molecule
     * @param charge        : charge of ions (masses are m/charge+1 if charge>0 or m if charge==0)
     * @param maxNrIsotopes : maximal number of additional neutrons considered
     * @param minProb       : only peaks with a probability larger than minProb are considered
     * @param massTol       : isotopic masses within massTol are merged and their probabilities summed up
     * @return Peaklist with mz of protonated ions (neutral masses if charge==0) and isotopic probabilities as intensities
     */

    public PeakList<PeakAnnotation> getIsotopicDistribution(Composition composition, int charge, int maxNrIsotopes, double minProb, double massTol) {

        this.maxNrIsotopes = maxNrIsotopes;
        this.minProb = minProb;
        this.massTol = massTol;

        Set<Atom> atoms = composition.getAtoms();

        Peak[] isoDist = null;
        for (Atom atom : atoms) {
            int cnt = composition.getCount(atom);

            if (cnt > maxNrAtoms) {
                throw new InvalidParameterException("Too many atoms. Increase maxNrAtoms using the FullResIsotopicDistributionCalculator constructor ");
            }

            if (isoDist == null) {
                isoDist = calcAtomIsoDist(atom, cnt);
            } else {
                isoDist = fold(isoDist, calcAtomIsoDist(atom, cnt));
            }
        }

        isoDist = discardSmallPeaks(isoDist);


        DoublePeakList<PeakAnnotation> peakList = new DoublePeakList<>(isoDist.length);

        if (isoDist != null) {
            for (int i = 0; i < isoDist.length; i++) {
                if (charge > 0)
                    peakList.add(isoDist[i].mz / charge + MassCalculator.PROTON_MASS, isoDist[i].prob);
                else
                    peakList.add(isoDist[i].mz, isoDist[i].prob);
            }
        }

        return peakList;
    }

    /**
     * Calculate isotopic distribution of an ion of a given mass with a chemical composition 'proportional' to typicalComp.
     * I.e. if the typical composition is "C40H72N12O13S2" (992.47832185 Th) and the mass is 500.000, the algorithm will
     * calculate the isotopic distribution of 'half' the composition ("C20H36N6O7S1") and set the monoisotopic peak to
     * 500.000.
     *
     * @param typicalComp  : typical chemical composition of molecule
     * @param charge : charge of ions (masses are m/charge+1 if charge>0 or m if charge==0)
     * @param maxNrIsotopes : maximal number of additional neutrons considered
     * @param minProb : only peaks with a probability larger than minProb are considered
     * @param massTol : isotopic masses within massTol are merged and their probabilities summed up
     * @return  Peaklist with mz of protonated ions (neutral masses if charge==0) and isotopic probabilities as intensities
     */
    public PeakList<PeakAnnotation> getIsotopicDistribution(double mass, Composition typicalComp, int charge, int maxNrIsotopes, double minProb, double massTol) {

        double typicalCompMass = typicalComp.getMolecularMass();

        double fact = mass/typicalCompMass;

        Composition.Builder compBuilder = new Composition.Builder();

        Set<Atom> atoms = typicalComp.getAtoms();

        for (Atom atom : atoms) {
            int cnt = typicalComp.getCount(atom);
            compBuilder.add(atom,(int) Math.round(cnt*fact));
        }

        return getIsotopicDistribution(compBuilder.build(),charge, maxNrIsotopes, minProb, massTol);
    }



    protected Peak[] calcAtomIsoDist(Atom atom, int totAtomCnt) {

        Peak[] isoPeaks;

        int i = atom.getSymbol().ordinal();

        if (nrNeutrons[i].length==1) {
            isoPeaks = new Peak[1];
            isoPeaks[0] = new Peak(totAtomCnt*isotopeMasses[i][0],1.0);

        } else  if (nrNeutrons[i].length==2) {

            isoPeaks = calcBinomialIsoDist(i,totAtomCnt);

        } else {
            logFactorialN = logFactorials[totAtomCnt];
            isoPeaks = calcMultinomialIsoDist(i,totAtomCnt);

        }

        return isoPeaks;
    }

    protected Peak[] calcBinomialIsoDist(int i, int totAtomCnt) {

        double isoMassDiff = isotopeMasses[i][nrNeutrons[i][1]]-isotopeMasses[i][0];
        int nr = (int)Math.round(maxNrIsotopes/isoMassDiff)+1;

        double[] binomProbs = calcBinomialProbs(totAtomCnt, isotopeProbs[i][nrNeutrons[i][1]], nr);

        Peak[] isoPeaks = new Peak[binomProbs.length];

        double monoIsoMass = isotopeMasses[i][0]*totAtomCnt;
        for (int j=0;j<binomProbs.length;j++) {
                isoPeaks[j] = new Peak(monoIsoMass+j*isoMassDiff,binomProbs[j]);
        }

        return isoPeaks;
    }

    protected Peak[] calcMultinomialIsoDist(int i, int totAtomCnt) {

        isotopicPeakBuffer = new Peak[maxNrIsotopes*(maxNrIsotopes+1)];  // initial size, grow if necessary later

        int[] isoCnts = new int[isotopeMasses[i].length];
        peakListIdx = 0;
        isoCnts[0] = totAtomCnt;
        calcCombinations(isotopeMasses[i], logIsotopeProbs[i], nrNeutrons[i], 0, totAtomCnt, isoCnts);

        isotopicPeakBuffer = Arrays.copyOfRange(isotopicPeakBuffer,0,peakListIdx);
        Arrays.sort(isotopicPeakBuffer,peakComparator);

        return mergePeaks(isotopicPeakBuffer);
    }

    protected void calcCombinations(double[] isoMasses, double[] isoLogProbs, int[] nrNeutrons, int isoIdx, int isoCnt, int[] isoCnts) {

        isoCnts[isoIdx] = isoCnt;

        int massdiff = 0;
        for (int i=1;i<=isoIdx;i++) massdiff += isoCnts[i]*nrNeutrons[i];

//        for (int i=0;i<isoCnts.length;i++) System.out.print(isoCnts[i]+",");
//        System.out.print(":"+isoIdx+" - "+isoCnt+" - "+totCnt+" - "+massdiff+" - "+peakListIdx+"\n");

        double mass = 0.0;
        for (int i=0;i<=isoIdx;i++) mass += isoMasses[i]*isoCnts[i];
        double prob = Math.exp(calcMultinomialLogProb(isoIdx, isoCnts, isoLogProbs));

        if (peakListIdx>=isotopicPeakBuffer.length)
            isotopicPeakBuffer = grow(isotopicPeakBuffer,peakArrayIncrease);

        isotopicPeakBuffer[peakListIdx] = new Peak(mass,prob);
        peakListIdx++;

        if (isoIdx==isoCnts.length-1) return;

        for (int i=1;i<=isoCnt;i++) {
            isoCnts[isoIdx] = isoCnt-i;

            if (massdiff+i>maxNrIsotopes) break;

            calcCombinations(isoMasses, isoLogProbs, nrNeutrons, isoIdx + 1, i, isoCnts);
        }
    }

    protected Peak[] fold(Peak[] isoDist1, Peak[] isoDist2) {

        Peak[] folded = new Peak[isoDist1.length*isoDist2.length];

        int k = 0;
        double monoIsotope = isoDist1[0].mz+isoDist2[0].mz;
        for (int i=0;i<isoDist1.length;i++) {
            for (int j=0;j<isoDist2.length;j++) {
                if (Math.round(isoDist1[i].mz+isoDist2[j].mz-monoIsotope)>maxNrIsotopes) break;

                folded[k] = new Peak(isoDist1[i].mz+isoDist2[j].mz,isoDist1[i].prob*isoDist2[j].prob);
                k++;
            }
        }

        folded = Arrays.copyOfRange(folded, 0, k);
        Arrays.sort(folded, peakComparator);

        return mergePeaks(folded);
    }


    protected Peak[] mergePeaks(Peak[] isoDist) {

        if (isoDist==null || isoDist.length==0) return null;

        boolean merged = false;
        int k=0;
        Peak prev = isoDist[k];
        double tmpMass = 0.0;
        for (int i=1;i<isoDist.length;i++) {
            Peak curr = isoDist[i];
            if (peakComparator.compare(curr,prev)!=0) {
                if (merged) {
                    tmpMass /= prev.prob;
                    prev.mz = tmpMass;
                    merged = false;
                }
                isoDist[k] = prev;
                k++;
                prev = curr;
            } else {
                if (!merged) {
                    tmpMass = prev.mz*prev.prob;
                }
                tmpMass += curr.prob*curr.mz;
                prev.prob += curr.prob;
                merged = true;
            }
        }
        if (merged)  {
            tmpMass /= prev.prob;
            prev.mz = tmpMass;
        }
        isoDist[k] = prev;

        return Arrays.copyOfRange(isoDist, 0, k + 1);
    }

    protected Peak[] discardSmallPeaks(Peak[] isoDist) {

        if (isoDist==null || isoDist.length==0) return null;

        int k = 0;
        for (int i=0;i<isoDist.length;i++) {
            Peak curr = isoDist[i];
            if (curr.prob>=minProb) {
                if (k<i) isoDist[k] = isoDist[i];
                k++;
            }
        }

        return Arrays.copyOfRange(isoDist, 0, k);
    }

    protected double[] calcBinomialProbs(int n, double p, int maxNrIsoPeaks) {

        double[] probs = new double[maxNrIsoPeaks];
        double p2 = 1.0-p;

        double tmp = p2;
        for (int i=1;i<n;i++) tmp *= p2;
        probs[0] = tmp;
        for (int k = 0; k < maxNrIsoPeaks-1; k++)
            probs[k+1] = (probs[k] * p * (n-k)) / ((k+1)*p2);

        return probs;
    }

    protected double calcMultinomialLogProb(int isoIdx, int[] isoCnts, double[] isoLogProbs){

        double logProb = 0.0;
        for (int i=0;i<=isoIdx;i++) logProb += isoCnts[i]*isoLogProbs[i]-logFactorials[isoCnts[i]];

        return logProb+logFactorialN;
    }

    protected double logfactorial(int n) {

        double res = 0.0;

        for (int i=1;i<=n;i++) res += logInts[i];

        return res;
    }

    protected Peak[] grow(Peak[] oldPeaks, int increase) {
        Peak[] newPeaks = new Peak[oldPeaks.length+increase];
        for (int i=0;i<oldPeaks.length;i++) newPeaks[i] = oldPeaks[i];

        return newPeaks;
    }

}
