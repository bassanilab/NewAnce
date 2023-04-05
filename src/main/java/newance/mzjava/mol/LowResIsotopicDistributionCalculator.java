package newance.mzjava.mol;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.apache.commons.math3.util.FastMath;

import java.util.Arrays;
import java.util.List;
import java.util.Set;

/**
 * Calculation of the isotopic probabilities of a Composition, but not the isotopic masses. This method groups all isotopic peaks with the same additional
 * number of neutrons into one peak and calculates the total probability for each peak. The method uses fast fourier transformation
 * to calculate the convolution products. The method performs significantly faster (3-5x) compared to the FullResIsotopicDistributionCalculator.
 * It can be used if the accurate (<10 ppm) isotopic masses are not important or are known otherwise.
 *
 * @author Markus Muller
 * @version 0.0
 */
public class LowResIsotopicDistributionCalculator {


    protected class ComplexPolar{
        // Complex.multiply needs 4 double multiplications and 2 additions. If the complex number is
        // represented in polar ccordinates (angle, radius) a multiplication only one complex and one
        // integer multiplication, which makes the algorithm perform significantly faster.

        private double angle;
        private double radius;

        public ComplexPolar(double angle, double radius) {
            this.angle = angle;
            this.radius = radius;
        }

        public ComplexPolar(Complex c) {
            this.angle = c.getArgument();
            this.radius = c.abs();
        }

        public Complex castToComplex() {
            return new Complex(radius*FastMath.cos(angle),radius* FastMath.sin(angle));
        }

        public void pow(int n) {
            double tmp = radius;
            for (int i=1;i<n;i++) tmp *= radius;
            radius = tmp;
            angle *= n;
        }

        public void multiply(ComplexPolar complexPolar) {
            radius *= complexPolar.radius;
            angle += complexPolar.angle;
        }

        public ComplexPolar copy(){
            return new ComplexPolar(angle,radius);
        }
    }

    protected double[][] isotopeMasses;
    protected double[][] isotopeFrequencies;
    protected final int maxNrIsotopes;
    protected Complex[][] ftIsotopeFrequencies;
    protected ComplexPolar[][] ftIsotopeFrequenciesPolar;
    protected final FastFourierTransformer fastFourierTransformer;


    public LowResIsotopicDistributionCalculator() {

        maxNrIsotopes = 16;
        fastFourierTransformer = new FastFourierTransformer(DftNormalization.STANDARD);

        fetchIsotopeData();

        FFTIsotopeFrequencies();
    }

    /**
     * Constructor
     * @param maxNrIsotopes is the maximal number of isotopes considered. The is set to the next power of 2
     */
    public LowResIsotopicDistributionCalculator(int maxNrIsotopes) {

        if (( maxNrIsotopes & (maxNrIsotopes - 1)) == 0 && maxNrIsotopes > 0)   // check if power of 2
            this.maxNrIsotopes = maxNrIsotopes;
        else {

            int i;
            for (i = 1;i < maxNrIsotopes; i *= 2);
            this.maxNrIsotopes = i;
        }

        fastFourierTransformer = new FastFourierTransformer(DftNormalization.STANDARD);

        fetchIsotopeData();

        FFTIsotopeFrequencies();
    }

    protected void fetchIsotopeData()
    {
        PeriodicTable periodicTable = PeriodicTable.getInstance();

        AtomicSymbol[] atoms = AtomicSymbol.values();
        isotopeMasses = new double[atoms.length][];
        isotopeFrequencies = new double[atoms.length][];
        ftIsotopeFrequencies = new Complex[atoms.length][];

        for (AtomicSymbol atomicSymbol : AtomicSymbol.values()){

            List<Atom> isotopes = periodicTable.getAtoms(atomicSymbol);
            Atom monoIsotope = periodicTable.getAtom(atomicSymbol);

            int i = atomicSymbol.ordinal();
            isotopeMasses[i] = new double[maxNrIsotopes];
            isotopeFrequencies[i] = new double[maxNrIsotopes];

            Arrays.fill(isotopeFrequencies[i], 0.0);

            for (Atom isotope : isotopes) {
                int nrAdditionalNeutrons = isotope.getMassNumber()-monoIsotope.getMassNumber();

                if (nrAdditionalNeutrons>=maxNrIsotopes) continue;

                isotopeMasses[i][nrAdditionalNeutrons] = isotope.getMass();
                isotopeFrequencies[i][nrAdditionalNeutrons] = isotope.getAbundance();
            }
        }
    }

    /**
     * Calculate isotopic distribution.
     * @param composition  : chemical composition of molecule
     * @return  Peaklist with mz of protonated ions (neutral masses if charge==0) and isotopic probabilities as intensities
     */
    public double[] getIsotopicDistribution(Composition composition) {

        return calcIsotopicDistribution(composition);

    }

    /**
     * Calculate isotopic distribution of an ion of a given mass with a chemical composition 'proportional' to typicalComp.
     * I.e. if the typical composition is "C40H72N12O13S2" (992.47832185 Th) and the mass is 500.000, the algorithm will
     * calculate the isotopic distribution of 'half' the composition ("C20H36N6O7S1").
     *
     * @param typicalComp  : typical chemical composition of molecule
     * @return  Peaklist with mz of protonated ions (neutral masses if charge==0) and isotopic probabilities as intensities
     */
    public double[] getIsotopicDistribution(double mass, Composition typicalComp) {

        double typicalCompMass = typicalComp.getMolecularMass();

        double fact = mass/typicalCompMass;

        Composition.Builder compBuilder = new Composition.Builder();

        Set<Atom> atoms = typicalComp.getAtoms();

        for (Atom atom : atoms) {
            int cnt = typicalComp.getCount(atom);
            compBuilder.add(atom, (int) Math.round(cnt * fact));
        }

        return getIsotopicDistribution(compBuilder.build());
    }

    protected void FFTIsotopeFrequencies() {

        for (int i=0;i<isotopeFrequencies.length;i++) {
            ftIsotopeFrequencies[i] = fastFourierTransformer.transform(isotopeFrequencies[i], TransformType.FORWARD);
        }

        ftIsotopeFrequenciesPolar = new ComplexPolar[ftIsotopeFrequencies.length][];

        for (int i=0;i<isotopeFrequencies.length;i++) {
            ftIsotopeFrequenciesPolar[i] = new ComplexPolar[ftIsotopeFrequencies[i].length];
            for (int j=0;j<ftIsotopeFrequencies[i].length;j++)
                ftIsotopeFrequenciesPolar[i][j] = new ComplexPolar(ftIsotopeFrequencies[i][j]);
        }
    }


    protected double[] calcIsotopicDistribution(Composition composition) {

        Set<Atom> atoms = composition.getAtoms();

        Complex[] fftIsoDist = new Complex[maxNrIsotopes];

        for (int j=0;j<maxNrIsotopes;j++) {
            ComplexPolar c = new ComplexPolar(0.0,1.0);
            for (Atom atom : atoms) {
                ComplexPolar cAtom = ftIsotopeFrequenciesPolar[atom.getSymbol().ordinal()][j].copy();
                cAtom.pow(composition.getCount(atom));
                c.multiply(cAtom);
            }

            fftIsoDist[j] = c.castToComplex();
        }

        fftIsoDist = fastFourierTransformer.transform(fftIsoDist,TransformType.INVERSE);
        double[] isoDist = new double[maxNrIsotopes];
        for (int j=0;j<maxNrIsotopes;j++) {
            isoDist[j] = fftIsoDist[j].getReal();
        }

        return isoDist;
    }

}
