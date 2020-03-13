/*
Copyright (C) SIB - Swiss Institute of Bioinformatics, Lausanne, Switzerland
Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
*/

package newance.psmconverter;


/**
 * @author Markus Müller
 */

public class SpectrumInfo {

    private final String spectrumFile;
    private final String spectrum;
    private int scanNumber;
    private double retentionTime;
    private double precursorNeutralMass;
    private double precursorIntensity;
    private double precursorMz;
    private int charge;
    private int index;

    public SpectrumInfo(String spectrum) {
        this.spectrum = spectrum;
        this.spectrumFile = spectrum.substring(0,spectrum.indexOf("."));
    }

    public String getSpectrum() {
        return spectrum;
    }

    public int getScanNumber() {
        return scanNumber;
    }

    public void setScanNumber(int scanNumber) {
        this.scanNumber = scanNumber;
    }

    public double getRetentionTime() {
        return retentionTime;
    }

    public void setRetentionTime(double retentionTime) {
        this.retentionTime = retentionTime;
    }

    public double getPrecursorNeutralMass() {
        return precursorNeutralMass;
    }

    public void setPrecursorNeutralMass(double precursorNeutralMass) {
        this.precursorNeutralMass = precursorNeutralMass;
    }

    public double getPrecursorIntensity() {
        return precursorIntensity;
    }

    public void setPrecursorIntensity(double precursorIntensity) {
        this.precursorIntensity = precursorIntensity;
    }

    public double getPrecursorMz() {
        return precursorMz;
    }

    public void setPrecursorMz(double precursorMz) {
        this.precursorMz = precursorMz;
    }

    public int getCharge() {
        return charge;
    }

    public void setCharge(int charge) {
        this.charge = charge;
    }

    public int getIndex() {
        return index;
    }

    public void setIndex(int index) {
        this.index = index;
    }

    public String getSpectrumFile() {
        return spectrumFile;
    }
}
