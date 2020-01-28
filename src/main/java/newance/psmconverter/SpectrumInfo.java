/**
 * Copyright (C) 2019, SIB/LICR. All rights reserved
 *
 * SIB, Swiss Institute of Bioinformatics
 * Ludwig Institute for Cancer Research (LICR)
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer. Redistributions in binary
 * form must reproduce the above copyright notice, this list of conditions and
 * the following disclaimer in the documentation and/or other materials provided
 * with the distribution. Neither the name of the SIB/LICR nor the names of
 * its contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL SIB/LICR BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
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
