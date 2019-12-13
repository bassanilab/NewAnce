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

package newance.psmcombiner;

import newance.psmconverter.PeptideMatchData;
import newance.util.NewAnceParams;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * @author Markus Müller
 */

public class CometScoreHistogram extends SmoothedScoreHistogram {

    protected static List<List<Integer>> nnIndexMap = null;

    protected final double minXCorr;
    protected final double maxXCorr;
    protected final int nrXCorrBins;
    protected final double minDeltaCn;
    protected final double maxDeltaCn;
    protected final int nrDeltaCnBins;
    protected final double minSpScore;
    protected final double maxSpScore;
    protected final int nrSpScoreBins;
    protected final double xCorrBinWidth;
    protected final double deltaCnBinWidth;
    protected final double spScoreBinWidth;
    protected final List<Float> xCorrMids;
    protected final List<Float> deltaCnMids;
    protected final List<Float> spScoreMids;

    public CometScoreHistogram(int[] nrBins) {
        super(nrBins);

        this.minXCorr = NewAnceParams.getInstance().getMinXCorr();
        this.maxXCorr = NewAnceParams.getInstance().getMaxXCorr();
        this.nrXCorrBins = NewAnceParams.getInstance().getNrXCorrBins();
        this.minDeltaCn = NewAnceParams.getInstance().getMinDeltaCn();
        this.maxDeltaCn = NewAnceParams.getInstance().getMaxDeltaCn();
        this.nrDeltaCnBins = NewAnceParams.getInstance().getNrDeltaCnBins();
        this.minSpScore = NewAnceParams.getInstance().getMinSpScore();
        this.maxSpScore = NewAnceParams.getInstance().getMaxSpScore();
        this.nrSpScoreBins = NewAnceParams.getInstance().getNrSpScoreBins();

        this.xCorrBinWidth = (maxXCorr-minXCorr)/nrXCorrBins;
        this.deltaCnBinWidth = (maxDeltaCn-minDeltaCn)/nrDeltaCnBins;
        this.spScoreBinWidth = (maxSpScore-minSpScore)/nrSpScoreBins;

        this.xCorrMids = calcMids(calcBreaks((float)minXCorr,(float)maxXCorr,nrXCorrBins));
        this.deltaCnMids = calcMids(calcBreaks((float)minDeltaCn,(float)maxDeltaCn,nrDeltaCnBins));
        this.spScoreMids = calcMids(calcBreaks((float)minSpScore,(float)maxSpScore,nrSpScoreBins));

        if (nnIndexMap==null) nnIndexMap = initNNIndexMap();
        this.smoothedHistogram = null;
    }

    public CometScoreHistogram(ScoreHistogram cometScoreHistogram) {

        super(cometScoreHistogram);

        this.minXCorr = NewAnceParams.getInstance().getMinXCorr();
        this.maxXCorr = NewAnceParams.getInstance().getMaxXCorr();
        this.nrXCorrBins = NewAnceParams.getInstance().getNrXCorrBins();
        this.minDeltaCn = NewAnceParams.getInstance().getMinDeltaCn();
        this.maxDeltaCn = NewAnceParams.getInstance().getMaxDeltaCn();
        this.nrDeltaCnBins = NewAnceParams.getInstance().getNrDeltaCnBins();
        this.minSpScore = NewAnceParams.getInstance().getMinSpScore();
        this.maxSpScore = NewAnceParams.getInstance().getMaxSpScore();
        this.nrSpScoreBins = NewAnceParams.getInstance().getNrSpScoreBins();

        this.xCorrBinWidth = (maxXCorr-minXCorr)/nrXCorrBins;
        this.deltaCnBinWidth = (maxDeltaCn-minDeltaCn)/nrDeltaCnBins;
        this.spScoreBinWidth = (maxSpScore-minSpScore)/nrSpScoreBins;

        this.xCorrMids = calcMids(calcBreaks((float)minXCorr,(float)maxXCorr,nrXCorrBins));
        this.deltaCnMids = calcMids(calcBreaks((float)minDeltaCn,(float)maxDeltaCn,nrDeltaCnBins));
        this.spScoreMids = calcMids(calcBreaks((float)minSpScore,(float)maxSpScore,nrSpScoreBins));

        if (nnIndexMap==null) nnIndexMap = initNNIndexMap();
        this.smoothedHistogram = null;
    }

    @Override
    protected int index(PeptideMatchData peptideMatchData) {

        int xcorrIdx = get1DIndex(peptideMatchData.getScore("xcorr"), minXCorr, xCorrBinWidth, nrXCorrBins-1);
        int deltacnIdx = get1DIndex(peptideMatchData.getScore("deltacn"), minDeltaCn, deltaCnBinWidth, nrDeltaCnBins-1);
        int spscoreIdx = get1DIndex(peptideMatchData.getScore("spscore"), minSpScore, spScoreBinWidth, nrSpScoreBins-1);

        int index = nrXCorrBins*(nrDeltaCnBins*spscoreIdx + deltacnIdx) + xcorrIdx;

        return index;
    }

    @Override
    protected Set<Integer> getNeighbourIndex(int bin) {

        List<Integer> scoreIdx = nnIndexMap.get(bin);

        int xcorrIdx = scoreIdx.get(0);
        int deltacnIdx = scoreIdx.get(1);
        int spscoreIdx = scoreIdx.get(2);

        Set<Integer> neighbours = new HashSet<>();

        if (xcorrIdx+1<nrXCorrBins) neighbours.add(nrXCorrBins*(nrDeltaCnBins*spscoreIdx + deltacnIdx) + xcorrIdx+1);
        if (xcorrIdx-1>=0) neighbours.add(nrXCorrBins*(nrDeltaCnBins*spscoreIdx + deltacnIdx) + xcorrIdx-1);
        if (deltacnIdx+1<nrDeltaCnBins) neighbours.add(nrXCorrBins*(nrDeltaCnBins*spscoreIdx + deltacnIdx+1) + xcorrIdx);
        if (deltacnIdx-1>=0) neighbours.add(nrXCorrBins*(nrDeltaCnBins*spscoreIdx + deltacnIdx-1) + xcorrIdx);
        if (spscoreIdx+1<nrSpScoreBins) neighbours.add(nrXCorrBins*(nrDeltaCnBins*(spscoreIdx+1) + deltacnIdx) + xcorrIdx);
        if (spscoreIdx-1>=0) neighbours.add(nrXCorrBins*(nrDeltaCnBins*(spscoreIdx-1) + deltacnIdx) + xcorrIdx);

        return neighbours;
    }


    protected List<List<Integer>> initNNIndexMap(){

        List<List<Integer>> indexMap = new ArrayList<>(totNrBins);
        for (int i=0;i<totNrBins;i++) indexMap.add(new ArrayList<>(3));

        for (int spscoreIdx = 0; spscoreIdx < nrSpScoreBins; spscoreIdx++) {
            for (int deltacnIdx = 0; deltacnIdx < nrDeltaCnBins; deltacnIdx++) {
                for (int xcorrIdx = 0; xcorrIdx < nrXCorrBins; xcorrIdx++) {
                    int index = nrXCorrBins*(nrDeltaCnBins*spscoreIdx + deltacnIdx) + xcorrIdx;

                    List<Integer> indexes = indexMap.get(index);
                    indexes.add(xcorrIdx);
                    indexes.add(deltacnIdx);
                    indexes.add(spscoreIdx);
                }
            }
        }
        return indexMap;
    }


    public void write(File outputFile) {

        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile));
            writer.write("XCorr\tDeltaCn\tSpScore\tValue\tType\n");

            for (int i=0;i<indexMap.size();i++) {

                List<Float> mids = getMids(i);
                String coords = "";
                for (Float mid : mids) coords+= String.format("%.3f",mid)+"\t";

                float value = 0f;
                if (indexMap.get(i)>=0) {
                    int idx = indexMap.get(i);
                    writer.write(coords+String.format("%.3f",targetCnts.get(idx))+"\ttarget\n");
                    writer.write(coords+String.format("%.3f",decoyCnts.get(idx))+"\tdecoy\n");
                    if (!lFDR.isEmpty()) value = lFDR.get(idx); // check whether lFDR is already calculated
                    writer.write(coords+String.format("%.3f",value)+"\tlFDR\n");
                } else {
                    writer.write(coords+String.format("%.3f",value)+"\ttarget\n");
                    writer.write(coords+String.format("%.3f",value)+"\tdecoy\n");
                    writer.write(coords+String.format("%.3f",value)+"\tlFDR\n");
                }
            }

            writer.close();
        } catch (IOException e) {

        }
    }

    protected List<Float> getMids(int bin) {

        int xCorrIdx = bin%nrXCorrBins;
        bin = (bin-xCorrIdx)/nrXCorrBins;
        int deltaCnIdx = bin%nrDeltaCnBins;
        bin = (bin-deltaCnIdx)/nrDeltaCnBins;
        int spScoreIdx = bin%nrSpScoreBins;

        return Arrays.asList(new Float[]{xCorrMids.get(xCorrIdx),deltaCnMids.get(deltaCnIdx),spScoreMids.get(spScoreIdx)});
    }

}
