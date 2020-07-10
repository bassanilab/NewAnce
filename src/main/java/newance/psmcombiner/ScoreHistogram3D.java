/*
Copyright (C) SIB - Swiss Institute of Bioinformatics, Lausanne, Switzerland
Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
*/

package newance.psmcombiner;

import newance.psmconverter.PeptideSpectrumMatch;
import newance.util.NewAnceParams;

import java.io.*;
import java.util.*;

/**
 * @author Markus Müller
 */

public class ScoreHistogram3D extends SmoothedScoreHistogram {

    protected static List<List<Integer>> nnIndexMap = null;

    protected final double minScore1;
    protected final double maxScore1;
    protected final int nrScore1Bins;
    protected final double minScore2;
    protected final double maxScore2;
    protected final int nrScore2Bins;
    protected final double minScore3;
    protected final double maxScore3;
    protected final int nrScore3Bins;
    protected final double score1BinWidth;
    protected final double score2BinWidth;
    protected final double score3BinWidth;
    protected final List<Float> score1Mids;
    protected final List<Float> score2Mids;
    protected final List<Float> score3Mids;
    protected final String score1;
    protected final String score2;
    protected final String score3;


    public ScoreHistogram3D(ScoreHistogram3D scoreHistogram3D) {
        super(scoreHistogram3D);

        this.score1 = scoreHistogram3D.score1;
        this.score2 = scoreHistogram3D.score2;
        this.score3 = scoreHistogram3D.score3;

        this.minScore1 = scoreHistogram3D.minScore1;
        this.maxScore1 = scoreHistogram3D.maxScore1;
        this.nrScore1Bins = scoreHistogram3D.nrScore1Bins;
        this.minScore2 = scoreHistogram3D.minScore2;
        this.maxScore2 = scoreHistogram3D.maxScore2;
        this.nrScore2Bins = scoreHistogram3D.nrScore2Bins;
        this.minScore3 = scoreHistogram3D.minScore3;
        this.maxScore3 = scoreHistogram3D.maxScore3;
        this.nrScore3Bins = scoreHistogram3D.nrScore3Bins;

        this.score1BinWidth = scoreHistogram3D.score1BinWidth;
        this.score2BinWidth = scoreHistogram3D.score2BinWidth;
        this.score3BinWidth = scoreHistogram3D.score3BinWidth;

        this.score1Mids = scoreHistogram3D.score1Mids;
        this.score2Mids = scoreHistogram3D.score2Mids;
        this.score3Mids = scoreHistogram3D.score3Mids;

        if (nnIndexMap==null) nnIndexMap = initNNIndexMap();
        this.smoothedHistogram = null;
    }

    public ScoreHistogram3D(int[] nrBins, double minScore1, double maxScore1, int nrScore1Bins,
                            double minScore2, double maxScore2, int nrScore2Bins,
                            double minScore3, double maxScore3, int nrScore3Bins,
                            String score1, String score2, String score3) {
        super(nrBins);

        this.score1 = score1;
        this.score2 = score2;
        this.score3 = score3;

        this.minScore1 = minScore1;
        this.maxScore1 = maxScore1;
        this.nrScore1Bins = nrScore1Bins;
        this.minScore2 = minScore2;
        this.maxScore2 = maxScore2;
        this.nrScore2Bins = nrScore2Bins;
        this.minScore3 = minScore3;
        this.maxScore3 = maxScore3;
        this.nrScore3Bins = nrScore3Bins;

        this.score1BinWidth = (maxScore1 - minScore1)/ nrScore1Bins;
        this.score2BinWidth = (maxScore2 - minScore2)/ nrScore2Bins;
        this.score3BinWidth = (maxScore3 - minScore3)/ nrScore3Bins;

        this.score1Mids = calcMids(calcBreaks((float) minScore1,(float) maxScore1, nrScore1Bins));
        this.score2Mids = calcMids(calcBreaks((float) minScore2,(float) maxScore2, nrScore2Bins));
        this.score3Mids = calcMids(calcBreaks((float) minScore3,(float) maxScore3, nrScore3Bins));

        if (nnIndexMap==null) nnIndexMap = initNNIndexMap();
        this.smoothedHistogram = null;
    }

    @Override
    protected int index(PeptideSpectrumMatch peptideSpectrumMatch) {

        int score1Idx = get1DIndex(peptideSpectrumMatch.getScore(score1), minScore1, score1BinWidth, nrScore1Bins-1);
        int score2Idx = get1DIndex(peptideSpectrumMatch.getScore(score2), minScore2, score2BinWidth, nrScore2Bins-1);
        int score3Idx = get1DIndex(peptideSpectrumMatch.getScore(score3), minScore3, score3BinWidth, nrScore3Bins-1);

        int index = nrScore1Bins *(nrScore2Bins *score3Idx + score2Idx) + score1Idx;

        return index;
    }

    protected int index(double score1, double score2, double score3) {

        int score1Idx = get1DIndex(score1, minScore1, score1BinWidth, nrScore1Bins-1);
        int score2Idx = get1DIndex(score2, minScore2, score2BinWidth, nrScore2Bins-1);
        int score3Idx = get1DIndex(score3, minScore3, score3BinWidth, nrScore3Bins-1);

        int index = nrScore1Bins *(nrScore2Bins *score3Idx + score2Idx) + score1Idx;

        return index;
    }

    public void add(double score1, double score2, double score3, float value, boolean isDecoy) {

        if (value<=0) return;

        int bin = index(score1,score2,score3);

        int idx = indexMap.get(bin);
        if (idx<0) {
            indexMap.set(bin,currIndex);
            psmBins.add(bin);
            currIndex++;
        }

        if (idx < 0) { // new bin
            if (!isDecoy) {
                targetCnts.add(value);
                decoyCnts.add(0f);
            } else {
                targetCnts.add(0f);
                decoyCnts.add(value);
            }
        } else { // already recorded bin
            if (!isDecoy) targetCnts.set(idx,targetCnts.get(idx)+value);
            else decoyCnts.set(idx,decoyCnts.get(idx)+value);
        }

        if (isDecoy) totTargetCnt += value;
        else totDecoyCnt += value;
    }

    @Override
    protected Set<Integer> getNeighbourIndex(int bin) {

        List<Integer> scoreIdx = nnIndexMap.get(bin);

        int score1Idx = scoreIdx.get(0);
        int score2Idx = scoreIdx.get(1);
        int score3Idx = scoreIdx.get(2);

        Set<Integer> neighbours = new HashSet<>();

        if (score1Idx+1<nrScore1Bins) neighbours.add(nrScore1Bins *(nrScore2Bins *score3Idx + score2Idx) + score1Idx+1);
        if (score1Idx-1>=0) neighbours.add(nrScore1Bins *(nrScore2Bins *score3Idx + score2Idx) + score1Idx-1);
        if (score2Idx+1<nrScore2Bins) neighbours.add(nrScore1Bins *(nrScore2Bins *score3Idx + score2Idx+1) + score1Idx);
        if (score2Idx-1>=0) neighbours.add(nrScore1Bins *(nrScore2Bins *score3Idx + score2Idx-1) + score1Idx);
        if (score3Idx+1<nrScore3Bins) neighbours.add(nrScore1Bins *(nrScore2Bins *(score3Idx+1) + score2Idx) + score1Idx);
        if (score3Idx-1>=0) neighbours.add(nrScore1Bins *(nrScore2Bins *(score3Idx-1) + score2Idx) + score1Idx);

        return neighbours;
    }


    protected List<List<Integer>> initNNIndexMap(){

        List<List<Integer>> indexMap = new ArrayList<>(totNrBins);
        for (int i=0;i<totNrBins;i++) indexMap.add(new ArrayList<>(3));

        for (int score3Idx = 0; score3Idx < nrScore3Bins; score3Idx++) {
            for (int score2Idx = 0; score2Idx < nrScore2Bins; score2Idx++) {
                for (int score1Idx = 0; score1Idx < nrScore1Bins; score1Idx++) {
                    int index = nrScore1Bins *(nrScore2Bins *score3Idx + score2Idx) + score1Idx;

                    List<Integer> indexes = indexMap.get(index);
                    indexes.add(score1Idx);
                    indexes.add(score2Idx);
                    indexes.add(score3Idx);
                }
            }
        }
        return indexMap;
    }


    public void write(String outputDir, String fileTag, String id) {

        if (smoothedHistogram!=null)
            smoothedHistogram.write(new File(outputDir+File.separatorChar+fileTag+"_smoothed_"+id+".txt"));

        write(new File(outputDir+File.separatorChar+fileTag+"_"+id+".txt"));
    }

    protected void write(File outputFile) {

        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile));
            writer.write("#"+String.format("%.3f\t%.3f\t%d\t%.3f\t%.3f\t%d\t%.3f\t%.3f\t%d\n",
                    minScore1, maxScore1, nrScore1Bins, minScore2, maxScore2, nrScore2Bins, minScore3, maxScore3, nrScore3Bins));
            writer.write(score1+"\t"+score2+"\t"+score3+"\tValue\tType\n");

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

    public static ScoreHistogram3D read(File inputFile) {

        ScoreHistogram3D scoreHistogram3D = null;
        try {
            BufferedReader reader = new BufferedReader(new FileReader(inputFile));

            String line = reader.readLine();

            if (!line.startsWith("#")) {
                System.out.println("Ivalid histogram input file "+inputFile+". First line must start with #.");
                return null;
            }

            line = line.substring(1);
            String[] fields = line.split("\t");
            double minScore1 = Double.parseDouble(fields[0]);
            double maxScore1 = Double.parseDouble(fields[1]);
            int nrScore1Bins = Integer.parseInt(fields[2]);
            double minScore2 = Double.parseDouble(fields[3]);
            double maxScore2 = Double.parseDouble(fields[4]);
            int nrScore2Bins = Integer.parseInt(fields[5]);
            double minScore3 = Double.parseDouble(fields[6]);
            double maxScore3 = Double.parseDouble(fields[7]);
            int nrScore3Bins = Integer.parseInt(fields[8]);

            int[] nrBins = new int[3];
            nrBins[0] = nrScore1Bins;
            nrBins[1] = nrScore2Bins;
            nrBins[2] = nrScore3Bins;

            line = reader.readLine();
            fields = line.split("\t");

            scoreHistogram3D = new ScoreHistogram3D(nrBins, minScore1, maxScore1, nrScore1Bins,
                    minScore2, maxScore2, nrScore2Bins, minScore3, maxScore3, nrScore3Bins,
                    fields[0], fields[1], fields[2]);

            while ((line = reader.readLine()) != null) {

                if (line.trim().isEmpty()) continue;

                fields = line.split("\t");

                double score1 = Double.parseDouble(fields[0]);
                double score2 = Double.parseDouble(fields[1]);
                double score3 = Double.parseDouble(fields[2]);
                float value = Float.parseFloat(fields[3]);
                boolean isDecoy = fields[4].equals("decoy");

                scoreHistogram3D.add(score1, score2, score3, value, isDecoy);
            }

            reader.close();
        } catch (IOException e) {

        }

        return scoreHistogram3D;
    }


    protected List<Float> getMids(int bin) {

        int score1Idx = bin% nrScore1Bins;
        bin = (bin-score1Idx)/nrScore1Bins;
        int score2Idx = bin% nrScore2Bins;
        bin = (bin-score2Idx)/nrScore2Bins;
        int score3Idx = bin% nrScore3Bins;

        return Arrays.asList(new Float[]{score1Mids.get(score1Idx), score2Mids.get(score2Idx),
                score3Mids.get(score3Idx)});
    }

    public void removeSpikeNoiseHistogram(boolean adjustTotalCounts) {
        if (smoothedHistogram==null) {
            smoothedHistogram = new ScoreHistogram3D(this);
        }

        smoothedHistogram.removeSpikeNoise(adjustTotalCounts);
    }


    public void smoothHistogram(boolean adjustTotalCounts) {
        if (smoothedHistogram==null) {
            smoothedHistogram = new ScoreHistogram3D(this);
        }

        smoothedHistogram.smooth(adjustTotalCounts);
    }



}
