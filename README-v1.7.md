# NewAnce Version 1.7

This version allows to calculate a group specific FDR for both MaxQuant and Comet. It also includes better variant annotations supporting complex variants. 
The code has been tested on Linux, Mac, and Windows using Java 1.8. Comet result previous to consensus filtering can be exported.
 
## TOC

1. [Command line options](#Command-line-options)
2. [Running NewAnce](#Running-NewAnce)
3. [Output format](#Output-format)
4. [NewAnce with Comet only](#NewAnce-with-Comet-only)
5. [Tests](#Tests)
6. [PDV export](#PDV-export)
7. [Extend MaxQuant features](#Extend-MaxQuant-features)
8. [Export score histograms](#Export-score-histograms)
9. [Spectronaut DIA export](#Spectronaut-DIA-export)
10. [Match peptides to fasta file](#Match-peptides-to-fasta-file)

## Command line options 

The NewAnce command line options can be obtained by typing

```
java -jar NewAnce-1.7.3-SNAPSHOT.jar -h
```
which write the command line option infor to standard output :

```
*************************************************************************************************************************
**                                            NewAnce command line help                                                **
*************************************************************************************************************************


usage: newance.psmcombiner.CometMaxQuantCombiner
 -coD,--cometPsmDir <arg>           Comet psm root directory (required)
 -coFDR,--cometFDR <arg>            FDR for filtering Comet PSMs before combination (default value 0.03)
 -coRE,--cometPsmRegex <arg>        Regular expression of Comet psm files (e.g. \.xml$) (required)
 -d,--debug                         Debug option
 -exclP,--excludeProts <arg>        Regular expression of proteins excluded from analysis. If not set no proteins are excluded.
 -fdrM,--fdrControlMethod <arg>     Method to control pFDR: combined or separate (default combined).
 -fH,--forceHistograms              Histograms are imported even if enough PSMs are available.
 -groupF,--groupProteinFile <arg>   Tab file with protein group assignments which will override assignment by groupRE
 -groupM,--groupingMethod <arg>     Method for PSM grouping: "fasta", "modif", "famo" (fasta&modif) or "none" (default none).
 -groupN,--groupNames <arg>         Comma separated list of names of sequence groups in fasta file (e.g. prot,lncRNA,TE ). Will be used as prefixes
                                    for output files.
 -groupRE,--groupRegEx <arg>        Comma separated list of regular expression defining sequence groups of fasta headers (e.g.
                                    "sp\||tr\|ENSP00","ENST00","SINE_|LINE_|LTR_|DNA_|Retroposon_" ). Will be used as prefixes for output files.
 -h,--help                          Help option for command line help
 -maxDC,--maxDeltaCn <arg>          Maximal Comet DeltaCn in histogram (default value 2500)
 -maxDS,--maxDeltaScore <arg>       Maximal MaxQuant DeltaScore in histogram (default value 1)
 -maxL,--maxLength <arg>            Maximal length of peptide (default value: 25)
 -maxPEP,--maxPEP <arg>             Maximal MaxQuant PEP in histogram (default value 2500)
 -maxR,--maxRank <arg>              Maximal rank of peptide in list of spectrum matches (rank 1 = best) (default value: 1)
 -maxScore,--maxScore <arg>         Maximal MaxQuant Score in histogram (default value 5)
 -maxSP,--maxSpScore <arg>          Maximal Comet SpScore in histogram (default value 1)
 -maxXC,--maxXCorr <arg>            Maximal Comet XCorr in histogram (default value 5)
 -maxZ,--maxCharge <arg>            Maximal charge of PSM (default value: 5)
 -minDC,--minDeltaCn <arg>          Minimal Comet DeltaCn in histogram (default value 0)
 -minDS,--minDeltaScore <arg>       Minimal MaxQuant DeltaScore in histogram (default value 0)
 -minL,--minLength <arg>            Minimal length of peptide (default value: 8)
 -minPEP,--minPEP <arg>             Minimal MaxQuant PEP in histogram (default value 0)
 -minPH,--minPsm4Histo <arg>        Minimal number of psms to calculate local FDR in histogram (default value: 100000).
 -minScore,--minScore <arg>             Minimal MaxQuant Score in histogram (default value 0)
 -minSP,--minSpScore <arg>              Minimal Comet SpScore in histogram (default value 0)
 -minXC,--minXCorr <arg>                Minimal Comet XCorr in histogram (default value 0)
 -minZ,--minCharge <arg>                Minimal charge of PSM (default value: 1)
 -mod,--modifications <arg>             Comma separated list of peptide modifications used in search (e.g. Cysteinyl:C3H5NO2S,Oxidation:O)
 -mqD,--maxquantPsmDir <arg>            MaxQuant psm root directory.
 -mqFDR,--maxquantFDR <arg>             FDR for filtering MaxQuant PSMs before combination (default value 0.03)
 -nrDCB,--nrDeltaCnBins <arg>           Number of Comet DeltaCn bins in histogram (default value 40)
 -nrDSB,--nrDeltaScoreBins <arg>        Number of MaxQuant DeltaScore bins in histogram (default value 40)
 -nrPEPB,--nrPEPBins <arg>              Number of MaxQuant PEP bins in histogram (default value 40)
 -nrScoreB,--nrScoreBins <arg>          Number of MaxQuant Score bins in histogram (default value 40)
 -nrSPB,--nrSpScoreBins <arg>           Number of Comet SpScore bins in histogram (default value 40)
 -nrTh,--nrThreads <arg>                Number of threads used by NewAnce (default value: nr of available processors - 2)
 -nrXCB,--nrXCorrBins <arg>             Number of Comet XCorr bins in histogram (default value 40)
 -outD,--outputDir <arg>                Output directory for results (required)
 -outT,--outputTag <arg>                Tag inserted into output file names after prefix.
 -ppG,--peptideProteinGrouping          Perform peptide protein grouping export.
 -rCoH,--readCometHistograms <arg>      Directory where Comet histograms files are placed.
 -repH,--reportHistogram                Report histograms to text files
 -rMqH,--readMaxQuantHistograms <arg>   Directory where MaxQuant histograms files are placed.
 -rP,--readParamFile <arg>              Name of file from which parameters should to read.
```

The comand line options in more detail:

```
-coD,--cometPsmDir <arg>           Comet psm root directory (required)

Directory where Comet pep.xml files are placed. All pep.xml files under this directory including sub-directories 
and matching the -coRE regular expression are considered for the analysis.
```

```
-coFDR,--cometFDR <arg>            FDR for filtering Comet PSMs before combination (required) (default value 0.03)

FDR for Comet search. If -fdrM is equal to "combined", the FDR is calculated for all PSM's of all groups that
pass the local lFDR threshold. If -fdrM is equal to "separate", the FDR is calculated seperately for all
groups defined by the -groupN option. The former method, which is the one discussed in the paper, is usually less 
conservative than the "separate" method.
```

```
-coRE,--cometPsmRegex <arg>        Regular expression of Comet psm files (e.g. \.xml$) (required)

Regular expression that defines the Comet pep.xml files used in the analysis.
 ```

```
-exclP,--excludeProts <arg>        Regular expression of proteins excluded from analysis. If not set no proteins are excluded.

This option can be used to exclude proteins (e.g. contaminant proteins) from the analysis. If the regular expression matches a 
protein string, the protein is excluded.
```

```
-fdrM,--fdrControlMethod <arg>     Method to control pFDR: global or groupwise (default global).

The FDR can be controlled in two ways: separate and combined. In both cases, the local FDR (lFDR) is calculated in the same way. The 
lFDR threshold is then adjusted to yield the target FDR, which is estimated for all groups together (combined) or for each 
group seperately (separate). Estimating the FDR seperately for each group usually yields less but more accurate PSMs.
```

```
-fH,--forceHistograms              Histograms are imported even if enough PSMs are available.

This option is only used if the -readH option is set. If both the -fH and -readH options are set, the histograms are imported 
even if enough PSMs for a specific histogram are available to calculate lFDR (as defined by the -minPH option).
```

```
-groupF,--groupProteinFile <arg>  Tab file with protein group assignments which will override assignment by groupRE
 
If -groupRE cannot capture all protein-group assignments, then the protein group assignments can be explicitely stated in the tab file provided in this option. A line in the file has the format 'protein\tgroup', whereas protein is a string that uniquely identifies a protein fasta header (can be a substring) and groups is one of the groups defined in the -groupN option.
```
 
```
-groupM,--groupingMethod <arg>    Method for PSM grouping: fasta, modif, famo (fasta&modif) or none (default none).
 
Protein grouping can be omitted ('none') or based on the fasta headers ('fasta') or PTMs of the peptides ('modif') or both fasta headers and modif ('famo').
If -groupM == 'famo', then each group defined by the -groupRE and -groupN options is split into modified and nonmodified peptides. -groupRE and -groupN options are only considered if -groupM is set to 'fasta' or 'famo'. 
```

```
-groupN,--groupNames <arg>        Comma separated list of names of sequence groups in fasta file (e.g. prot,lncRNA,TE ). 
                                  Will be used as prefixes for output files.
 
If -groupM is set to 'fasta', then this option defines the names assigned to the groups defined by the -groupRE option. The 
comma separated list of names must contain n+1 names if the -groupRE list contains n regular expressions. Proteins matching 
the first regular expression will be assigned to the first group, all remaining proteins matching the second regular 
expression will be assigned to the second group and so on. All remaining proteins will be assigned to the last group. 
```

```
-groupRE,--groupRegEx <arg>       Comma separated list of regular expression defining sequence groups of fasta headers (e.g.
                                  "sp\||tr\|ENSP00","ENST00","SINE_|LINE_|LTR_|DNA_|Retroposon_" ). Will be used as prefixes
                                  for output files.
                                  
If -groupM is set to 'fasta', then this option defines the regular expressions to assign proteins to groups. Proteins matching 
the first regular expression will be assigned to the first group, all remaining proteins matching the second regular 
expression will be assigned to the second group and so on. All remaining proteins will be assigned to the last group. 
```

```
-h,--help                          Help option for command line help
```

```
-maxDC,--maxDeltaCn <arg>          Maximal Comet DeltaCn in histogram (default value 2500)

Maximal value of Comet deltaCn score in histogram.
```

```
-maxDS,--maxDeltaScore <arg>          Maximal MaxQuant DeltaScore in histogram (default value 300)

Maximal MaxQuant DeltaScore in histogram.
```

```
-maxL,--maxLength <arg>            Maximal length of peptide (default value: 25)

Maximal length for peptides considered by NewAnce.
```

```
-maxPEP,--maxPEP <arg>          Maximal MaxQuant PEP value in histogram (default value 0.1)

Maximal MaxQuant PEP value in histogram.
```

```
-maxR,--maxRank <arg>              Maximal rank of peptide in list of spectrum matches (rank 1 = best) (default value: 1)

Maximal rank for PSMs considered in NewAnce.
```

```
-maxScore,--maxScore <arg>          Maximal MaxQuant Score in histogram (default value 500)

Maximal MaxQuant PEP value in histogram.
```
 
```
-maxSP,--maxSpScore <arg>          Maximal Comet SpScore in histogram (default value 1)

Maximal value for Comet SPScore in histogram.
```

```
-maxXC,--maxXCorr <arg>            Maximal Comet XCorr in histogram (default value 5)

Maximal value for Comet XCorr in histogram.
```

```
-maxZ,--maxCharge <arg>            Maximal charge of PSM (default value: 5)

Maximal charge for PSMs considered by NewAnce.
```

```
-minDC,--minDeltaCn <arg>          Minimal Comet DeltaCn in histogram (default value 0)

Minimal value for Comet DeltaCn in histogram.
```

```
-minDS,--minDeltaScore <arg>          Minimal MaxQuant DeltaScore in histogram (default value 0)

Minimal MaxQuant DeltaScore in histogram.
```

```
-minL,--minLength <arg>            Minimal length of peptide (default value: 8)

Minimal length for peptides considered by NewAnce.
```

```
-minPEP,--minPEP <arg>          Minimal MaxQuant PEP value in histogram (default value 0)

Minimal MaxQuant PEP value in histogram.
```

```
-minPH,--minPsm4Histo <arg>        Minimal number of psms to calculate local FDR in histogram (default value: 100000).

Minimal number of psms to calculate local lFDR in a histogram. If less data points are available, the lFDR estimate is 
considered unreliable and a precalculated default histogram is used instead for charge 1-3. For charge states higher than 3 
the default histogram for charge 3 is used. The default histograms were obtained with immunopeptidomics MS/MS spectra in an 
OrbiTrap at a resolution of 15'000. For other MS/MS acquistion methods, the user can  import his/her own default histograms 
(see -readH option).
```

```
-minScore,--minScore <arg>          Minimal MaxQuant Score in histogram (default value 0)

Minimal MaxQuant PEP value in histogram.
```

```
-minSP,--minSpScore <arg>          Minimal Comet SpScore in histogram (default value 0)

Minimal value for Comet SpScore in histogram.
```

```
-minXC,--minXCorr <arg>            Minimal Comet XCorr in histogram (default value 0)

Minimal value for Comet XCorr in histogram.
```

```
-minZ,--minCharge <arg>            Minimal charge of PSM (default value: 1)

Minimal charge for PSMs considered by NewAnce.
```

```
-mod,--modifications <arg>         Comma separated list of peptide modifications used in search (e.g. Cysteinyl:C3H5NO2S,Oxidation:O)
 
The user can define a comma separated list of molecules (format name:formula e.g. Cysteinyl:C3H5NO2S, Oxidation:O, 
Carbamidomethyl:C2H3NO, Phospho:HO3P). If the mass of the molecule matches the mass of the modification, the name defined in 
for the modification is used in the output tables. This is useful if the user wants to change the name of a modification or if 
the modification used in the search is not in UniMod (default are UniMod names).
```

```
-mqD,--maxquantPsmDir <arg>        MaxQuant psm root directory. If not provided only Comet is used.

Directory under which MaxQuant msms.txt and peptides.txt files are placed. All msms.txt files under this directory are 
considered for the analysis, i.e. there can be several msms.txt files. NewAnce also requires the peptides.txt files, which 
have to be placed in the same directory as the corresponding msms.txt file.
```

```
-mqFDR,--maxquantFDR <arg>            FDR for filtering MaxQuant PSMs before combination (required) (default value 0.03)

FDR for MaxQuant search. If -fdrM is equal to "combined", the FDR is calculated for all PSM's of all groups that
pass the local lFDR threshold. If -fdrM is equal to "groupwise", the FDR is calculated seperately for all the
groups defined by the -groupN option. The former method, which is the one discussed in the paper, is usually less 
conservative than the "separate" method.
```

```
-nrDCB,--nrDeltaCnBins <arg>       Number of Comet DeltaCn bins in histogram (default value 40)
 
Number of Comet DeltaCn bins in histogram. The higher this number, the finer the grid cells for FDR calculation and the more 
data is needed. 
```

```
-nrDSB,--nrDeltaScoreBins <arg>       Number of MaxQuant DeltaScore bins in histogram (default value 40)
 
Number of MaxQuant DeltaScore bins in histogram. The higher this number, the finer the grid cells for FDR calculation and the more 
data is needed. 
```

```
-nrPEPB,--nrPEPBins <arg>       Number of MaxQuant PEP bins in histogram (default value 40)
 
Number of MaxQuant PEP bins in histogram. The higher this number, the finer the grid cells for FDR calculation and the more 
data is needed. 
```

```
-nrScoreB,--nrScoreBins <arg>       Number of MaxQuant Score bins in histogram (default value 40)
 
Number of MaxQuant Score bins in histogram. The higher this number, the finer the grid cells for FDR calculation and the more 
data is needed. 
```

```
-nrSPB,--nrSpScoreBins <arg>       Number of Comet SpScore bins in histogram (default value 40)
 
Number of Comet SpScore bins in histogram. The higher this number, the finer the grid cells for FDR calculation and the more 
data is needed. 
```

```
-nrTh,--nrThreads <arg>            Number of threads used by NewAnce (default value: nr of available processors - 2)
 
Number of threads assigned for the execution of NewAnce. 
```

```
-nrXCB,--nrXCorrBins <arg>         Number of Comet XCorr bins in histogram (default value 40)

Number of Comet XCorr bins in histogram. The higher this number, the finer the grid cells for FDR calculation and the more 
data is needed. 
```

```
-outD,--outputDir <arg>            Output directory for results (required)

Directors where output files are written to. If output files with the same name already exist, they will be overwritten.
```

```
-outT,--outputTag <arg>            Tag inserted into output file names after prefix.
 
PSM ouput file have the fomat group_tag_NewAncePSMs.txt. The tag can be specified here. 
```

```
-ppG,--peptideProteinGrouping      Perform peptide protein grouping export.

This option has to be set if the protein grouping has to be performed. NewAnce implements a greedy protein grouping algorithm 
described in the publication, where it groups proteins, that share peptides together.
```

```
-rCoH,--readCometHistograms <arg>      Directory where Comet histograms files are placed.

Instead of being calculated with the available Comet PSMs, precaluclated score histograms are imported. This is useful in case there are not 
enough data to build them. NewAnce provides default histograms for high resolution OrbiTrap immunopeptidomics data. Otherwise 
a user can create histograms by analysing a larger dataset (20+ runs) and exporting the histograms (CometHistogramCalculator, see below). When importing histograms the histogram settings -maxDC,-maxSP,-maxXC,-maxZ,-minDC,-minSP,-minXC,-minZ,-nrDCB,-nrSPB,-nrXCB are overwritten.
```

```
-repH,--reportHistogram            Report histograms to text files

Comet score histograms are written to text files. These histograms can be imported with the -readH option and used as score 
histograms in case there is not enough data to build them.
```

```
-rMqH,--readMaxQuantHistograms <arg>      Directory where MaxQuant histograms files are placed.

Instead of being calculated with the available MaxQuant PSMs, precaluclated score histograms are imported. This is useful in case there are not 
enough data to build them. NewAnce provides default histograms for high resolution OrbiTrap immunopeptidomics data. Otherwise 
a user can create histograms by analysing a larger dataset (20+ runs) and exporting the histograms (MaxQuantHistogramCalculator). When importing histograms the histogram settings -maxDC,-maxSP,-maxXC,-maxZ,-minDC,-minSP,-minXC,-minZ,-nrDCB,-nrSPB,-nrXCB are overwritten.
```

```
-rP,--readParamFile <arg>          Name of file from which parameters should to read.

NewAnce parameters are read from the specified parameter file.
```

```
-seFa,--searchFastaFile <arg>      Fasta file that was used for the search (required for protein grouping export and 
                                   annotation of variants)

Fasta file that was used for Comet search. This is required if the protein grouping option -ppG is set. Also, for variants 
defined in the PEFF format in the fasta header (e.g \VariantSimple=(141|A|rs75062661_0) ) Comet does not report the comment 
field (rs75062661_0) in the pepxml file. If the -seFa option is provided, NewAnce parses the search fasta file and retrieves 
the variant annotation and adds it to the results files.
```

```
-smD,--smoothDegree <arg>          Degree of smoothing (0: no smoothing, n: n x smoothing) (default value 1)

Comet score histograms can be smoothed, where NewAnce implements a simple nearest neighbour smoothing (i.e. a cell's psm count 
is replaced by the average psm counts of the cell and its 6 nearest neighbours). The smoothing algorithm can be run n times 
(n=0: no smoothing). The higher n the higher the degree of smoothing.
```

```
-spRE,--spectrumFilter <arg>       If this option is set, only spectrum ids that match this regexp are used.  If not set no filtering is performed.

In case NewAnce should be limited to certain spectra, these spectra can be defined with this regular expression. The spectrum 
format used in NewAnce is: spectrum_raw_file_name.scan_number.scan_number.charge 
```

```
-upFa,--uniProtFastaFile <arg>     Fasta file with coding or canonical proteins (e.g. UniProt fasta file)

All peptides matching to a protein sequence in this file are considered as peptides of protein coding genes. Exceptions can be 
defined with the -noncP option.
```

```
-v,--version                       Version of NewAnce software

Show NewAnce version.
```

```
-wAll,--writeFullPSMExport         If flag is set, all Comet and MaxQuant PSMs are written to a tab file.

All PSMs identified by Comet and MaxQuant are added to a tab file ending with _CometPSMs.txt or _MaxQuantPSMs.txt. The local FDR values and 
whether the PSM passes the FDR thresholf defined in the -coFDR or -mqFDR options are added for each PSM.
```

```
-wP,--write2ParamFile <arg>        Filename where parameters should be written to.

Write NewAnce parameters to specified parameter file.
```


## Running NewAnce

The NewAnce command line options can be obtained by typing

```
java -jar NewAnce-1.7.3-SNAPSHOT.jar -h
```

The NewAnce version can be obtained by typing

```
java -jar NewAnce-1.7.3-SNAPSHOT.jar -v
```

Running NewAnce with 12GB memory

```
java -Xmx12G -jar NewAnce-1.4.0-SNAPSHOT.jar -coD 0D5P/lncRNA/Comet -coRE .*pep.xml$ -mqD 0D5P/lncRNA/MaxQuant -coFDR 0.03 -outD 0D5P/lncRNA/NewAnce -outT 0D5P -groupM fasta -groupRE "sp\||tr\||ENSP00,ENST00" -groupN prot,lncRNA,ere -upFa SeqDBs/human_proteome.fasta -maxR 1 -minZ 1 -maxZ 3 -minL 8 -maxL 15
```


## Output format

NewAnce writes the PSMs for the protein coding and non-canonical groups to tab delimited text files . The file names are groupName_outputTag_NewAncePSMs.txt. The columnnames in these files are:

* Spectrum : spectrum id in the format raw_file_name.scan_number.scan_number.charge . Leading 0's are removed from the scan_number tag. 
* ScanNr : scan number of MS/MS acquisition 
* Charge : precursor ion charge 
* RT : retention time of MS/MS acquisition in minutes
* NeutralMass : neutral mass of precursor ion
* Peptide : amino acid sequence of matching peptide with PTM annotation
* Sequence : bare amino acid sequence of matching peptide
* PeptideMass : neutral mass of matching peptide
* ModifName : name of PTM (if present, NA otherwise). Comma separated if several PTMs are present.
* ModifPosition : position of PTM (if present, NA otherwise). Comma separated if several PTMs are present.
* ModifMass : mass of PTM (if present, NA otherwise). Comma separated if several PTMs are present.   
* ModifAA : amino acid of PTM (if present, NA otherwise). NT (CT) for n(c)-terminal modification. Comma separated if several PTMs are present. 
* Proteins : NewAnce protein ids in format: \[protein1,protein2,protein3]     
* IsVariant : true if peptide contains sequence variant annotated in fasta header of Comet or MaxQuant searches. false otherwise.   
* VariantStart : start position of variant in WT sequence (if present, NA otherwise). Comma separated if several variants are present.
* VariantEnd : end position of variant in WT sequence (if present, NA otherwise). Comma separated if several variants are present.
* VariantSeq : variant sequence (if present, NA otherwise). Comma separated if several variants are present.  
* WTSeq : WT sequence(if present, NA otherwise). 10 AA are added to N- and C-term.
* IsDecoy : true if peptide originates from reversed protein sequence. false otherwise.
* Rank : rank of PSM in Comet search
* Group : PSM group assigned by NewAnce (one of the groups defined in -groupN option)
* XCorr : Comet xcorr score  
* DeltaCn : Comet deltacn score  
* SpScore : Comet spscore score    
* Expect : Comet expect score
* Massdiff : Comet mass difference between neutral precursor and peptide masses (atomic units)
* Tot_num_ions : Comet total number of theoretical fragment ions considered for the match  
* Num_matched_ions : Comet number of theoretical fragment ions that produced a match to experimental fragment masses
* Comet.lFDR : local FDR for Comet PSMs calculated by NewAnce 
* Comet.passFDR : true if PSMs passed the Comet FDR filter, false otherwise
* Mass.Error[ppm] : MaxQuant relative mass error between neutral precursor and peptide masses (ppm)  
* Score : MaxQuant Score
* Delta.score : MaxQuant Delta score
* PEP : MaxQuant PEP value
* Localization.prob : MaxQuant Localization probability of PTM (NA if no PTM is present)
* MaxQuant.lFDR : local FDR for MaxQuant PSMs calculated by NewAnce 
* MaxQuant.passFDR : true if PSMs passed the MaxQuant FDR filter, false otherwise



## NewAnce with Comet only
The group-wise FDR calculation is the same as in NewAnce, but only based on Comet Results and no consensus filtering with MaxQuant results is performed. 

```
The command line options can be obtained by

java -Xmx12G -cp NewAnce-1.7.3-SNAPSHOT.jar newance.scripts.CometGroupFDRFilter -h

usage: newance.scripts.CometGroupFDRFilter
 -alpha,--combineHistoWeight <arg>   Histograms are alpha*data_histo + (1-alpha)*prior_histo. 0 <0 alpha <= 1. Only used if rCoH option is set.
 -coD,--cometPsmDir <arg>            Comet psm root directory (required)
 -coFDR,--cometFDR <arg>             FDR for filtering Comet PSMs before combination (default value 0.03)
 -coRE,--cometPsmRegex <arg>         Regular expression of Comet psm files (e.g. \.xml$) (required)
 -d,--debug                          Debug option
 -exclP,--excludeProts <arg>         Regular expression of proteins excluded from analysis. If not set no proteins are excluded.
 -fdrM,--fdrControlMethod <arg>      Method to control pFDR: combined or separate (default combined).
 -fH,--forceHistograms               Histograms are imported even if enough PSMs are available.
 -groupF,--groupProteinFile <arg>    Tab file with protein group assignments which will override assignment by groupRE
 -groupM,--groupingMethod <arg>      Method for PSM grouping: fasta, modif, famo (fasta&modif) or none (default none).
 -groupN,--groupNames <arg>          Comma separated list of names of sequence groups in fasta file (e.g. prot,lncRNA,TE ). Will be used as prefixes
                                     for output files.
 -groupRE,--groupRegEx <arg>         Comma separated list of regular expression defining sequence groups of fasta headers (e.g.
                                     "sp\||tr\|ENSP00","ENST00","SINE_|LINE_|LTR_|DNA_|Retroposon_" ). Will be used as prefixes for output files.
 -h,--help                           Help option for command line help
 -maxDC,--maxDeltaCn <arg>           Maximal Comet DeltaCn in histogram (default value 2500)
 -maxL,--maxLength <arg>             Maximal length of peptide (default value: 25)
 -maxR,--maxRank <arg>               Maximal rank of peptide in list of spectrum matches (rank 1 = best) (default value: 1)
 -maxSP,--maxSpScore <arg>           Maximal Comet SpScore in histogram (default value 1)
 -maxXC,--maxXCorr <arg>             Maximal Comet XCorr in histogram (default value 5)
 -maxZ,--maxCharge <arg>             Maximal charge of PSM (default value: 5)
 -minDC,--minDeltaCn <arg>           Minimal Comet DeltaCn in histogram (default value 0)
 -minL,--minLength <arg>             Minimal length of peptide (default value: 8)
 -minPH,--minPsm4Histo <arg>         Minimal number of psms to calculate local FDR in histogram (default value: 100000).
 -minSP,--minSpScore <arg>           Minimal Comet SpScore in histogram (default value 0)
 -minXC,--minXCorr <arg>             Minimal Comet XCorr in histogram (default value 0)
 -minZ,--minCharge <arg>             Minimal charge of PSM (default value: 1)
 -mod,--modifications <arg>          Comma separated list of peptide modifications used in search (e.g. Cysteinyl:C3H5NO2S,Oxidation:O)
 -nrDCB,--nrDeltaCnBins <arg>        Number of Comet DeltaCn bins in histogram (default value 40)
 -nrSPB,--nrSpScoreBins <arg>        Number of Comet SpScore bins in histogram (default value 40)
 -nrTh,--nrThreads <arg>             Number of threads used by NewAnce (default value: nr of available processors - 2)
 -nrXCB,--nrXCorrBins <arg>          Number of Comet XCorr bins in histogram (default value 40)
 -outD,--outputDir <arg>             Output directory for results (required)
 -outT,--outputTag <arg>             Tag inserted into output file names after prefix.
 -ppG,--peptideProteinGrouping       Perform peptide protein grouping export.
 -rCoH,--readCometHistograms <arg>   Directory where Comet histograms files are placed.
 -repH,--reportHistogram             Report histograms to text files
 -rP,--readParamFile <arg>           Name of file from which parameters should to read.
 -seFa,--searchFastaFile <arg>       Fasta file that was used for the search (required for protein grouping export and annotation of variants in the
                                     Comet results)
 -smD,--smoothDegree <arg>           Degree of smoothing (0: no smoothing, n: n x smoothing) (default value 1)
 -spRE,--spectrumFilter <arg>        If this option is set, only spectrum ids that match this regexp are used.  If not set no filtering is performed.
 -upFa,--uniProtFastaFile <arg>      Fasta file with coding or canonical proteins (e.g. UniProt fasta file)
 -v,--version                        Version of NewAnce software
 -wAll,--writeFullCometExport        If flag is set, all Comet PSMs are written to a tab file.
 -wP,--write2ParamFile               This option is set if parameters should be written to file.

```

```
Running NewAnce/CometGroupFDRFilter with 12GB memory


java -Xmx12G -cp NewAnce-1.7.3-SNAPSHOT.jar newance.scripts.CometGroupFDRFilter -coD 0D5P/lncRNA/Comet -coRE .HLAIp.*pep\.xml$ -fdrM combined -coFDR 0.03 -outD 0D5P/lncRNA/NewAnce -outT 0D5P  -groupM famo -groupRE "sp\||tr\||ENSP00,ENST00" -groupN "prot,lncRNA,ere" -upFa SeqDBs/human_proteome.fasta -seFa SeqDBs/OD5P_peff.fasta -smD 5 -maxR 1 -minZ 1 -maxZ 3 -minL 8 -maxL 15 -wAll -rCoH histos/Comet/Type_I'
```



## Tests

### 1 Testing Comet pep.xml reader

The CometXMLReaderTest tests how NewAnce reads the Comet pep.xml file and outputs the PSMs in tab format to standard output. This can be used to check whether the Comet pep.xml file is read correctly.

Option for CometXMLReaderTest additional to NewAnce options:

```
-maxP,--maxDisplayedPsms <arg>          Maximal number of psms written to standard output

Maximal number of psms written to standard output. If option is not set, all PSMs in the pep.xml files are written.
```

Printing CometXMLReaderTest options:

```
java -jar -cp NewAnce-1.4.0-SNAPSHOT.jar newance.testers.CometXMLReaderTest -h
```

Running CometXMLReaderTest:

```
java -Xmx12G -cp NewAnce-1.4.0-SNAPSHOT.jar newance.testers.CometXMLReaderTest -coD 0D5P/lncRNA/Comet -coRE .*pep.xml$ -maxR 1 -minZ 1 -maxZ 3 -minL 8 -maxL 15 -maxP 10
```

### 2 Testing MaxQuant msms.txt reader

The MaxQuantPsmReaderTest tests how NewAnce reads the MaxQuant msms.txt and peptides.txt files and outputs the PSMs in tab format to standard output. This can be used to check whether the MaxQuant files are read correctly.

Option for MaxQuantPsmReaderTest additional to NewAnce options:

```
-maxP,--maxDisplayedPsms <arg>          Maximal number of psms written to standard output

Maximal number of psms written to standard output. If option is not set, all PSMs in the pep.xml files are written.
```

Printing MaxQuantPsmReaderTest options:

```
java -cp NewAnce-1.4.0-SNAPSHOT.jar newance.testers.MaxQuantPsmReaderTest -h
```

Running MaxQuantPsmReaderTest:

```
java -Xmx12G -cp NewAnce-1.4.0-SNAPSHOT.jar newance.testers.MaxQuantPsmReaderTest -mqD 0D5P/lncRNA/MaxQuant -maxR 1 -minZ 1 -maxZ 3 -minL 8 -maxL 15 -maxP 10
```

### 3 Testing fasta file reader

The FastaReaderTest test reads a fasta file and the fasta entry in tab format to standard output. This can be used to check whether the fasta files are read correctly. This also outputs the protein id used in NewAnce (NewAnceID), which is the protein id matched by the regular expression in the NewAnce -protRE option.

Options for FastaReaderTest:

```
-maxP,--maxDisplayedPsms <arg>          Maximal number of psms written to standard output

Maximal number of fasta entries written to standard output. If option is not set, all entries in the fasta file are written.
```

```
-fa,--fastaFile <arg>                   Fasta file to be tested

This test reveals how NewAnce reads and interprets the given fasta file
```

Printing FastaReaderTest options:

```
java -cp NewAnce-1.4.0-SNAPSHOT.jar newance.testers.FastaReaderTest -h
```

Running FastaReaderTest:

```
java -Xmx12G -cp NewAnce-1.4.0-SNAPSHOT.jar newance.testers.FastaReaderTest -fa SeqDBs/uniprot_human_reviewed_proteome_18_12_2018.fasta -maxP 10
```

### 4 Testing specific peptides

To obtain more information about specific peptides, e.g. whether they were identified by MaxQuant or Comet before consensus 
and FDR filtering, you can run the CometMaxQuantScoreCombiner_PeptideTest class. It will print all PSMs related to these 
peptides to standard output.

Printing CometMaxQuantScoreCombiner_PeptideTest options:

```
java -cp NewAnce-1.4.0-SNAPSHOT.jar newance.testers.CometMaxQuantScoreCombiner_PeptideTest -h
```

Option for PeptideTest additional to NewAnce options:

```
-pept,--peptides <arg>          Comma separated list of peptides to be printed (e.g. [TPAPRPLGI,VIDYPPIAY,AQFRVTEA]).

Peptide sequences for which more information is requested.
```

Running PeptideTest for peptides TPAPRPLGI,VIDYPPIAY, and AQFRVTEA:

```
java -Xmx12G -cp NewAnce-1.4.0-SNAPSHOT.jar newance.testers.CometMaxQuantScoreCombiner_PeptideTest -coD 0D5P/lncRNA/Comet -coRE .*pep.xml$ -mqD 0D5P/lncRNA/MaxQuant -coFDR 0.03 -outD 0D5P/lncRNA/NewAnce -outT 0D5P -protRE sp\||tr\| -protG prot -noncG lncRNA -upFa SeqDBs/human_proteome.fasta -maxR 1 -minZ 1 -maxZ 3 -minL 8 -maxL 15 -pepts [TPAPRPLGI,VIDYPPIAY,AQFRVTEA]
```

## PDV export

The [PDV viewer](http://pdv.zhang-lab.org) is a MS/MS search result viewer. In order to export the peptide-spectrum matches to PDV format, the .mgf files containing the spectra and NewAnce result files are required.

Options for CreatePDVExport:

```
-mgf,--mgfDir <arg>                    Directory containing mgf files (required)

Directors where mgf files of spectra used for Comet/MaxQuant searches are found.
```

```
-naf,--newAnceResultFile <arg>         Result file from NewAnce analysis (required)

NewAnce result file containing peptide spectrum matches.
```

Printing CreatePDVExport options:

```
java -cp NewAnce-1.4.0-SNAPSHOT.jar newance.scripts.CreatePDVExport -h
```

Running CreatePDVExport:

```
java -Xmx12G -cp NewAnce-1.4.0-SNAPSHOT.jar newance.scripts.CreatePDVExport -mgf 0D5P/mgf -naf 0D5P/lncRNA/NewAnce/lncRNA_0D5P_NewAncePSMs.txt
```

## Extend MaxQuant features

The default NewAnce output files contain the features described above. MaxQuant msms.txt result files contain many very useful PSM features, which can be added to the NewAnce result files. This will create a new NewAnce result file with \_extend tag inserted at the end of the filename.  

Options for AddMaxQuantFeatures:

```
-mqD,--maxquantPsmDir <arg>                    MaxQuant psm root directory (required)

Directory under which MaxQuant msms.txt files are found.
```

```
-mqF,--maxquantFeatures <arg>                    Comma separated list [feature1,feature2,feature3] of MaxQuant features to be added to NewAnce result file.

Comma separated list in format [feature1,feature2,feature3] of MaxQuant features to be added to NewAnce result file.
```

```
-naf,--newAnceResultFile <arg>         Result file from NewAnce analysis (required)

NewAnce result file containing peptide spectrum matches.
```

Printing AddMaxQuantFeatures options:

```
java -cp NewAnce-1.4.0-SNAPSHOT.jar newance.scripts.AddMaxQuantFeatures -h
```

Running AddMaxQuantFeatures:

```
java -Xmx12G -cp NewAnce-1.4.0-SNAPSHOT.jar newance.scripts.AddMaxQuantFeatures -mqD 0D5P/lncRNA/MaxQuant -mqF "[Matches,Intensities,Mass Deviations [Da]]" -naf 0D5P/lncRNA/NewAnce/lncRNA_0D5P_NewAncePSMs.txt
```

## Export score histograms

This class can be used to calculate and export prior score histograms for NewAnce. It is implemented in a memory efficient 
way, so it can be run with many (100's) Comet .pep.xml files. The prior score histograms can then be imported by NewAnce with 
the -readH option in case NewAnce is run with only a few .pep.xml files, which do not contain sufficient PSMs to calculate the score histograms accurately. Histograms will be written to the output directory specified with the -outD option. There will be a histogram file for each charge x, which will be named outputTag_Zx.txt.

Printing CometHistogramCalculator options:

```
java -cp NewAnce-1.4.0-SNAPSHOT.jar newance.scripts.CometHistogramCalculator -h
```

Running CometHistogramCalculator:

```
java -Xmx12G -cp NewAnce-1.4.0-SNAPSHOT.jar newance.scripts.CometHistogramCalculator -coD 0D5P/lncRNA/Comet -coRE .*pep.xml$ -outD 0D5P/lncRNA/NewAnce/histos -outT prior_histo -spRE .* -maxR 1 -minZ 1 -maxZ 3 -minL 8 -maxL 15 -nrTh 5
```

## Spectronaut DIA export

This class can be used to export the NewAnce PSMs to the Spectronaut format for spectrum library generation and DIA analysis 
(version 1.7.1 or higher).  Information for precursor ion intensities and scan event numbers can be obtained either from 
MaxQuant result files or from mgf files. 

```
-mqD,--maxquantPsmDir <arg>            MaxQuant psm root directory

Root directory under which MaxQuant msms.txt file(s) is (are) found. All msms.txt file in under this 
root directory are considered. The 'Scan Event' and 'MS1 Intensity' intensity values in the Spectronaut file
are obtained from the values in the msms.txt files.
```

```
-mgf,--mgfDir <arg>            Directory containing mgf files

Root directory under which mgf files are found. The 'Scan Event' and 'MS1 Intensity' intensity values in the Spectronaut file
are calculated from the mgf files. If -mqD option is set too, the MaxQuant msms.txt file(s) are used and the mgf files are 
ignored. If neither -mgf nor -mqD options are set, the 'Scan Event' and 'MS1 Intensity' intensity values are set to NA.
```

```
-naf,--newAnceResultFile <arg>         Result file from NewAnce analysis (required)

NewAnce result file containing peptide spectrum matches (PSMs). Each PSM in this file is exported to Spectronaut format. 
"_Spectronaut" tage will be added to the end of the NewAnce file name.
```

Printing CreateSpectronautExport options:

```
java -cp NewAnce-1.7.1-SNAPSHOT.jar newance.scripts.CreateSpectronautExport -h
```

Running CreateSpectronautExport:

```
java -Xmx12G -cp NewAnce-1.7.1-SNAPSHOT.jar newance.scripts.CreateSpectronautExport -mqD 0D5P/lncRNA/MaxQuant -naf 0D5P/lncRNA/NewAnce/lncRNA_0D5P_NewAncePSMs.txt 
```

## Match peptides to fasta file

This class can be used to match peptide sequences from a tab file to protein sequences in a fasta file (matches treat I/L as equal). 
This can be useful if you want to check whether ceratin peptides are part of UniProt for example.


Printing AnnotatePeptides options:

```
java -cp NewAnce-1.4.0-SNAPSHOT.jar newance.scripts.AnnotatePeptides -h
```

Running AnnotatePeptides:

```
java -Xmx12G -cp NewAnce-1.4.0-SNAPSHOT.jar newance.scripts.AnnotatePeptides -fa SeqDBs/human_proteome.fasta -col PeptideSequence -p 0D5P/tmp/peptides.txt 
```

