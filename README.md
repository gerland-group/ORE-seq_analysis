# ORE-seq Analysis Script and Functions

## Requirements

Tested with R versions v3.4.3 and v3.6.3. Needs packages detailed in `restriction_enzyme/RE_Rprofile.R`

## How to use

1. Map fragments with bowtie2 using the combined S. cerevisiae and S. pombe reference genome (see `reference_genome/ScerAndSpomWithMT.fsa`: with chromosomes named as follows: chr01 – chr16 for the sixteen S. cerevisiae and chrI, chrII, chrIII for the three S. pombe chromosomes). Use alignment parameters: `-X 500 --no-discordant --no-mixed --no-unal`.
3. Remove multiply mapped reads using `samtools view -hf 0x2`
4. Index BAM file using `samtools index`
5. Download this repository
6. Install R & packages detailed in `restriction_enzyme/RE_Rprofile.R`
7. Make new folder `<Example>` in folder 'restriction_enzyme' for your analysis
8. Make new folders `data/bam` within `<Example>`
9. Put bam and bai files into `<Example>/data/bam`. Name files according to rules below.
10. Start R and set working directory to `<Example>`
11. `source("../RE_analysis.R")` or run `RE_analysis.R` step by step in RStudio from within `<Example>`
  
## Sample naming rules

Bam files within `data/bam` need to follow these naming conventions:
* The script needs bam files for both samples (X% cut and 100% cut) with identical file name except the ending: Samples with one RE digest end with `_X.bam`, while samples with second RE digest end with `_1.bam`
* If there was no second digest and only cut-uncut analysis is wanted, the `_1.bam` file can be a copy / hard link of the `_X.bam` file and the cut-all cut results should then be ignored.
* File names must contain the RE name of the enzyme present in the sample after a "_" sign e.g. `<Strain>_BamHI-HF_<RE units>_X.bam`, where the information in `<Strain>` and `<RE units>` is not used by the script and `<RE units>` could be ommitted.
* If the spike-in (if present) used a different RE, then add `<spike-in RE>-norm`, e.g. `<Strain>_AluI_EcoRI-norm_<RE units>_X.bam`
* Usable RE names can be checked and added in the `RE_info.txt` file
* Multiple REs can be used on the main genome (not the spike-in): `<Strain>_BamHI-HF_<RE units>_KpnI_<RE units>_X.bam`, which will be analysed accordingly.
  * To only get results of one RE if others where present in the same sample, set parentheses to ignore REs: `<Strain>_BamHI-HF_<RE units>_(KpnI_400)_X.bam`. This will be analysed similarly to `<Strain>_BamHI-HF_<RE units>_X`, with the following difference:
  * The sites of the ignored RE (and their neighbourhoods) will still be excluded when calculating the background
  * The sites of the ignored RE might exclude sites of the "main" enzyme when they are close to each other
* For calibration samples to be used for fitting the uncut correction factor, add the cut percentage with `X_pct_cut` as in this example: `<Strain>_AluI_10_pct_cut.bam`
      
## How to use other REs

Add the required information to `RE_info.txt`, see `RE_info_README.txt`.
  
## How to use other genomes

Other genomes than S. cerevisiae (and S. pombe for spike-in) are not supported by default, as there are unfortunately several references to the chromosome names within the code. If these are treated properly, the script should run for other genomes as well.
  
## How to fit uncut correction factors

Run section `3.1.1 Calc and plot deviation from calibration samples` in the script that is skipped by default.

## ToDo
  
* Upload reference genome in reference_genome/ScerAndSpomWithMT.fsa
* Upload code files: restriction_enzyme/RE_Rprofile.R, RE_info.txt and all others
