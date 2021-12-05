# ORE-seq Analysis Script and Functions

## Requirements

Tested with R versions v3.4.3 and v3.6.3. Needs packages detailed in `restriction_enzyme/RE_Rprofile.R`

## How to use

1. Map fragments with bowtie2 using the combined S. cerevisiae and S. pombe reference genome (see `reference_genome/ScerAndSpomWithMT.fsa`: with chromosomes named as follows: chr01 – chr16 for the sixteen S. cerevisiae and chrI, chrII, chrIII for the three S. pombe chromosomes). Use alignment parameters: `-X 500 --no-discordant --no-mixed --no-unal`.
3. Remove multiply mapped reads using `samtools view -hf 0x2`
4. Index BAM file using `samtools index` to get `sample_name.bai` index files
5. Download this repository
6. Install R & packages detailed in `restriction_enzyme/RE_Rprofile.R`
7. Make new folder `<Example>` in folder 'restriction_enzyme' for your analysis
8. Make new folders `data/bam` within `<Example>`
9. Put bam and bai files into `<Example>/data/bam`. Name files according to rules below.
10. Start R and set working directory to `<Example>`
11. `source("../RE_analysis.R")` or run `RE_analysis.R` step by step in RStudio from within `<Example>`
12. Find desired output files and plots (see section below)
  
## Sample naming rules

Bam files within `data/bam` need to follow these naming conventions:
* The script needs bam files for both samples (X% cut and 100% cut) with identical file name except the ending: Samples with one RE digest end with `_X.bam`, while samples with second RE digest end with `_1.bam`.
* If there was no second digest and only cut-uncut analysis is wanted, the `_1.bam` file still needs to be provided, but can be a copy or soft link of the `_X.bam` file and the cut-all cut results should then be ignored.
* File names must contain the RE name of the enzyme present in the sample after a "_" sign e.g. `<Strain>_BamHI-HF_<RE units>_X.bam`, where the information in `<Strain>` and `<RE units>` is not used by the script and `<RE units>` could be ommitted.
* If the spike-in (if present) used a different RE, then add `<spike-in RE>-norm`, e.g. `<Strain>_AluI_EcoRI-norm_<RE units>_X.bam`
* Usable RE names can be checked and added in the `RE_info.txt` file
* Multiple REs can be used on the main genome (not the spike-in): `<Strain>_BamHI-HF_<RE units>_KpnI_<RE units>_X.bam`, which will be analysed accordingly.
  * To only get results of one RE if others where present in the same sample, set parentheses to ignore REs: `<Strain>_BamHI-HF_<RE units>_(KpnI_400)_X.bam`. This will be analysed similarly to `<Strain>_BamHI-HF_<RE units>_X`, with the following difference:
  * The sites of the ignored RE (and their neighbourhoods) will still be excluded when calculating the background
  * The sites of the ignored RE might exclude sites of the "main" enzyme when they are close to each other
* For calibration samples to be used for fitting the uncut correction factor, add the cut percentage with `X_pct_cut` as in this example: `<Strain>_AluI_10_pct_cut.bam`

## Analysis output

* The script creates a folder structure beginning with the folder `analysis_results` which contains different subfolders for each type of plot and result files.
* Depending on the parameters chosen in the script the main results path is `analysis_results/window_limit_times_1_max_length_500/close_distances_200_300/background_corrected` (in the following called `MAIN`) with plot folders for certain intermediate results along the way.
* Genomic mean accessibilities are saved in `MAIN/acc_site_means_simple.txt` with all_mean = cut-all cut results, cut_uncut_all_3 = uncorrected cut-uncut result, cut_uncut_4 = corrected cut-uncut result
* Histograms of site accessibilities are saved in `MAIN/accessibility_histograms/` for plus/minus strand and starting/ending fragments as well as combined results (last column) with all_mean = cut-all cut results, cut_uncut_all_3 = uncorrected cut-uncut result, cut_uncut_4 = corrected cut-uncut result
* Individual site results are saved in an R dataframe in `MAIN/occs_df_list.RData` with columns chr, enzyme, pos, eff_coverage, eff_cuts, occ_X_1 (= cut-all cut), occ_cut_uncut_2 (= uncorrected cut-uncut) and occ_cut_uncut_4 (=corrected cut-uncut)
      
## How to use other REs

Add the required information to `RE_info.txt`, see `RE_info_README.txt`.
  
## How to use other genomes

Other genomes than S. cerevisiae (and S. pombe for spike-in) are not supported by default, as there are unfortunately some references to the chromosome names within the code. If these are treated properly, the script should run for other genomes as well.
  
## How to fit uncut correction factors

Run section `3.1.1 Calc and plot deviation from calibration samples` in the script that is skipped by default.

## ToDo
  
* Rename columns in output files
* Add more output explanations
* Improve readability and structure of calibration
