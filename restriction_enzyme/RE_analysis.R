# Run this script from inside a project folder with .Rprofile sourcing RE_profile.R #

source("../../my_Rprofile.R")


#### 0 Initialization ####

if(grepl("_REL606$", getwd())) {
  
  bam_folder <- paste0(sub("_REL606$", "", getwd()), "/data/bam_files/")
  RData_folder <- paste0(sub("_REL606$", "", getwd()), "/data/RData_files/")  # uses results for the normal non-REL606 run
  
  dir.create("data/RData_files/", recursive = TRUE)  # to save some REL606 results here
  
  bam_names <- gsub(".bam", "", grep("\\.bam$", dir(bam_folder), value = TRUE))
  if(length(bam_names) == 0) {
    bam_names <- gsub(".RData", "", grep(".RData", dir(RData_folder), value = TRUE))
    bam_names <- bam_names[!grepl("^cut_uncut_fragment_lengths", bam_names)]
  }
  
  main_genome_name = "REL606"
  main_genome <- read_genome_fasta("../../reference_genomes/Ecoli_spike-in/")
  names(main_genome) <- main_genome_name
  
  bad_regions <- data.frame(chr = character(0), start = numeric(0), end = numeric(0))
} else {
  
  bam_folder <- "data/bam_files/"
  RData_folder <- "data/RData_files/"
  dir.create(RData_folder)
  
  bam_names <- gsub(".bam", "", grep("\\.bam$", dir(bam_folder), value = TRUE))
  if(length(bam_names) == 0) {
    bam_names <- gsub(".RData", "", grep(".RData", dir(RData_folder), value = TRUE))
    bam_names <- bam_names[!grepl("^cut_uncut_fragment_lengths", bam_names)]
  }
  
  main_genome_name = "cerevisiae"
  main_genome <- read_genome_fasta("../../reference_genomes/S288C_reference_genome_R64-1-1_20110203/")
  
  bad_regions <- read.table("../../external_data/bad_regions.tab", header = TRUE)
}

norm_genome_name <- "pombe"
norm_genome <- read_genome_fasta("../../reference_genomes/GSE28839/")

genome_names <- c("main" = main_genome_name, "norm" = norm_genome_name)
genome <- c(main_genome, norm_genome)

chr_sizes <- width(genome)
names(chr_sizes) <- names(genome)

chr_name_list <- list("main" = names(main_genome), "norm" = names(norm_genome))

RE_info <- read.table("../RE_info.txt", stringsAsFactors = FALSE)

if(length(bam_names) != sum(grepl("([a-zA-Z0-9\\(\\)]_1$|[a-zA-Z0-9\\(\\)]_X$)", bam_names))) {
  stop("all bam file names need to end in _X or _1")
}

# For sample naming rules, see the README file
sample_to_RE_table <- data.frame(RE = sapply(bam_names, function(name) paste(rownames(RE_info)[sapply(paste0("(^", rownames(RE_info), "|_", rownames(RE_info), ")"), grepl, name)], collapse = "+")))
sample_to_RE_table$RE <- as.character(sample_to_RE_table$RE)

if(any(sample_to_RE_table$RE=="")) {
  stop("all bam file need to contain at least one restriction enzyme name")
}

analysis_folder <- "analysis_results/"
dir.create(analysis_folder)


#### 1.1 Load bam and save as RData  ####

bam_duplicates_flag <- NA  # NA: keep duplicates, FALSE: drop duplicates, TRUE: only duplicates (BWA does not set this flag anyway, so this approach can not find duplicates)

bam_names_todo <- character(length = 0)
for(name in bam_names) {
  if(isTRUE(bam_duplicates_flag)) {
    RData_name <- paste0(name, "_bam_duplicates_only")
  } else if(isTRUE(!bam_duplicates_flag)) {
    RData_name <- paste0(name, "_no_bam_duplicates")
  } else {
    RData_name <- name  # this is what we will assume the rest of the script
  }
  if(!file.exists(paste0(RData_folder, RData_name, ".RData"))) {
    bam_names_todo <- rbind(bam_names_todo, name)
  }
}

temp_function <- function(bam_name) save_PE_bam_as_chr_df_list(bam_name, bam_folder = bam_folder, RData_folder = RData_folder, bam_duplicates_flag = bam_duplicates_flag)
run_function_parallel(bam_names_todo, temp_function, core_factor = 1/4, single_core_to_debug = FALSE)

#### 1.2 Save sample data without duplicates (in addition to the normal sample data) ####

add_sample_analysis_without_duplicates <- FALSE

if(add_sample_analysis_without_duplicates) {
  bam_names_todo <- character(length = 0)
  for(name in bam_names) {
    if(!file.exists(paste0(RData_folder, name, "_no_duplicates.RData"))) {
      bam_names_todo <- rbind(bam_names_todo, name)
    }
  }
  
  temp_function <- function(bam_name) remove_duplicate_fragments(bam_name, RData_folder = RData_folder)
  run_function_parallel(bam_names_todo, temp_function, core_factor = 1, single_core_to_debug = FALSE)
  
  # Alternative: read in bam and save the duplicates as RData
  # for(name in bam_names) {
  #   save_bam_as_chr_df_list(paste0(name, "_duplicate"), bam_folder, duplicate=TRUE)
  # }
  # any_duplicate_comparison(data_names, genome, chr_name_list, sample_to_RE_table, RE_info,
  #                          bad_regions, multi_position_window=multi_position_window_directly_at_cut_site, output_file=paste0(plot_folder, "stats_duplicates.txt"))  
}


#### 1.3 Plot fragment length histograms ####

temp_list <- plot_PE_fragment_length_histograms(bam_names, chr_name_list, plot_folder = analysis_folder, RData_folder = RData_folder)
bam_names_low_coverage_X_1 <- temp_list[[1]]
mean_fragment_length_df <- temp_list[[2]]
rm(temp_list)

if(length(bam_names_low_coverage_X_1)>0) {
  print("WARNING: The following files have low coverage in main or norm genome, consider removing them from the analysis with the code below, if you find unexpected low coverage.")
  print(bam_names_low_coverage_X_1)
}

if(TRUE) {
  data_names <- bam_names  # do not remove low coverage files by default
} else {
  bam_names_low_coverage <- unique(sub("(_X$|_1$)", "", bam_names_low_coverage_X_1))
  data_names <- setdiff(bam_names, paste0(bam_names_low_coverage, c(rep("_X", length(bam_names_low_coverage)), rep("_1", length(bam_names_low_coverage)))))
}


#### 1.4 Calc counts near average site ####

large_window_limit <- 200
small_window_limits <- c(-10, 40)
close_distances_for_distr_near_cut_sites <- list("dist_ign_site" = 200, "dist_ign_half" = 300)

if(add_sample_analysis_without_duplicates) {
  sample_to_RE_table_normal <- sample_to_RE_table
  sample_to_RE_table_no_duplicates <- sample_to_RE_table
  rownames(sample_to_RE_table_no_duplicates) <- paste0(rownames(sample_to_RE_table_no_duplicates), "_no_duplicates")
  sample_to_RE_table_all <- rbind(sample_to_RE_table, sample_to_RE_table_no_duplicates)
  
  sample_to_RE_table <- sample_to_RE_table_all
  
  data_names_normal <- data_names
  data_names <- sort(c(data_names_normal, paste0(data_names_normal, "_no_duplicates")))
}

temp_function <- function(data_names) calc_cuts_of_average_site_and_plot_genomic_cuts(data_names, genome, chr_name_list, sample_to_RE_table,
                                                                                      RE_info, bad_regions, large_window_limit, small_window_limits, plot_folder = analysis_folder,
                                                                                      close_distances = close_distances_for_distr_near_cut_sites, RData_folder = RData_folder, max_length = 500)
# cuts_near_average_site_list_test <- temp_function(data_names[c(1, 2)])
# cuts_near_average_site_list <- cuts_near_average_site_list_test

cuts_near_average_site_list <- run_function_parallel(data_names, temp_function, core_factor = 1/2, single_core_to_debug = FALSE)


#### 1.5 Plot read distribution around average site ####

count_window_df <- calc_count_limits_and_plot_average_site(cuts_near_average_site_list, chr_name_list, genome_names, large_window_limit, small_window_limits, plot_folder = analysis_folder)

# hack to shorten names in table
count_window_df$name <- sapply(as.character(count_window_df$name), function(x){strsplit(x, "\\.")[[1]][1]})

write.table(count_window_df, file = paste0(analysis_folder, "count_window_limits_and_resection_lengths.tsv"), sep="\t", row.names = FALSE)

count_window_df

save(data_names, cuts_near_average_site_list, count_window_df, RE_info, file = paste0("counts_near_average_site.RData"))


#### 2.1 Count cuts in / outside window ####

use_both_strands <- TRUE
window_limit_factor <- 1
max_length <- 500

if(!use_both_strands){
  window_name <- paste0("window_limit_times_", window_limit_factor, "_not_both_strands")
  stop("Implementation not complete for use_both_strands == FALSE")
} else {
  window_name <- paste0("window_limit_times_", window_limit_factor, "_max_length_", max_length)
}

window_folder <- paste0(analysis_folder, window_name, "/")
dir.create(window_folder)


# Single threaded test:
# count_list_test <- count_fragments_at_cut_sites(data_names, genome, chr_name_list, sample_to_RE_table, RE_info,
#                                                 subset(count_window_df, close_dist_str == paste0(close_distances_for_distr_near_cut_sites, collapse = "_")),
#                                                 max_length, bad_regions, window_limit_factor = window_limit_factor, use_both_strands = use_both_strands)
# count_list <- count_list_test

temp_function <- function(data_names) count_fragments_at_cut_sites(data_names, genome, chr_name_list, sample_to_RE_table, RE_info,
                                                                   subset(count_window_df, close_dist_str == paste0(close_distances_for_distr_near_cut_sites, collapse = "_")),
                                                                   max_length, bad_regions, window_limit_factor = window_limit_factor, use_both_strands = use_both_strands)

parallel_count_list <- run_function_parallel(data_names, temp_function, core_factor = 1/2, single_core_to_debug = FALSE)

# parallel run needs moving around the samples inside the different list, which is done here
str(parallel_count_list, max.level=2)
count_list <- list()
for(list_name in unique(names(parallel_count_list))) {  # restructuring of output list (move all samples into one list type)
  count_list[[list_name]] <- do.call(c, parallel_count_list[names(parallel_count_list) == list_name])
  names(count_list[[list_name]]) <- sub(paste0(list_name, "."), "",  names(count_list[[list_name]]))
}
str(count_list, max.level=2)
rm(parallel_count_list)

list2env(count_list, .GlobalEnv)
rm(count_list)

# filter samples with nan in cut counts
data_names_nan_in_cut_counts <- character(0)
for(name in data_names) {
  if(any(is.nan(do.call(c, cut_counts_nf_df_list[[name]][, 3:9])))) {
    data_names_nan_in_cut_counts <- c(data_names_nan_in_cut_counts, name)
  }
}
data_names_nan_in_cut_counts <- unique(data_names_nan_in_cut_counts)
print(data_names_nan_in_cut_counts)
cut_counts_nf_df_list <- cut_counts_nf_df_list[setdiff(data_names, data_names_nan_in_cut_counts)]
bg_counts_df_list <- bg_counts_df_list[setdiff(data_names, data_names_nan_in_cut_counts)]
data_names <- setdiff(data_names, data_names_nan_in_cut_counts)

# new 2021-09-19
for(name in data_names) {
  for(this_enzyme in rownames(RE_info)) {  # assume name contains "EcoRI-norm" and this_enzyme = "EcoRI", then
    if(grepl(paste0(this_enzyme,"-norm"), name)) {
      # remove EcoRI sites in the main genome
      cut_counts_nf_df_list[[name]] = cut_counts_nf_df_list[[name]][with(cut_counts_nf_df_list[[name]], enzyme != this_enzyme | chr %in% chr_name_list[["norm"]]), ]
    }
  }
}

save(data_names, cut_counts_nf_df_list, bg_counts_df_list, fragment_counts_list, count_window_list, window_name, window_limit_factor, window_folder, RE_info,
     file = paste0("counts_", window_name, ".RData"))


#### 2.2 Plot cut and uncut counts agains distance to next neighbour ####

plot_counts_vs_nn_distance(cut_counts_nf_df_list, chr_name_list, chr_sizes, genome_names, plot_folder = window_folder, outlier_quantile=0.99)
#plot_counts_vs_nn_distance_alternative(cut_counts_nf_df_list, chr_name_list, chr_sizes, genome_names, plot_folder = window_folder, outlier_quantile=0.99)


#### 2.3 Filter out close sites ####

close_distances <- list("dist_ign_site" = 200, "dist_ign_half" = 300)

filter_name <- paste0(window_name, "_close_distances_", paste0(close_distances, collapse="_"))
filter_folder <- paste0(window_folder, "close_distances_", paste0(close_distances, collapse="_"), "/")
dir.create(filter_folder)

cut_counts_df_list <- set_close_site_cut_counts_NA(cut_counts_nf_df_list, genome, sample_to_RE_table, RE_info, bad_regions, close_distances)

# lapply(cut_counts_df_list, function(df) aggregate(df[, 4:ncol(df)], by = list(df[, "enzyme"]), mean, na.rm = TRUE))

# new 2021-09-19
for(name in data_names) {
  for(this_enzyme in rownames(RE_info)) {  # assume name contains "EcoRI-norm" and this_enzyme = "EcoRI", then
    if(grepl(paste0(this_enzyme,"-norm"), name)) {
      # remove non-EcoRI sites in the norm genome i.e. in S. pombe, only after the close site filtering
      cut_counts_df_list[[name]] = cut_counts_df_list[[name]][with(cut_counts_df_list[[name]], (chr %in% chr_name_list[["main"]]) | enzyme == this_enzyme), ]  # this needs to be improved to work with more than one "-norm" enzyme if necessary
      cut_counts_nf_df_list[[name]] = cut_counts_nf_df_list[[name]][with(cut_counts_nf_df_list[[name]], (chr %in% chr_name_list[["main"]]) | enzyme == this_enzyme), ]
    }
  }
}


#### 2.4 Cut statistics and accs/occs from mean counts ####

cut_statistics_fc <- print_cut_count_statistics(cut_counts_df_list, bg_counts_df_list, fragment_counts_list, count_window_list, chr_name_list, output_file=paste0(window_folder, "cut_count_statistics_fc.txt"))

print_occs_from_mean_counts(cut_statistics_fc, normalization_method="norm_cuts", output_file=paste0(filter_folder, "occs_from_mean_counts_fc_normed_by_norm_cuts.txt"))
print_occs_from_mean_counts(cut_statistics_fc, normalization_method="norm_reads", output_file=paste0(filter_folder, "occs_from_mean_counts_fc_normed_by_norm_reads.txt"))
print_occs_from_mean_counts(cut_statistics_fc, normalization_method="main_reads", output_file=paste0(filter_folder, "occs_from_mean_counts_fc_normed_by_main_reads.txt"))


#### 2.5 Plot count histograms and comparisons ####

plot_count_histograms(cut_counts_df_list, chr_name_list, genome_names, plot_folder = filter_folder, outlier_quantile=0.99)
#plot_count_histograms_alternative(cut_counts_df_list, chr_name_list, genome_names, plot_folder = filter_folder, outlier_quantile=0.99)
plot_count_comparisons(cut_counts_df_list, chr_name_list, genome_names, plot_folder = filter_folder)
plot_cut_uncut_comparisons(cut_counts_df_list, chr_name_list, genome_names, plot_folder = filter_folder)


plot_count_histograms(cut_counts_nf_df_list, chr_name_list, genome_names, plot_folder = filter_folder, outlier_quantile=0.99, name_add = "_not_filtered")
plot_count_comparisons(cut_counts_nf_df_list, chr_name_list, genome_names, plot_folder = filter_folder, name_add = "_not_filtered")
plot_cut_uncut_comparisons(cut_counts_nf_df_list, chr_name_list, genome_names, plot_folder = filter_folder, name_add = "_not_filtered")


#### 2.6 Plot cut and uncut fragments lengths ####

temp_function <- function(name) plot_cut_uncut_length_distributions(name, cut_counts_df_list, count_window_list, RE_info, max_length = max_length, plot_folder = filter_folder)
run_function_parallel(data_names, temp_function, core_factor = 1/2, single_core_to_debug = FALSE)

temp_function <- function(name) plot_cut_uncut_length_distributions(name, cut_counts_df_list, count_window_list, RE_info, max_length = max_length, plot_folder = filter_folder, name_add = "_not_filtered")
run_function_parallel(data_names, temp_function, core_factor = 1/2, single_core_to_debug = FALSE)


#### 2.7 Whats happening at close sites ####

if(FALSE) {
  # uncut counts of one-side close sites / non-close sites
  sapply(cut_counts_df_list, function(df) with(df, mean((uncut_plus + uncut_minus)[is.na(starts_plus) | is.na(ends_plus)]) /
                                                 mean((uncut_plus + uncut_minus)[!(is.na(starts_plus) | is.na(ends_plus))])))
  
  # uncut counts of two-side close sites / non-close sites
  sapply(cut_counts_df_list, function(df) with(df, mean((uncut_plus + uncut_minus)[is.na(starts_plus) & is.na(ends_plus)]) /
                                                 mean((uncut_plus + uncut_minus)[!(is.na(starts_plus) | is.na(ends_plus))])))
  
  # start cut counts of close sites / non-close sites
  sapply(data_names, function(name) with(cut_counts_nf_df_list[[name]], mean((starts_plus + starts_minus)[is.na(cut_counts_df_list[[name]]$starts_plus)]) /
                                           mean((starts_plus + starts_minus)[!is.na(cut_counts_df_list[[name]]$starts_plus)])))
  
  # end cut counts of close sites / non-close sites
  sapply(data_names, function(name) with(cut_counts_nf_df_list[[name]], mean((ends_plus + ends_minus)[is.na(cut_counts_df_list[[name]]$ends_plus)]) /
                                           mean((ends_plus + ends_minus)[!is.na(cut_counts_df_list[[name]]$ends_plus)])))
  
  # start coverage of close sites / non-close sites
  sapply(data_names, function(name) with(cut_counts_nf_df_list[[name]], mean((starts_plus + starts_minus + uncut_plus + uncut_minus)[is.na(cut_counts_df_list[[name]]$starts_plus)]) /
                                           mean((starts_plus + starts_minus + uncut_plus + uncut_minus)[!is.na(cut_counts_df_list[[name]]$starts_plus)])))
  
  # end coverage of close sites / non-close sites
  sapply(data_names, function(name) with(cut_counts_nf_df_list[[name]], mean((ends_plus + ends_minus + uncut_plus + uncut_minus)[is.na(cut_counts_df_list[[name]]$ends_plus)]) /
                                           mean((ends_plus + ends_minus + uncut_plus + uncut_minus)[!is.na(cut_counts_df_list[[name]]$ends_plus)])))
}


#### 3.1 Calc accessibilities ####

background_method <- "corrected"  # options: ("untreated", "corrected"), applies only to cut-allcut method, cut-uncut is always corrected!

min_coverage <- 40  # only average over sites with at least this effective coverage (also used for optimizing the cut uncut correction)

acc_max <- 1.5  # cap value for accessibilities greater than 1
acc_min <- -0.25  # lower limit for plotting (values could still be lower, but are bound by -shear_prob/(1-shear_prob))

background_name <- paste0(filter_name, "_background_", background_method)
background_folder <- paste0(filter_folder, "background_", background_method, "/")
dir.create(background_folder)

save.image()

#### 3.1.1 Calc and plot deviation from calibration samples ####

if(grepl("RE_calibration", getwd())) {  # run this for new calibration in RE_calibration folder
  
  REs_to_calibrate = c("BamHI")
  
  calibration_folder <- paste0(background_folder, "calibration/")
  dir.create(calibration_folder)
  
  uncut_loss_correction_values <- seq(1.0, 2.0, 0.01)  # for finding the best uncut_loss_correction in RE_calibration

  calc_accs_and_deviation <- function(uncut_loss_correction) {
    temp_list <- calc_accessibilities(cut_counts_df_list, bg_counts_df_list, count_window_list, background_method, chr_name_list, acc_max, uncut_loss_corrections = rep(uncut_loss_correction, 2))
    accs_df_list <- temp_list[[1]]
    uncut_loss_correction_list <- lapply(c("combined", REs_to_calibrate),
                                         function(enzyme) calc_deviation_from_test_runs(accs_df_list, calibration_folder,
                                                                                        name_add = paste0("_", as.character(1000*uncut_loss_correction)),
                                                                                        min_coverage = min_coverage, enzyme = enzyme))
    deviation_df <- bind_rows(lapply(uncut_loss_correction_list, "[[", "deviation_df"))
    sample_df <- bind_rows(lapply(uncut_loss_correction_list, "[[", "sample_df"))
    deviation_df$uncut_loss_correction <- uncut_loss_correction
    sample_df$uncut_loss_correction <- uncut_loss_correction
    
    return(list("deviation_df" = deviation_df, "sample_df" = sample_df))
  }
  
  uncut_loss_correction_list <- run_function_parallel(uncut_loss_correction_values, calc_accs_and_deviation, core_factor = 1/2, single_core_to_debug = FALSE, unlist_after = FALSE)

  sample_df <- bind_rows(lapply(uncut_loss_correction_list, "[[", "sample_df"))
  deviation_df <- bind_rows(lapply(uncut_loss_correction_list, "[[", "deviation_df"))
  
  deviation_df$rel_error_2 <- deviation_df$cut_uncut_all_2 / deviation_df$cut_uncut_all_1
  deviation_df$rel_error_4 <- deviation_df$cut_uncut_all_4 / deviation_df$cut_uncut_all_3
  
  save(deviation_df, sample_df, file = paste0(calibration_folder, "uncut_loss_correction_dfs.RData"))
  
  uncut_loss_corrections_per_enzyme <- list()
  for(this_error_method in c("sqrt_sample_mean_squ_diff_site_mean", "sqrt_sample_mean_site_mean_squ_diff")) {
    uncut_loss_corrections_per_enzyme[[this_error_method]] <- list()
    for(this_enzyme in unique(deviation_df$enzyme)) {
      error_method_df <- subset(deviation_df, error_method == this_error_method & enzyme == this_enzyme)
      min_2 <- with(error_method_df, uncut_loss_correction[min(rel_error_2) == rel_error_2])
      min_4 <- with(error_method_df, uncut_loss_correction[min(rel_error_4) == rel_error_4])
      
      error_method_df_combined <- combine_cols_for_plot(error_method_df, c("rel_error", "occ_method"), c("rel_error_2" = "occ_method_2", "rel_error_4" = "occ_method_4"))
      
      ggplot(error_method_df_combined, aes_string(x="uncut_loss_correction", y="rel_error", color="occ_method")) + geom_point() + ylim(c(0, 1)) +
        ggtitle(paste0("y-axis: Ratio of corrected error over uncorrected error\nMinima: method 2 at ", min_2, ", method 4 at ", min_4))
      ggsave(file = paste0(calibration_folder, "optimization_curve_", this_enzyme, "_", this_error_method, "_min_coverage_", min_coverage, ".pdf"))
      
      uncut_loss_corrections_per_enzyme[[this_error_method]][[this_enzyme]] <- c(min_2, min_4)
    }
  }
  
  temp_df <- subset(deviation_df, error_method == "sqrt_sample_mean_squ_diff_site_mean" & enzyme != "combined")
  temp_df$RE <- factor(temp_df$enzyme, levels = REs_to_calibrate)
  ggplot(temp_df, aes(x=uncut_loss_correction, y=rel_error_4, color=RE)) + geom_line() + ylab("relative fit error") + xlab("uncut correction factor") +
    theme_paper + ggplot_colours + ylim(c(0, 1))
  ggsave(paste0(calibration_folder, "calibration_curve.pdf"), height = 5, width = 7, units = "cm")
  
  save(uncut_loss_corrections_per_enzyme, file="uncut_loss_corrections_per_enzyme.RData")
  
  this_error_method <- "sqrt_sample_mean_squ_diff_site_mean"
  good_correction_sample_df <- bind_rows(lapply(c("combined", REs_to_calibrate),
                                                function(this_enzyme) subset(sample_df, enzyme == this_enzyme & uncut_loss_correction == uncut_loss_corrections_per_enzyme[[this_error_method]][[this_enzyme]])))
  
  good_correction_sample_df$fit_method <- ifelse(good_correction_sample_df$enzyme == "combined", "combined", "individual")

  good_correction_sample_df$enzyme <- sapply(strsplit(good_correction_sample_df$file_name, split = "_"), function(strs) REs_to_calibrate[REs_to_calibrate %in% strs])
  
  temp <- subset(good_correction_sample_df, fit_method == "combined")[, c("enzyme", "calibration_samples_accs", "all_mean", "sd_all_mean")]
  colnames(temp) <- c("enzyme", "sample_occ", "mean_occ", "sd_occ")
  temp$method <- "X-1 method"
  calibration_results_df <- temp
  
  temp <- subset(good_correction_sample_df, fit_method == "combined")[, c("enzyme", "calibration_samples_accs", "cut_uncut_all_3", "sd_cut_uncut_all_3")]
  colnames(temp) <- c("enzyme", "sample_occ", "mean_occ", "sd_occ")
  temp$method <- "cut-uncut unfitted"
  calibration_results_df <- rbind(calibration_results_df, temp)
  
  temp <- subset(good_correction_sample_df, fit_method == "combined")[, c("enzyme", "calibration_samples_accs", "cut_uncut_all_4", "sd_cut_uncut_all_4")]
  colnames(temp) <- c("enzyme", "sample_occ", "mean_occ", "sd_occ")
  temp$method <- "cut-uncut combined fit"
  calibration_results_df <- rbind(calibration_results_df, temp)
  
  temp <- subset(good_correction_sample_df, fit_method == "individual")[, c("enzyme", "calibration_samples_accs", "cut_uncut_all_4", "sd_cut_uncut_all_4")]
  colnames(temp) <- c("enzyme", "sample_occ", "mean_occ", "sd_occ")
  temp$method <- "cut-uncut individual fit"
  calibration_results_df <- rbind(calibration_results_df, temp)
  
  calibration_results_df$sample_occ <- 1-calibration_results_df$sample_occ
  calibration_results_df$mean_occ <- 1-calibration_results_df$mean_occ
  
  if(length(REs_to_calibrate) > 1) {
    optional_facet_wrap <- facet_wrap(~enzyme)
  } else {
    optional_facet_wrap <- list()
  }
  
  ggplot(calibration_results_df, aes(x=sample_occ, y=mean_occ, color = method, shape = method)) + geom_point(size=1) +
    xlim(c(0, 1)) + ylim(c(-0.05, 1.05)) + optional_facet_wrap + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey")
  ggsave(paste0(calibration_folder, "mean_occ.pdf"), width=10, height=3)
  
  ggplot(calibration_results_df, aes(x=sample_occ, y=sd_occ, color = method, shape = method)) + geom_point(size=1) + geom_line(size=0.3) +
    xlim(c(0, 1)) + ylim(c(0, 0.2)) + optional_facet_wrap
  ggsave(paste0(calibration_folder, "sd_occ.pdf"), width=10, height=3)
  
  ggplot(calibration_results_df, aes(x=sample_occ, y=mean_occ, color = method, shape = method)) + geom_point(size=1) +
    optional_facet_wrap + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
    geom_errorbar(aes(ymin=mean_occ-sd_occ, ymax=mean_occ+sd_occ, width=.03)) + xlim(c(-0.05, 1.05)) + ylim(c(-0.2, 1.2))
  ggsave(paste0(calibration_folder, "errorbar_occ.pdf"), width=10, height=3)
  
  calibration_results_df$method[calibration_results_df$method=="X-1 method"] <- "cut-all cut"
  calibration_results_df$method <- factor(calibration_results_df$method, levels=c("cut-uncut unfitted", "cut-uncut individual fit", "cut-uncut combined fit", "cut-all cut"))
  
  if(length(REs_to_calibrate) > 1) {
    optional_facet_wrap <- facet_wrap(~method + enzyme, ncol = 3, nrow = 4)
  } else {
    optional_facet_wrap <- facet_wrap(~method, ncol = 3, nrow = 4)
  }
  
  ggplot(calibration_results_df, aes(x=sample_occ, y=mean_occ)) + optional_facet_wrap +
    diag_line + geom_errorbar(aes(ymin=mean_occ-sd_occ, ymax=mean_occ+sd_occ), size=0.25, width=0.05, alpha=0.7) + geom_point(size=0.5) + 
    ylab("measured mean absolute occupancy / %") + xlab("prepared fraction of uncut DNA molecules before shearing / %") + labs(color = "enzyme") + ggplot_colours +
    theme_paper + theme(legend.position="none") + 
    scale_y_continuous(labels = percent_no_sign, limit = c(-0.2, 1.1), breaks = seq(0, 1, 0.25)) + scale_x_continuous(labels = percent_no_sign, limit = c(-0.05, 1.05))
  
  ggsave(paste0(calibration_folder, "occ_errobars_facet_method.pdf"), units="cm", scale = 1, width=15, height=15)
}


#### 3.1.2 Calc accs for best uncut_loss_correction ####

if(TRUE) {  # load uncut_loss_corrections_per_enzyme from file
  load("../RE_calibration/uncut_loss_corrections_per_enzyme_OLD.RData")
} else {    # or set manually here (the two values per RE correspond to two different acc. calculation methods, with the second method being the one to be used by default. Replace with actual values!)
  uncut_loss_corrections_per_enzyme <- list()
  uncut_loss_corrections_per_enzyme[["sqrt_sample_mean_squ_diff_site_mean"]]$AluI = c(1, 1)
  uncut_loss_corrections_per_enzyme[["sqrt_sample_mean_squ_diff_site_mean"]]$BamHI = c(1, 1)
  uncut_loss_corrections_per_enzyme[["sqrt_sample_mean_squ_diff_site_mean"]]$HindIII = c(1, 1)
  uncut_loss_corrections_per_enzyme[["sqrt_sample_mean_squ_diff_site_mean"]]$KpnI = c(1, 1)
  uncut_loss_corrections_per_enzyme[["sqrt_sample_mean_squ_diff_site_mean"]]$combined = c(1, 1)
}

temp_list <- calc_accessibilities(cut_counts_df_list, bg_counts_df_list, count_window_list, background_method, chr_name_list, acc_max,
                                  uncut_loss_corrections_per_enzyme = uncut_loss_corrections_per_enzyme[["sqrt_sample_mean_squ_diff_site_mean"]])
list2env(temp_list, .GlobalEnv)
rm(temp_list)


#### 3.2 Calc accessibility site means and sds ####

high_cov_fracs <- sapply(accs_df_list, function(accs_df) sum(accs_df$eff_coverage >= min_coverage, na.rm = TRUE) / sum(!is.na(accs_df$eff_coverage)))

if(any(high_cov_fracs < 0.8)) {
  warning("Some high_cov_fracs < 0.9")
  print(high_cov_fracs)
  #if(any(high_cov_fracs < 0.25)) stop("Some high_cov_fracs < 0.25")
}

occ_columns <- which("starts_plus" == colnames(accs_df_list[[1]])):ncol(accs_df_list[[1]])

acc_site_means <- data.frame(sapply(occ_columns, function(col) lapply(accs_df_list, function(accs_df) mean(accs_df[accs_df$eff_coverage >= min_coverage, col], na.rm = TRUE))))
acc_site_sds <- data.frame(sapply(occ_columns, function(col) lapply(accs_df_list, function(accs_df) sd(accs_df[accs_df$eff_coverage >= min_coverage, col], na.rm = TRUE))))
colnames(acc_site_means) <- colnames(accs_df_list[[1]][occ_columns])
colnames(acc_site_sds) <- colnames(accs_df_list[[1]][occ_columns])

cat(kable(acc_site_means, digits = 3), file=paste0(background_folder, "acc_site_means.txt"), sep="\n")
cat(kable(acc_site_sds, digits = 3), file=paste0(background_folder, "acc_site_sds.txt"), sep="\n")

if(FALSE & grepl("RE_calibration", getwd())) {  # Test site average with exp(mean(log x), for all sites (low coverage sites are not ignored)
  log_average <- function(mus) {
    log_mus <- log(mus)
    min_log_mu <- min(log_mus[is.finite(log_mus)])
    log_mus[is.infinite(log_mus)] <- min_log_mu
    exp(mean(log_mus, na.rm = TRUE))
  }
  mus_av_list <- lapply(mus_df_list, function(mus_df) lapply(mus_df[, -c(1, 2, 3)], log_average))
  acc_site_logmeans <- data.frame(t(sapply(names(mus_av_list), function(name) sapply(names(mus_av_list[[name]]), function(type) calc_accs_from_mus(mus_av_list[[name]][[type]], shear_prob_per_bp_list[[name]]*(length(count_window_list[[paste0(name,"_X")]][[1]])))))))
  cat(kable(acc_site_logmeans, digits = 3), file=paste0(background_folder, "acc_site_logmeans.txt"), sep="\n")
}


#### 3.3 Print simple summaries ####

index_X_1 = which(names(acc_site_means) %in% "all_mean")
index_cut_uncut <- which(grepl("cut_uncut_all", names(acc_site_means)))
cat(kable(acc_site_means[, c(index_X_1, index_cut_uncut)], digits = 3), file=paste0(background_folder, "acc_site_means_simple.txt"), sep="\n")
cat(kable(acc_site_sds[, c(index_X_1, index_cut_uncut)], digits = 3), file=paste0(background_folder, "acc_site_sds_simple.txt"), sep="\n")

acc_means_with_sd <- acc_site_means[, c("cut_uncut_all_3", "cut_uncut_all_4")]
colnames(acc_means_with_sd) <- c("mean_acc_uncorr", "mean_acc_correc")
acc_means_with_sd <- cbind(acc_means_with_sd, acc_site_sds[, c("cut_uncut_all_3", "cut_uncut_all_4")])
colnames(acc_means_with_sd)[c(3,4)] <- c("sd_acc_uncorr", "sd_acc_correc")
acc_means_with_sd<- acc_means_with_sd[, c(1,3,2,4)]
cat(kable(acc_means_with_sd, digits = 3), file=paste0(background_folder, "acc_cut_uncut_means_with_sd.txt"), sep="\n")


#### 3.4 Plot accessibility histograms and comparisons ####

plot_accessibility_histograms(accs_df_list, plot_folder = background_folder, acc_limits = c(acc_min, acc_max), min_coverage = 0)
plot_accessibility_histograms(accs_df_list, plot_folder = background_folder, acc_limits = c(acc_min, acc_max), min_coverage = min_coverage)

plot_accessibility_comparisons(accs_df_list, plot_folder = background_folder, acc_limits = c(acc_min, acc_max), min_coverage = 0)
plot_accessibility_comparisons(accs_df_list, plot_folder = background_folder, acc_limits = c(acc_min, acc_max), min_coverage = min_coverage)


#### 3.5 Calc occupancies ####

occs_df_list <- calc_occs_df_list(accs_df_list, min_coverage = 0)
save(occs_df_list, file=paste0(background_folder, "occs_df_list.RData"))  # save with low coverage sites

occs_df_list <- calc_occs_df_list(accs_df_list, min_coverage = min_coverage)  # continue plotting without low coverage sites

lapply(occs_df_list, function(df) aggregate(df[, 4:ncol(df)], by = list(df[, 2]), mean, na.rm = TRUE))  # for checking different enzyme results in the same sample


#### 3.6 Plot genomic occupancies ####

plot_genomic_occs(occs_df_list, plot_folder = background_folder, occ_limits = c(1-acc_max, 1-acc_min))


if(main_genome_name == "cerevisiae") {
  
  
  #### 4.1 Plot single genes/promoters ####
  
  single_genes <- read.table("../../external_data/single_genes.tab", header=TRUE, stringsAsFactors=FALSE)
  
  plot_single_genes(single_genes, occs_df_list, plot_folder = background_folder, occ_limits = c(1-acc_max, 1-acc_min))
  
  
  #### 4.2 Align occupancies ####
  
  alignment_window <- 1000
  
  genes_to_align <- read.table("../../external_data/Xu_2009_ORF-Ts_V64.bed", header=FALSE, sep="\t", quote="")[, c(1, 2, 3, 6)]
  chrs_main <- chr_name_list$main
  levels(genes_to_align[, 1]) <- chrs_main
  
  # occs_aligned is very memory inefficient (consists mostly of NAs), needs a lot of RAM
  occs_aligned <- align_occupancies(occs_df_list, genome, genes_to_align, chrs_main, alignment_window)
  
  save(occs_aligned, alignment_window, file=paste0("aligned_occupancies_", background_name, ".RData"))
  
  
  #### 4.3 Plot average gene ####
  
  avg_occs_list <- plot_promoter_regions(occs_aligned, alignment_window, plot_folder = background_folder, occ_limits = c(1-acc_max, 1-acc_min))
  list2env(avg_occs_list, .GlobalEnv)
  save(avg_gene_occs, avg_gene_occs_smoothed, window, file=paste0("average_gene_occupancies_", background_name, ".RData"))
  rm(avg_occs_list)
  
  rm(occs_aligned)  # don't save these again, it takes too long
  
}


#### 5 Save results ####

save(list = setdiff(ls(), lsf.str()), file = paste0(background_name, ".RData"))  # save all data excluding functions
backup_folder <- "used_scripts/"
dir.create(backup_folder)
file.copy("../RE_analysis.R", paste0(backup_folder, "RE_analysis_", background_name, ".R"), overwrite = TRUE)  # save scripts used to calculate this data
file.copy("../RE_functions.R", paste0(backup_folder, "RE_functions_", background_name, ".R"), overwrite = TRUE)

save(list = ls(), file = ".RData")
