library('MAUDE')
# library('dplyr')
# library('tidyr')
library('tidyverse')
library('survcomp')

source('common.R')

unequal_truth = '/oak/stanford/groups/pritch/users/jake_harold/simulation/results/truth/unequal_bins_9_4_1000_100_0.10_3_20000_10000_3_410000_equalmixture_200.00_4.tsv'
filename = unequal_truth

# parse_unequal_bins_fname(unequal_truth)

# tmp = s_counts[!complete.cases(s_counts), ]


# bin_sizes = fname_to_bins(filename)

# filenames = Sys.glob('results/truth/*.tsv')
# filename = 'results/truth/vary_coverage_4_4_1000_100_0.10_3_20000_10000_1_205000_difficult_200.00_4.tsv'
# filenames = filenames[grepl('vary_coverage_[0-9]+_[0-9]+_', filenames)]
# filename = 'results/truth/nimble_vary_coverage_4_4_1000_100_0.10_3_20000_10000_1_205000_difficult_200.00_4.rds'

filename = '/oak/stanford/groups/pritch/users/jake_harold/simulation/results/truth/vary_coverage_2_bins_10_4_1000_100_0.10_3_20000_10000_1_1025000_equalmixture_200.00_4.tsv'

filename = '/oak/stanford/groups/pritch/users/jake_harold/simulation/results/maude/vary_coverage_2_bins_10_4_1000_100_0.10_3_20000_10000_1_1025000_equalmixture_200.00_4.csv'

filename = '/oak/stanford/groups/pritch/users/jake_harold/simulation/results/truth/vary_coverage_3_4_1000_100_0.10_3_20000_10000_1_1025000_difficult_200.00_4.tsv'

filename = '/oak/stanford/groups/pritch/users/jake_harold/simulation/results/truth/unimodal_moi_9_4_1000_100_0.10_3_20000_10000_1_7_50000_equalmixture_200.00_4.tsv'

# true_nimble = sub('truth', 'nimble', filename)
# true_nimble = sub('tsv', 'rds', true_nimble)
# fname = filename
# true_nimble = readRDS(filename)
# true_nimble = readRDS(true_nimble)

# filename = filenames[1]

args <- commandArgs(trailingOnly = TRUE)

filename = args[1]
out_file = args[2]


counts = read.table(filename, sep = '\t', header = TRUE)
counts[is.na(counts)] = 0


sim_type = fname_to_simulation_type(filename)

# debugonce(fname_to_bins)
bin_sizes = fname_to_bins(filename)[[1]]

bin_labels <- c('A', 'B', 'C', 'D')
bin_names <- c('Low', 'Low_Mid', 'Mid_High', 'High')

truth = readRDS(paste0('results/truth/sim_', sub('.tsv', '.rds', basename(filename))))

truth$counts[which(is.na(truth$counts))] = 0

# mutate counts table
g_counts = gather(counts, key = 'sample', value = 'count', -sgRNA, -Gene) %>%
  mutate(donor = str_extract(sample, 'D[0-9]'),
         bin = gsub('D[0-9]_', '', sample)) %>%
  dplyr::select(-sample)

unique(g_counts$bin)

s_counts = spread(g_counts, bin, count) %>%
  mutate(negControl = grepl('non-targeting', Gene, ignore.case = TRUE)) %>%
  # dplyr::rename(notSorted = input) %>%
  rename_at(all_of(bin_names), ~bin_labels)

s_counts = s_counts[c('sgRNA', 'Gene', 'donor', bin_labels, 'negControl')]

truth_df = data.frame(sgRNA = 1:length(unique(s_counts$sgRNA)), notSorted = truth$p_guide)

s_counts = left_join(s_counts, truth_df, by = 'sgRNA')

# Make bins
if (length(bin_sizes) == 1) {
  bin_bounds <- makeBinModel(data.frame(Bin = bin_labels, fraction = bin_sizes[[1]])) %>%
    filter(Bin %in% bin_labels)
} else {
  bin_bounds = lapply(bin_sizes,
    function(b) {
      makeBinModel(data.frame(Bin = bin_labels, fraction = b)) %>%
        filter(Bin %in% bin_labels)
    })
}


if (length(bin_sizes) == 1) {
  all_bounds = lapply(unique(s_counts$donor), function(d) {
    mutate(bin_bounds, donor = d)
    }) %>%
    bind_rows()
} else {
  donor_names = unique(s_counts$donor)
  all_bounds = lapply(seq_along(donor_names), function(i) {
    d = donor_names[i]
    mutate(bin_bounds[[i]], donor = d)
    }) %>%
    bind_rows()
}

# all_bounds = lapply(unique(s_counts$donor), function(d) {
#   mutate(bin_bounds, donor = d)
# }) %>%
#   bind_rows()

# get guide-level stats
guideLevelStats <- findGuideHitsAllScreens(experiments = unique(s_counts["donor"]),
                                           countDataFrame = s_counts,
                                           binStats = all_bounds,
                                           unsortedBin = "notSorted",
                                           negativeControl = "negControl",
                                           sortBins = bin_labels)

#get element level stats
elementLevelStats <- getElementwiseStats(unique(guideLevelStats["donor"]),
                                         guideLevelStats,
                                         elementIDs = "Gene",
                                         negativeControl = "negControl")


maude_results = summarize(group_by(elementLevelStats, Gene),
  p_value = combine.test(p.value))
maude_results = mutate(maude_results, fdr = p.adjust(p_value, method = 'BH'))

# # Call things sig if they are sig in 2 screens
# maude_results <- elementLevelStats %>%
#   mutate(donor_sig = FDR < 0.05) %>%
#   group_by(Gene) %>%
#   summarise(num_sig = sum(donor_sig),
#             avg_z = mean(meanZ)) %>%
#   mutate(screen_sig = num_sig >= 2) %>%
#   mutate(screen_sig = case_when(
#     screen_sig == TRUE ~ 'Maude sig',
#     screen_sig == FALSE ~ 'Maude not sig'
#   ))

write_csv(maude_results, out_file)
