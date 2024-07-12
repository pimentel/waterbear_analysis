library('dplyr')
library('tidyr')
library('ggplot2')
library('cowplot')

fname = '/oak/stanford/groups/pritch/users/jake_harold/simulation/results/nimble/vary_coverage_4_1_1000_100_0.25_2_20000_10000_1_2050000_equalmixture_200.00_4.rds'

# fname = metadata_filter$filename[1]
fname = '/oak/stanford/groups/pritch/users/jake_harold/simulation/results/truth/vary_coverage_2_4_1000_100_0.50_3_50000_25000_3_410000_easy_200.00_4.rds'

maude_fname = '/oak/stanford/groups/pritch/users/jake_harold/simulation/results/maude/vary_coverage_3_4_1000_100_0.10_3_20000_10000_1_410000_difficult_200.00_4.csv'

fname = sub('maude', 'nimble', maude_fname)
fname = sub('csv', 'rds', fname)

argv = commandArgs(trailingOnly = TRUE)
fname = argv[1]

mageck_fname = sub('.rds$', '.tsv', basename(fname))

cbb <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

theme_set(theme_cowplot())

sim = readRDS(paste0('results/truth/nimble_', basename(fname)))
# sim = readRDS(paste0('results/coverage/truth/', basename(fname)))

counts = sim$data$x

# TODO: remove debug
# which(is.na(counts[1, 2373, ]))
# end TODO:

n_donors = dim(counts)[1]
n_guides = dim(counts)[2]
n_bins = dim(counts)[3]

bin_names = NULL
if (n_bins == 4) {
  bin_names = c('Low', 'Low_Mid', 'Mid_High', 'High')
}

donor_names = paste0('D', 1:n_donors)

tmp = lapply(1:n_donors,
  function(i) {
    tmp = lapply(1:n_bins,
      function(j) {
        counts[i, , j]
      })
    names(tmp) = paste0(donor_names[i], '_', bin_names)
    tmp
  })

df = data.frame(do.call('c', tmp))
df = dplyr::mutate(df,
  sgRNA = c(sim$const$guide_index, sim$const$nt_data_index),
  Gene = c(sim$const$guide_to_gene, rep('Non-Targeting', length(sim$const$nt_data_index))))
df = dplyr::select(df, sgRNA, Gene, everything())

dir.create('results/truth/', recursive = TRUE, showWarnings = FALSE)

write.table(df, paste0('results/truth/', mageck_fname),
  quote = FALSE, row.names = FALSE, sep = '\t')

# mageck test -k results/truth/unimodal_moi_2_4_1000_100_0.75_3_50000_25000_1_10_1000000_difficult_200.00_4.tsv \
#   -t D1_Low,D2_Low,D3_Low -c D1_High,D2_High,D3_High --sort-criteria pos -n \
#   mageck_tmp

# mageck_results = read.table('mageck_tmp.gene_summary.txt', sep = '\t', header = TRUE)

# fdr = pmin(mageck_results$pos.fdr, mageck_results$neg.fdr)

# mageck_results = dplyr::mutate(mageck_results, value = pmin(pos.fdr, neg.fdr))
# mageck_results = dplyr::rename(mageck_results, gene_id = id)

# tmp = inner_join(mageck_results, oracle_mapping, by = 'gene_id')

# res = compute_benchmark_statistics(tmp)

# mean(fdr < 0.10)



# compute_benchmark_statistics(mageck_results)

# fdr_table = get_fdr_table_nimble(metadata_filter[1, ])
