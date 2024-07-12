library('dplyr')
library('tidyr')
library('ggplot2')
library('cowplot')

cbb <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

theme_set(theme_cowplot())

strip_text = function(x) {
  y = sub('^.+\\[', '', x)
  sub('\\]', '', y)
}

compute_benchmark_statistics = function(fdr_table) {
  # assumed you have grouped on appropriate columns (at least the method)
  stopifnot(!is.null(fdr_table))

  fdr_table = dplyr::arrange(fdr_table, value)
  fdr_table = dplyr::mutate(fdr_table,
    fp = cumsum(!true_effect), tp = cumsum(true_effect),
    fn = sum(true_effect) - tp, rank = 1:length(true_effect),
    sensitivity = tp / (fn + tp), fdr = fp / rank)

  fdr_table
}

generate_oracle_mapping = function(guide_to_gene, p_effects, n_total_genes) {
  p_effects = as.numeric(p_effects)
  n_total_genes = as.numeric(n_total_genes)
  oracle_mapping = data.frame(sim_id = unique(guide_to_gene))
  oracle_mapping = mutate(oracle_mapping,
    gene_id = as.character(sim_id))
  oracle_mapping = mutate(oracle_mapping,
    true_effect = FALSE)
  oracle_mapping$true_effect[1:(p_effects * n_total_genes)] = TRUE
  oracle_mapping
}

process_mageck_results = function(metadata_row) {
  info = metadata_row
  fname = nimble_fname_to_mageck(info$filename)
  fname = paste0('results/mageck/', fname)
  mageck = read_mageck(fname)
  mageck = inner_join(mageck, oracle_mapping, by = 'gene_id')
  benchmark = compute_benchmark_statistics(mageck)
  dplyr::mutate(benchmark, bin_configuration = info$bin_configuration[1],
    seed = info$seed,
    moi = info$moi, effect_size_type = info$effect_size_type,
    n_cells = info$n_cells)
}

# now process mageck
read_mageck = function(fname) {
  df = read.table(fname, sep = '\t', header = TRUE)
  df = dplyr::mutate(df, value = pmin(pos.fdr, neg.fdr), id = as.character(id))
  df = dplyr::select(df, gene_id = id, value)
  df
}

nimble_fname_to_mageck = function(fname) {
  sub('.rds$', '.gene_summary.txt', basename(fname))
}

wb_recode2 = function (wb_summary, wo, extract_regex = "gene_inclusion")
{
    if (!is.null(wb_summary$all.chains)) {
        s = wb_summary$all.chains
    }
    gi = data.frame(s, mapping = rownames(s))
    gi = dplyr::filter(gi, grepl(extract_regex, mapping))
    gi = dplyr::mutate(gi, mapping = sub(paste0(extract_regex,
        "\\["), "", mapping))
    gi = dplyr::mutate(gi, mapping = sub("\\]", "", mapping))
    gi = dplyr::mutate(gi, mapping = as.integer(mapping))
    gm = dplyr::distinct(dplyr::select(wo$test_guide_names, gene,
        mapping))
    gi = inner_join(gi, gm, by = "mapping")
    gi
}

read_waterbear_results_2_bins = function(metadata_row) {
  info = metadata_row
  # x = readRDS(info$filename)
  x = readRDS(info$summary_filename)

  wo = readRDS(paste0('results/truth/wo_', basename(info$filename)[1]))
  # gene_inclusion = waterbear:::wb_(x, wo)
  gene_inclusion = wb_recode2(x, wo)
  gene_inclusion = dplyr::mutate(gene_inclusion, value = 1 - Mean)
  gene_inclusion = dplyr::rename(gene_inclusion, gene_id = gene)

  gene_inclusion = inner_join(gene_inclusion, oracle_mapping, by = c('gene_id'))
  gene_inclusion = arrange(gene_inclusion, value)

  benchmark = compute_benchmark_statistics(gene_inclusion)
  dplyr::mutate(benchmark, bin_configuration = info$bin_configuration[1],
    seed = info$seed,
    moi = info$moi, effect_size_type = info$effect_size_type,
    n_cells = info$n_cells)
}


########################################################################
# 2 bin mode
########################################################################

fnames = Sys.glob('results/nimble/vary_coverage_2_bins*.rds')


parse_coverage_2_bins_fname = function(x) {
  x = sub('vary_coverage_2_bins_', '', x)
  x = sub('.rds', '', x)
  x = strsplit(x, '_')
  ret = lapply(x,
    function(s) {
      data.frame(
        seed = s[1],
        n_guides_per_target = s[2],
        n_total_genes = s[3],
        n_control_guides = s[4],
        p_gene_effects = s[5],
        n_replicates = s[6],
        n_mcmc_samples = s[7],
        n_mcmc_warmup = s[8],
        bin_configuration = s[9],
        n_cells = s[10],
        effect_size_type = s[11],
        dispersion = s[12],
        n_mcmc_chains = s[13],
        stringsAsFactors = FALSE)
    })
  dplyr::bind_rows(ret)
}

metadata = parse_coverage_2_bins_fname(basename(fnames))
metadata = mutate(metadata, filename = fnames)

########################################################################
# p 0.10
########################################################################

metadata = dplyr::filter(metadata, !grepl('_summary.rds', filename))
metadata = dplyr::mutate(metadata,
  summary_filename = paste0(dirname(filename), '/', sub('.rds', '_summary.rds', basename(filename))))
metadata = dplyr::filter(metadata, file.exists(summary_filename))

metadata1000 = dplyr::filter(metadata, n_total_genes == 1000,
  bin_configuration == 1, p_gene_effects == '0.10')
metadata1000 = filter(metadata1000, effect_size_type == 'equalmixture')

# metadata1000 = dplyr::filter(metadata, n_mcmc_samples == 50000, n_total_genes == 1000,
#   effect_size_type == 'easy', bin_configuration == 1)

p_gene_effects = 0.10
n_total_genes = 1000

# sim_truth = readRDS(paste0('results/truth/sim_', basename(metadata1000$filename)[1]))
sim = readRDS(paste0('results/truth/nimble_', basename(metadata1000$filename)[1]))
# sim = readRDS('results/truth/nimble_vary_coverage_2_bins_2_4_1000_100_0.75_3_50000_25000_3_1025000_moderate_200.00_4.rds')


oracle_mapping = generate_oracle_mapping(sim$const$guide_to_gene,
  p_gene_effects, n_total_genes)

p_gene_effects = '0.10'
n_total_genes = '1000'
bin_configuration = '1'
n_replicates = '3'

tmp_filter = data.frame(
  p_gene_effects = p_gene_effects,
  n_total_genes = n_total_genes,
  bin_configuration = bin_configuration,
  n_replicates = n_replicates
  )

metadata_filter = dplyr::semi_join(metadata1000, tmp_filter)

# metadata_filter = filter(metadata_filter, n_cells == '205000000')


res = lapply(1:nrow(metadata_filter),
  function(i) {
    print(i)
    info = metadata_filter[i, ]
    read_waterbear_results_2_bins(info)
  })
res = bind_rows(res)
res = mutate(res, method = 'waterbear')
res = mutate(res, coverage = as.numeric(n_cells) / (as.integer(n_total_genes) * 4 + 100))

mageck_results = lapply(1:nrow(metadata_filter),
  function(i) {
    print(i)
    info = metadata_filter[i, ]
    process_mageck_results(info)
  })
mageck_results = bind_rows(mageck_results)
mageck_results = mutate(mageck_results, method = 'MAGeCK')
mageck_results = mutate(mageck_results, coverage = as.numeric(n_cells) / (as.integer(n_total_genes) * 4 + 100))

cbb <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

tmp = dplyr::filter(res, coverage == '500')
tmp$true_effect %>% sum

s_res = group_by(res, effect_size_type, n_cells, rank)
s_res = summarize(s_res, mean_fdr = mean(fdr), sd_fdr = sd(fdr),
  mean_sensitivity = mean(sensitivity), sd_sensitivity = sd(sensitivity))
s_res = mutate(s_res, coverage = as.numeric(n_cells) / (as.integer(n_total_genes) * 4 + 100))

# XXX: this is it!
saveRDS(s_res, '2_bins_coverage.rds')


base = paste0('vary_coverage_2_bins_', p_gene_effects, '_', n_total_genes, '_',
  bin_configuration, '_', n_replicates)

# tmp = dplyr::filter(res, bin_configuration == 1, effect_size_type == 'easy')
# p = ggplot(tmp, aes(fdr, tp, color = moi), alpha = 0.2)
p = ggplot(res, aes(fdr, sensitivity, color = factor(coverage)), alpha = 0.2)
p = p + geom_point()
# p = p + geom_line(size = 1.5)
p = p + geom_path(size = 1.5)
p = p + xlim(0, 0.15)
p = p + xlab('true false discovery rate')
p = p + scale_color_manual(values = cbb)
p = p + facet_wrap(~effect_size_type)
p1 = p
# save_plot(paste0('img/', base, '.pdf'), p1, base_height = 7)
save_plot(paste0('img/ajkljskl', base, '.pdf'), p1, base_height = 7)

saveRDS(res, file = paste0('results/', base, '.rds'))

# tmp = dplyr::filter(res, bin_configuration == 1, effect_size_type == 'easy')
# p = ggplot(tmp, aes(fdr, tp, color = moi), alpha = 0.2)
p = ggplot(s_res, aes(mean_fdr, mean_sensitivity, color = factor(coverage)), alpha = 0.2)
p = p + geom_point()
# p = p + geom_line(size = 1.5)
p = p + geom_path(size = 1.5)
p = p + xlim(0, 0.15)
p = p + xlab('true false discovery rate')
p = p + scale_color_manual(values = cbb)
p = p + facet_wrap(~effect_size_type)
p1 = p
save_plot('tmp.pdf', p1, base_height = 7)
# save_plot(paste0('img/summary_', base, '.pdf'), p1, base_height = 7)

########################################################################
# p 0.25
########################################################################

metadata1000 = dplyr::filter(metadata, n_total_genes == 1000,
  bin_configuration == 1, p_gene_effects == '0.25')

p_gene_effects = 0.25
n_total_genes = 1000

# sim_truth = readRDS(paste0('results/truth/sim_', basename(metadata1000$filename)[1]))
sim = readRDS(paste0('results/truth/nimble_', basename(metadata1000$filename)[1]))
# sim = readRDS('results/truth/nimble_vary_coverage_2_bins_2_4_1000_100_0.75_3_50000_25000_3_1025000_moderate_200.00_4.rds')

oracle_mapping = generate_oracle_mapping(sim$const$guide_to_gene,
  p_gene_effects, n_total_genes)

p_gene_effects = as.character(p_gene_effects)
n_total_genes = '1000'
bin_configuration = '1'
n_replicates = '3'

tmp_filter = data.frame(
  p_gene_effects = p_gene_effects,
  n_total_genes = n_total_genes,
  bin_configuration = bin_configuration,
  n_replicates = n_replicates
  )

metadata_filter = dplyr::semi_join(metadata1000, tmp_filter)


res = lapply(1:nrow(metadata_filter),
  function(i) {
    print(i)
    info = metadata_filter[i, ]
    read_waterbear_results(info)
  })
res = bind_rows(res)
res = mutate(res, method = 'waterbear')
res = mutate(res, coverage = as.numeric(n_cells) / (as.integer(n_total_genes) * 4 + 100))

mageck_results = lapply(1:nrow(metadata_filter),
  function(i) {
    print(i)
    info = metadata_filter[i, ]
    process_mageck_results(info)
  })
mageck_results = bind_rows(mageck_results)
mageck_results = mutate(mageck_results, method = 'MAGeCK')
mageck_results = mutate(mageck_results, coverage = as.numeric(n_cells) / (as.integer(n_total_genes) * 4 + 100))

cbb <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

base = paste0('vary_coverage_2_bins_', p_gene_effects, '_', n_total_genes, '_',
  bin_configuration, '_', n_replicates)

# tmp = dplyr::filter(res, bin_configuration == 1, effect_size_type == 'easy')
# p = ggplot(tmp, aes(fdr, tp, color = moi), alpha = 0.2)
p = ggplot(res, aes(fdr, sensitivity, color = factor(coverage)), alpha = 0.2)
p = p + geom_point()
# p = p + geom_line(size = 1.5)
p = p + geom_path(size = 1.5)
p = p + xlim(0, 0.15)
p = p + xlab('true false discovery rate')
p = p + scale_color_manual(values = cbb)
p = p + facet_wrap(~effect_size_type)
p1 = p
save_plot(paste0('img/', base, '.pdf'), p1, base_height = 7)

saveRDS(res, file = paste0('results/', base, '.rds'))





########################################################################

metadata1000 = dplyr::filter(metadata, n_mcmc_samples == 50000, n_total_genes == 1000,
  effect_size_type == 'difficult', bin_configuration == 1)

p_gene_effects = 0.10
n_total_genes = 1000

sim = readRDS(paste0('results/coverage/truth/nimble_vary_coverage_2_4_1000_100_0.50_3_50000_25000_3_410000_easy_200.00_4.rds'))

oracle_mapping = generate_oracle_mapping(sim$const$guide_to_gene,
  p_gene_effects, n_total_genes)

p_gene_effects = '0.10'
n_total_genes = '1000'
n_mcmc_samples = '50000'
n_mcmc_warmup = '25000'
bin_configuration = '1'

tmp_filter = data.frame(
  p_gene_effects = p_gene_effects,
  n_total_genes = n_total_genes,
  n_mcmc_samples = n_mcmc_samples,
  n_mcmc_warmup = n_mcmc_warmup,
  bin_configuration = bin_configuration
  )

metadata_filter = dplyr::semi_join(metadata1000, tmp_filter)

# mageck_results = lapply(1:nrow(metadata_filter),
#   function(i) {
#     print(i)
#     info = metadata_filter[i, ]
#     process_mageck_results(info)
#   })
# mageck_results = bind_rows(mageck_results)
# mageck_results = mutate(mageck_results, method = 'MAGeCK')

res = lapply(1:nrow(metadata_filter),
  function(i) {
    print(i)
    info = metadata_filter[i, ]
    read_nimble_results(info)
  })
res = bind_rows(res)
res = mutate(res, method = 'waterbear')

res = mutate(res, coverage = as.numeric(n_cells) / (as.integer(n_total_genes) * 4 + 100))

cbb <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")


# tmp = dplyr::filter(res, bin_configuration == 1, effect_size_type == 'easy')
# p = ggplot(tmp, aes(fdr, tp, color = moi), alpha = 0.2)
p = ggplot(res, aes(fdr, sensitivity, color = factor(coverage)), alpha = 0.2)
p = p + geom_point()
# p = p + geom_line(size = 1.5)
p = p + geom_path(size = 1.5)
p = p + xlim(0, 0.15)
p = p + xlab('true false discovery rate')
p = p + scale_color_manual(values = cbb)
p1 = p
save_plot('img/coverage_difficult_2_bins.pdf', p1, base_height = 7)

# tmp = dplyr::filter(res, bin_configuration == 1, effect_size_type == 'easy')
# p = ggplot(tmp, aes(fdr, tp, color = moi), alpha = 0.2)
p = ggplot(res, aes(rank, sensitivity, color = factor(coverage)), alpha = 0.2)
p = p + geom_point()
# p = p + geom_line(size = 1.5)
p = p + geom_path(size = 1.5)
p = p + xlim(0, 0.15)
p = p + xlab('rank')
p = p + scale_color_manual(values = cbb)
p1 = p
save_plot('img/coverage_difficult_2_bins_rank.pdf', p1, base_height = 7)


# update: the new 2 bin mode with configuration 4

metadata = dplyr::filter(metadata, !grepl('_summary.rds', filename))
metadata = dplyr::mutate(metadata,
  summary_filename = paste0(dirname(filename), '/', sub('.rds', '_summary.rds', basename(filename))))
metadata = dplyr::filter(metadata, file.exists(summary_filename))

metadata1000 = dplyr::filter(metadata, n_total_genes == 1000,
  bin_configuration == 4, p_gene_effects == '0.10')
metadata1000 = filter(metadata1000, effect_size_type == 'equalmixture')

# metadata1000 = dplyr::filter(metadata, n_mcmc_samples == 50000, n_total_genes == 1000,
#   effect_size_type == 'easy', bin_configuration == 1)

p_gene_effects = 0.10
n_total_genes = 1000

# sim_truth = readRDS(paste0('results/truth/sim_', basename(metadata1000$filename)[1]))
sim = readRDS(paste0('results/truth/nimble_', basename(metadata1000$filename)[1]))
# sim = readRDS('results/truth/nimble_vary_coverage_2_bins_2_4_1000_100_0.75_3_50000_25000_3_1025000_moderate_200.00_4.rds')


oracle_mapping = generate_oracle_mapping(sim$const$guide_to_gene,
  p_gene_effects, n_total_genes)

p_gene_effects = '0.10'
n_total_genes = '1000'
bin_configuration = '4'
n_replicates = '3'

tmp_filter = data.frame(
  p_gene_effects = p_gene_effects,
  n_total_genes = n_total_genes,
  bin_configuration = bin_configuration,
  n_replicates = n_replicates
  )

metadata_filter = dplyr::semi_join(metadata1000, tmp_filter)

# metadata_filter = filter(metadata_filter, n_cells == '205000000')


res = lapply(1:nrow(metadata_filter),
  function(i) {
    print(i)
    info = metadata_filter[i, ]
    read_waterbear_results_2_bins(info)
  })
res = bind_rows(res)
res = mutate(res, method = 'waterbear')
res = mutate(res, coverage = as.numeric(n_cells) / (as.integer(n_total_genes) * 4 + 100))

mageck_results = lapply(1:nrow(metadata_filter),
  function(i) {
    print(i)
    info = metadata_filter[i, ]
    process_mageck_results(info)
  })
mageck_results = bind_rows(mageck_results)
mageck_results = mutate(mageck_results, method = 'MAGeCK')
mageck_results = mutate(mageck_results, coverage = as.numeric(n_cells) / (as.integer(n_total_genes) * 4 + 100))

cbb <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

tmp = dplyr::filter(res, coverage == '500')
tmp$true_effect %>% sum

s_res = group_by(res, effect_size_type, n_cells, rank)
s_res = summarize(s_res, mean_fdr = mean(fdr), sd_fdr = sd(fdr),
  mean_sensitivity = mean(sensitivity), sd_sensitivity = sd(sensitivity))
s_res = mutate(s_res, coverage = as.numeric(n_cells) / (as.integer(n_total_genes) * 4 + 100))

# XXX: this is it!
saveRDS(s_res, '2_bins_coverage_bin4.rds')

# update: the new 2 bin mode with configuration 5

metadata = dplyr::filter(metadata, !grepl('_summary.rds', filename))
metadata = dplyr::mutate(metadata,
  summary_filename = paste0(dirname(filename), '/', sub('.rds', '_summary.rds', basename(filename))))
metadata = dplyr::filter(metadata, file.exists(summary_filename))

metadata1000 = dplyr::filter(metadata, n_total_genes == 1000,
  bin_configuration == 5, p_gene_effects == '0.10')
metadata1000 = filter(metadata1000, effect_size_type == 'equalmixture')

# metadata1000 = dplyr::filter(metadata, n_mcmc_samples == 50000, n_total_genes == 1000,
#   effect_size_type == 'easy', bin_configuration == 1)

p_gene_effects = 0.10
n_total_genes = 1000

# sim_truth = readRDS(paste0('results/truth/sim_', basename(metadata1000$filename)[1]))
sim = readRDS(paste0('results/truth/nimble_', basename(metadata1000$filename)[1]))
# sim = readRDS('results/truth/nimble_vary_coverage_2_bins_2_4_1000_100_0.75_3_50000_25000_3_1025000_moderate_200.00_4.rds')


oracle_mapping = generate_oracle_mapping(sim$const$guide_to_gene,
  p_gene_effects, n_total_genes)

p_gene_effects = '0.10'
n_total_genes = '1000'
bin_configuration = '5'
n_replicates = '2'

tmp_filter = data.frame(
  p_gene_effects = p_gene_effects,
  n_total_genes = n_total_genes,
  bin_configuration = bin_configuration,
  n_replicates = n_replicates
  )

metadata_filter = dplyr::semi_join(metadata1000, tmp_filter)
# metadata_filter = filter(metadata_filter, n_cells == '205000000')


res = lapply(1:nrow(metadata_filter),
  function(i) {
    print(i)
    info = metadata_filter[i, ]
    read_waterbear_results_2_bins(info)
  })
res = bind_rows(res)
res = mutate(res, method = 'waterbear')
res = mutate(res, coverage = as.numeric(n_cells) / (as.integer(n_total_genes) * 4 + 100))

mageck_results = lapply(1:nrow(metadata_filter),
  function(i) {
    print(i)
    info = metadata_filter[i, ]
    process_mageck_results(info)
  })
mageck_results = bind_rows(mageck_results)
mageck_results = mutate(mageck_results, method = 'MAGeCK')
mageck_results = mutate(mageck_results, coverage = as.numeric(n_cells) / (as.integer(n_total_genes) * 4 + 100))

cbb <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

tmp = dplyr::filter(res, coverage == '500')
tmp$true_effect %>% sum

s_res = group_by(res, effect_size_type, n_cells, rank)
s_res = summarize(s_res, mean_fdr = mean(fdr), sd_fdr = sd(fdr),
  mean_sensitivity = mean(sensitivity), sd_sensitivity = sd(sensitivity))
s_res = mutate(s_res, coverage = as.numeric(n_cells) / (as.integer(n_total_genes) * 4 + 100))

# XXX: this is it!
saveRDS(s_res, '2_bins_coverage_bin5.rds')
