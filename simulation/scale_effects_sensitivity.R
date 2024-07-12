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

get_fdr_table_nimble = function(metadata_row) {
  info = metadata_row
  x = readRDS(info$filename)
  all_chains = x$summary[['all.chains']]
  var_names = rownames(all_chains)
  gene_inclusion = all_chains[grepl('gene_inclusion', var_names), ]
  gene_inclusion = data.frame(variable = rownames(gene_inclusion), gene_inclusion)
  gene_inclusion = dplyr::mutate(gene_inclusion, gene_id = strip_text(variable),
    value = 1 - Mean)
  gene_inclusion = inner_join(gene_inclusion, oracle_mapping, by = c('gene_id'))
  gene_inclusion = arrange(gene_inclusion, value)
  gene_inclusion
}

read_nimble_results = function(metadata_row) {
  info = metadata_row
  x = readRDS(info$filename)
  all_chains = x$summary[['all.chains']]
  var_names = rownames(all_chains)
  gene_inclusion = all_chains[grepl('gene_inclusion', var_names), ]
  gene_inclusion = data.frame(variable = rownames(gene_inclusion), gene_inclusion)
  gene_inclusion = dplyr::mutate(gene_inclusion, gene_id = strip_text(variable),
    value = 1 - Mean)
  gene_inclusion = inner_join(gene_inclusion, oracle_mapping, by = c('gene_id'))
  gene_inclusion = arrange(gene_inclusion, value)
  benchmark = compute_benchmark_statistics(gene_inclusion)
  dplyr::mutate(benchmark, bin_configuration = info$bin_configuration[1],
    seed = info$seed,
    moi = info$moi, effect_size_type = info$effect_size_type,
    n_cells = info$n_cells)
}

########################################################################
# looking at coverage
fnames = Sys.glob('results/nimble/scale_effects_*.rds')

parse_scale_effects_fname = function(x) {
  x = sub('scale_effects_', '', x)
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
        # moi = s[10],
        n_cells = s[10],
        effect_size_type = s[11],
        dispersion = s[12],
        n_mcmc_chains = s[13],
        stringsAsFactors = FALSE)
    })
  dplyr::bind_rows(ret)
}

metadata = parse_scale_effects_fname(basename(fnames))
metadata = mutate(metadata, filename = fnames)

metadata1000 = dplyr::filter(metadata, n_mcmc_samples == 50000, n_total_genes == 1000,
  seed == 2)


p_gene_effects = 0.10
n_total_genes = 1000
sim = readRDS(paste0('results/truth/nimble_', basename(dplyr::filter(metadata1000, p_gene_effects == p_gene_effects)$filename[1])))
oracle_mapping = generate_oracle_mapping(sim$const$guide_to_gene,
  p_gene_effects, n_total_genes)

p_gene_effects = '0.10'
n_total_genes = '1000'
n_mcmc_samples = '50000'
n_mcmc_warmup = '25000'
# bin_configuration = '1'

tmp_filter = data.frame(
  p_gene_effects = p_gene_effects,
  n_total_genes = n_total_genes,
  n_mcmc_samples = n_mcmc_samples,
  n_mcmc_warmup = n_mcmc_warmup
  # bin_configuration = bin_configuration
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

all_results = res
# all_results = bind_rows(res, mageck_results)

p = ggplot(all_results, aes(fdr, sensitivity, color = seed, linetype = method), alpha = 0.2)
p = p + geom_point()
p = p + geom_path(size = 1.5)
p = p + xlim(0, 0.15)
p = p + xlab('true false discovery rate')
p = p + scale_color_manual(values = cbb)
p = p + facet_wrap(~n_cells, labeller = label_both)
p1 = p

save_plot('img/tmp.pdf', p1, base_height = 7)

all_results = dplyr::group_by(all_results, method, seed, effect_size_type, n_cells,
  bin_configuration)
alpha = 0.10

summary_statistics = summarize(all_results,
  sensitivity = sum((value < alpha) & true_effect) / sum(true_effect),
  fpr = sum((value < alpha) & !true_effect) / sum(!true_effect),
  fdr = sum((value < alpha) & !true_effect) / sum(value < alpha))

p = ggplot(summary_statistics, aes(effect_size_type, sensitivity, group = bin_configuration, color = bin_configuration))
p = p + geom_path(size = 1.25)
p = p + scale_color_manual(values = cbb)
save_plot('img/scale_effects_sensitivity.pdf', p, base_height = 7)

p = ggplot(summary_statistics, aes(effect_size_type, fdr, group = bin_configuration, color = bin_configuration))
p = p + geom_path(size = 1.25)
p = p + scale_color_manual(values = cbb)
save_plot('img/scale_effects_fdr.pdf', p, base_height = 7)

save_plot('img/tmp.pdf', p1, base_height = 7)

p_gene_effects = 0.50
n_total_genes = 1000
sim = readRDS(paste0('results/truth/nimble_', basename(dplyr::filter(metadata1000, p_gene_effects == p_gene_effects)$filename[1])))
oracle_mapping = generate_oracle_mapping(sim$const$guide_to_gene,
  p_gene_effects, n_total_genes)

p_gene_effects = '0.50'
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

mageck_results = lapply(1:nrow(metadata_filter),
  function(i) {
    print(i)
    info = metadata_filter[i, ]
    process_mageck_results(info)
  })
mageck_results = bind_rows(mageck_results)
mageck_results = mutate(mageck_results, method = 'MAGeCK')

res = lapply(1:nrow(metadata_filter),
  function(i) {
    print(i)
    info = metadata_filter[i, ]
    read_nimble_results(info)
  })
res = bind_rows(res)
res = mutate(res, method = 'waterbear')

all_results = bind_rows(res, mageck_results)

p = ggplot(all_results, aes(fdr, sensitivity, color = seed, linetype = method), alpha = 0.2)
p = p + geom_point()
p = p + geom_path(size = 1.5)
p = p + xlim(0, 0.15)
p = p + xlab('true false discovery rate')
p = p + scale_color_manual(values = cbb)
p = p + facet_wrap(~n_cells, labeller = label_both)
p1 = p

all_results = dplyr::group_by(all_results, method, seed, effect_size_type, n_cells)
alpha = 0.10

summary_statistics = summarize(all_results,
  sensitivity = sum((value < alpha) & true_effect) / sum(true_effect),
  fpr = sum((value < alpha) & !true_effect) / sum(!true_effect),
  fdr = sum((value < alpha) & !true_effect) / sum(value < alpha))


# effect_size_type = 'easy'
res = lapply(1:nrow(metadata500),
# res = lapply(1,
  function(i) {
    print(i)
    info = metadata500[i, ]
    x = readRDS(info$filename)
    all_chains = x$summary[['all.chains']]
    var_names = rownames(all_chains)
    gene_inclusion = all_chains[grepl('gene_inclusion', var_names), ]
    gene_inclusion = data.frame(variable = rownames(gene_inclusion), gene_inclusion)
    gene_inclusion = dplyr::mutate(gene_inclusion, gene_id = strip_text(variable),
      value = 1 - Mean)
    gene_inclusion = inner_join(gene_inclusion, oracle_mapping, by = c('gene_id'))
    gene_inclusion = arrange(gene_inclusion, value)
    benchmark = compute_benchmark_statistics(gene_inclusion)
    dplyr::mutate(benchmark, bin_configuration = info$bin_configuration[1],
      seed = info$seed,
      moi = info$moi, effect_size_type = info$effect_size_type,
      n_cells = info$n_cells)
  })
res = bind_rows(res)

tmp = dplyr::filter(res, n_cells == '500000', bin_configuration == 1)
p = ggplot(tmp, aes(fdr, sensitivity, color = moi), alpha = 0.2)
p = p + geom_point()
p = p + geom_line()
p = p + facet_wrap(~effect_size_type)
p
dev.off()

p = ggplot(tmp, aes(rank, sensitivity, color = moi), alpha = 0.2)
p = p + geom_point()
p = p + geom_line()
p = p + facet_wrap(~effect_size_type)
p = p + xlim(0, 100)
p
dev.off()

########################################################################

metadata1000 = dplyr::filter(metadata, n_mcmc_samples == 50000, n_total_genes == 1000,
  n_cells == '50000', effect_size_type == 'easy', bin_configuration == 1)


sim = readRDS('results/truth/nimble_unimodal_moi_2_4_1000_100_0.10_3_50000_25000_3_2_1000000_moderate_200.00_4.rds')
oracle_mapping = data.frame(sim_id = unique(sim$const$guide_to_gene))
oracle_mapping = mutate(oracle_mapping,
  gene_id = as.character(sim_id))
oracle_mapping = mutate(oracle_mapping,
  true_effect = FALSE)
oracle_mapping$true_effect[1:(0.10 * 1000)] = TRUE

res = lapply(1:nrow(metadata1000),
# res = lapply(1,
  function(i) {
    print(i)
    info = metadata1000[i, ]
    x = readRDS(info$filename)
    all_chains = x$summary[['all.chains']]
    var_names = rownames(all_chains)
    gene_inclusion = all_chains[grepl('gene_inclusion', var_names), ]
    gene_inclusion = data.frame(variable = rownames(gene_inclusion), gene_inclusion)
    gene_inclusion = dplyr::mutate(gene_inclusion, gene_id = strip_text(variable),
      value = 1 - Mean)
    gene_inclusion = inner_join(gene_inclusion, oracle_mapping, by = c('gene_id'))
    gene_inclusion = arrange(gene_inclusion, value)
    benchmark = compute_benchmark_statistics(gene_inclusion)
    dplyr::mutate(benchmark, bin_configuration = info$bin_configuration[1],
      seed = info$seed,
      moi = info$moi, effect_size_type = info$effect_size_type,
      n_cells = info$n_cells)
  })
res = bind_rows(res)

cbb <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")


tmp = dplyr::filter(res, bin_configuration == 1, effect_size_type == 'easy')
# p = ggplot(tmp, aes(fdr, tp, color = moi), alpha = 0.2)
p = ggplot(tmp, aes(fdr, sensitivity, color = moi), alpha = 0.2)
p = p + geom_point()
# p = p + geom_line(size = 1.5)
p = p + geom_path(size = 1.5)
p = p + xlim(0, 0.15)
p = p + xlab('true false discovery rate')
p = p + theme(legend.position = c(0.2, 0.8))
p = p + scale_color_manual(values = cbb)
p1 = p

dev.off()

p = ggplot(tmp, aes(rank, sensitivity, color = moi), alpha = 0.2)
p = p + geom_point()
p = p + geom_line(size = 1.5)
p = p + xlim(0, 200)
p = p + scale_color_manual(values = cbb)
p = p + theme(legend.position = c(0.2, 0.8))
p2 = p

p = plot_grid(p1, p2, nrow = 1)
save_plot('img/easy_effect_moi_low_cells.pdf', p, base_height = 7)

dev.off()


#####
### look at many moi

metadata1000 = dplyr::filter(metadata, n_mcmc_samples == 50000, n_total_genes == 1000,
  effect_size_type == 'easy', bin_configuration == 1)

sim = readRDS('results/truth/nimble_unimodal_moi_2_4_1000_100_0.10_3_50000_25000_3_2_1000000_moderate_200.00_4.rds')
oracle_mapping = data.frame(sim_id = unique(sim$const$guide_to_gene))
oracle_mapping = mutate(oracle_mapping,
  gene_id = as.character(sim_id))
oracle_mapping = mutate(oracle_mapping,
  true_effect = FALSE)
oracle_mapping$true_effect[1:(0.10 * 1000)] = TRUE

res = lapply(1:nrow(metadata1000),
# res = lapply(1,
  function(i) {
    print(i)
    info = metadata1000[i, ]
    x = readRDS(info$filename)
    all_chains = x$summary[['all.chains']]
    var_names = rownames(all_chains)
    gene_inclusion = all_chains[grepl('gene_inclusion', var_names), ]
    gene_inclusion = data.frame(variable = rownames(gene_inclusion), gene_inclusion)
    gene_inclusion = dplyr::mutate(gene_inclusion, gene_id = strip_text(variable),
      value = 1 - Mean)
    gene_inclusion = inner_join(gene_inclusion, oracle_mapping, by = c('gene_id'))
    gene_inclusion = arrange(gene_inclusion, value)
    benchmark = compute_benchmark_statistics(gene_inclusion)
    dplyr::mutate(benchmark, bin_configuration = info$bin_configuration[1],
      seed = info$seed,
      moi = info$moi, effect_size_type = info$effect_size_type,
      n_cells = info$n_cells)
  })
res = bind_rows(res)

res = group_by(res, rank, bin_configuration, n_cells, effect_size_type, moi)

s_res = summarize(res, fdr = mean(fdr), fdr_sd = sd(fdr),
  sensitivity = mean(sensitivity), sensitivity_sd = sd(sensitivity))

i = 1
tmp = dplyr::filter(s_res, bin_configuration == i)
p = ggplot(tmp, aes(fdr, sensitivity, color = factor(moi)))
p = p + geom_point()
p = p + geom_line()
p = p + facet_wrap(~n_cells, 'label_both')
p = p + scale_color_brewer(type = 'div', palette = 'Spectral')
p = p + theme(legend.position = c(0.8, 0.2))
p = p + xlim(0, 0.2)
p = p + ylim(0.2, 1)
p
dev.off()



########################################################################

debugonce(parse_coverage)
coverage_metadata = parse_coverage(basename(fnames))

coverage500 = dplyr::filter(coverage_metadata, n_iter == 50000, n_genes == 500)

res = lapply(1:nrow(coverage500),
# res = lapply(1,
  function(i) {
    print(i)
    info = coverage500[i, ]
    x = readRDS(paste0('results/nimble/', info$filename))
    all_chains = x$summary[['all.chains']]
    var_names = rownames(all_chains)
    gene_inclusion = all_chains[grepl('gene_inclusion', var_names), ]
    gene_inclusion = data.frame(variable = rownames(gene_inclusion), gene_inclusion)
    gene_inclusion = dplyr::mutate(gene_inclusion, gene_id = strip_text(variable),
      value = 1 - Mean)
    gene_inclusion = inner_join(gene_inclusion, oracle_mapping, by = c('gene_id'))
    gene_inclusion = arrange(gene_inclusion, value)
    benchmark = compute_benchmark_statistics(gene_inclusion)
    dplyr::mutate(benchmark, bin_configuration = info$bin_configuration[1],
      n_cells_per_guide = info$n_cells_per_guide[1], n_samples = info$n_samples[1])
  })
res = bind_rows(res)

p = ggplot(res, aes(fdr, sensitivity, color = factor(n_cells_per_guide)))
p = p + geom_point()
p = p + geom_line()
p = p + facet_wrap(~bin_configuration + n_samples)
p = p + theme(legend.position = c(0.8, 0.2))
p = p + xlim(0, 0.2)
p
dev.off()


p = ggplot(res, aes(rank, sensitivity, color = factor(n_cells_per_guide)))
p = p + geom_point()
p = p + geom_line()
# p = p + facet_wrap(~bin_configuration)
p = p + facet_wrap(~bin_configuration + n_samples)
p = p + theme(legend.position = c(0.8, 0.2))
p = p + xlim(0, 100)
p
dev.off()

# TODO: given a filename, compute sensitivity
tmp = readRDS(fnames[14])

all_chains = tmp$summary[['all.chains']]
var_names = rownames(all_chains)
gene_inclusion = all_chains[grepl('gene_inclusion', var_names), ]
gene_inclusion = data.frame(variable = rownames(gene_inclusion), gene_inclusion)
gene_inclusion = dplyr::mutate(gene_inclusion, gene_id = strip_text(variable),
  value = 1 - Mean)
gene_inclusion = dplyr::arrange(gene_inclusion, value)

dim(gene_effect_sizes)

source('dm_simple.R')

simple_parameters = list(
  sigma_guide = 0.08732188,
  sigma_gene = 0.2381384,
  # sigma_gene = 0.96647413,
  dispersion = 200,
  n_samples = 3,
  p_effects = 0.05,
  n_genes = 500,
  n_guides_per_gene = 4,
  # bin_sizes = c(0.15, 0.7, 0.15),
  bin_sizes = c(0.15, 0.35, 0.35, 0.15),
  n_cells_per_guide = 2000,
  n_reporter_cells = 1e6,
  # guide_rate_dispersion = 6.3
  # guide_rate_dispersion = 200,
  guide_rate_dispersion = guide_rate_total_mass,
  n_control_guides = 100
  )

sim = simulate_equal_bin_sizes(simple_parameters)



oracle_mapping = data.frame(sim_id = unique(sim$guide_to_gene))
oracle_mapping = mutate(oracle_mapping,
  gene_id = as.character(as.integer(factor(as.integer(sim_id)))))
oracle_mapping = mutate(oracle_mapping,
  true_effect = FALSE)
oracle_mapping$true_effect[1:length(sim$gene_effect_sizes)] = TRUE



gene_inclusion = inner_join(gene_inclusion, oracle_mapping, by = c('gene_id'))
gene_inclusion = arrange(gene_inclusion, value)

tmp = compute_benchmark_statistics(gene_inclusion)




# TODO: look at convergence of each sample
# TODO: look at bias in parameters

# TODO: look at the different effect sizes

# looking at coverage from the observed effects
fnames = Sys.glob('results/nimble/coverage_observed_effect*.rds')
coverage_metadata = parse_coverage_observed_effect(basename(fnames))
coverage_metadata = mutate(coverage_metadata, path = fnames)

coverage500 =  dplyr::filter(coverage_metadata, n_iter == 50000, n_genes == 500)

res = lapply(1:nrow(coverage500),
# res = lapply(1,
  function(i) {
    print(i)
    info = coverage500[i, ]
    x = readRDS(paste0('results/nimble/', info$filename))
    all_chains = x$summary[['all.chains']]
    var_names = rownames(all_chains)
    gene_inclusion = all_chains[grepl('gene_inclusion', var_names), ]
    gene_inclusion = data.frame(variable = rownames(gene_inclusion), gene_inclusion)
    gene_inclusion = dplyr::mutate(gene_inclusion, gene_id = strip_text(variable),
      value = 1 - Mean)
    gene_inclusion = inner_join(gene_inclusion, oracle_mapping, by = c('gene_id'))
    gene_inclusion = arrange(gene_inclusion, value)
    benchmark = compute_benchmark_statistics(gene_inclusion)
    dplyr::mutate(benchmark, bin_configuration = info$bin_configuration[1],
      n_cells_per_guide = info$n_cells_per_guide[1], n_samples = info$n_samples[1],
      which_effects = info$which_effects[1])
  })
res = bind_rows(res)

# get oracle
sim = readRDS(paste0('results/truth/sim_', coverage500$filename[1])

oracle_mapping = data.frame(sim_id = unique(sim$guide_to_gene))
oracle_mapping = mutate(oracle_mapping,
  gene_id = as.character(as.integer(factor(as.integer(sim_id)))))
oracle_mapping = mutate(oracle_mapping,
  true_effect = FALSE)
oracle_mapping$true_effect[1:length(sim$gene_effect_sizes)] = TRUE


for (i in 1:3) {
  tmp = dplyr::filter(res, bin_configuration == i)
  p = ggplot(tmp, aes(fdr, sensitivity, color = factor(n_cells_per_guide)))
  p = p + geom_point()
  p = p + geom_line()
  p = p + facet_wrap(~n_samples + which_effects, labeller = 'label_both')
  p = p + theme(legend.position = c(0.8, 0.2))
  p = p + xlim(0, 0.2)
  save_plot(p, file = paste0('img/fdr_sensitivity_observed_effect_bin_', i, '.pdf'),
    base_height = 10)
}

p = ggplot(res, aes(rank, sensitivity, color = factor(n_cells_per_guide)))
p = p + geom_point()
p = p + geom_line()
p = p + facet_wrap(~bin_configuration + n_samples + which_effects, labeller = 'label_both')
p = p + theme(legend.position = c(0.8, 0.2))
p = p + xlim(0, 100)

save_plot(p, file = 'img/rank_sensitivity_observed_effect.pdf',
  base_height = 10)

# coverage at 1000
coverage1000 = dplyr::filter(coverage_metadata, n_iter == 50000, n_genes == 1000,
  which_effects == 'difficult', n_samples == 3, bin_configuration == 1)

sim = readRDS('results/truth/nimble_coverage_observed_effect_easy_2_4_50000_25000_2_0.05_1000_4_1_1000.rds')

oracle_mapping = data.frame(sim_id = unique(sim$const$guide_to_gene))
oracle_mapping = mutate(oracle_mapping,
  gene_id = as.character(sim_id))
oracle_mapping = mutate(oracle_mapping,
  true_effect = FALSE)
oracle_mapping$true_effect[1:(0.05 * 1000)] = TRUE

res = lapply(1:nrow(coverage1000),
# res = lapply(1,
  function(i) {
    print(i)
    info = coverage1000[i, ]
    x = readRDS(paste0('results/nimble/', info$filename))
    all_chains = x$summary[['all.chains']]
    var_names = rownames(all_chains)
    gene_inclusion = all_chains[grepl('gene_inclusion', var_names), ]
    gene_inclusion = data.frame(variable = rownames(gene_inclusion), gene_inclusion)
    gene_inclusion = dplyr::mutate(gene_inclusion, gene_id = strip_text(variable),
      value = 1 - Mean)
    gene_inclusion = inner_join(gene_inclusion, oracle_mapping, by = c('gene_id'))
    gene_inclusion = arrange(gene_inclusion, value)
    benchmark = compute_benchmark_statistics(gene_inclusion)
    dplyr::mutate(benchmark, bin_configuration = info$bin_configuration[1],
      n_cells_per_guide = info$n_cells_per_guide[1], n_samples = info$n_samples[1],
      which_effects = info$which_effects[1])
  })
res = bind_rows(res)

tmp = dplyr::filter(res, bin_configuration == 1, which_effects == 'difficult', n_samples == 3)
p = ggplot(tmp, aes(rank, sensitivity, color = factor(n_cells_per_guide)), alpha = 0.2)
p = p + geom_point()
p = p + geom_path(size = 1.5)
p = p + xlim(0, 200)
# p = p + xlab('true false discovery rate')
p = p + theme(legend.position = c(0.8, 0.2))
p = p + scale_color_manual(values = cbb)
p1 = p
p1
# dev.off()

tmp = dplyr::filter(res, bin_configuration == 1, which_effects == 'difficult', n_samples == 3)
p = ggplot(tmp, aes(fdr, sensitivity, color = factor(n_cells_per_guide)), alpha = 0.2)
# p = ggplot(tmp, aes(rank, sensitivity, color = factor(n_cells_per_guide)), alpha = 0.2)
p = p + geom_point()
p = p + geom_path(size = 1.5)
p = p + theme(legend.position = c(0.2, 0.8))
p = p + scale_color_manual(values = cbb)
p = p + theme(legend.position = c(0.7, 0.2))
p = p + xlab('true false discovery rate')
p = p + xlim(0, 0.15)
p2 = p
# dev.off()

p = plot_grid(p2, p1, nrow = 1)
save_plot('img/difficult_effect_coverage.pdf', p, base_height = 7)


########################################################################
# starting with simple proportion of 10
########################################################################


p_gene_effects = '0.10'
n_total_genes = '1000'
n_mcmc_samples = '50000'
n_mcmc_warmup = '25000'
bin_configuration = '1'
moi = as.character(c(1:5, 7, 10))
# effect_size_type = 'easy'

tmp_filter = data.frame(
  p_gene_effects = p_gene_effects,
  n_total_genes = n_total_genes,
  n_mcmc_samples = n_mcmc_samples,
  n_mcmc_warmup = n_mcmc_warmup,
  bin_configuration = bin_configuration,
  # effect_size_type = effect_size_type,
  moi = moi)

metadata_filter = dplyr::semi_join(metadata, tmp_filter)

tmp_sim = readRDS(paste0('results/truth/nimble_', basename(metadata_filter$filename[1])))

oracle_mapping = generate_oracle_mapping(
  tmp_sim$const$guide_to_gene,
  p_gene_effects, n_total_genes)


res = lapply(1:nrow(metadata_filter),
  function(i) {
    print(i)
    info = metadata_filter[i, ]
    read_nimble_results(info)
  })
res = bind_rows(res)

p = ggplot(res, aes(fdr, sensitivity, color = moi), alpha = 0.2)
p = p + geom_point()
p = p + geom_path(size = 1.5)
p = p + xlim(0, 0.15)
p = p + xlab('true false discovery rate')
p = p + scale_color_manual(values = cbb)
p = p + facet_wrap(~n_cells, labeller = label_both)
p1 = p

img_name = paste0('p_gene_', p_gene_effects, '_n_genes_', n_total_genes,
  '_mcmc_', n_mcmc_samples, '_', n_mcmc_warmup, '_bin_', bin_configuration,
  '_estype_', effect_size_type, '.pdf')

save_plot(paste0('img/sensitivity_moi_', img_name),
  p, base_height = 7)

########################################################################
# 75% of gene effects
########################################################################

p_gene_effects = '0.75'
n_total_genes = '1000'
n_mcmc_samples = '50000'
n_mcmc_warmup = '25000'
bin_configuration = '1'
moi = as.character(c(1:5, 7, 10))
# effect_size_type = 'easy'

tmp_filter = data.frame(
  p_gene_effects = p_gene_effects,
  n_total_genes = n_total_genes,
  n_mcmc_samples = n_mcmc_samples,
  n_mcmc_warmup = n_mcmc_warmup,
  bin_configuration = bin_configuration,
  # effect_size_type = effect_size_type,
  moi = moi)

metadata_filter = dplyr::semi_join(metadata, tmp_filter)

tmp_sim = readRDS(paste0('results/truth/nimble_', basename(metadata_filter$filename[1])))

oracle_mapping = generate_oracle_mapping(
  tmp_sim$const$guide_to_gene,
  p_gene_effects, n_total_genes)

res = lapply(1:nrow(metadata_filter),
  function(i) {
    print(i)
    info = metadata_filter[i, ]
    read_nimble_results(info)
  })
res = bind_rows(res)
res = dplyr::mutate(res, method = 'waterbear')

p = ggplot(res, aes(fdr, sensitivity, color = moi), alpha = 0.2)
p = p + geom_point()
p = p + geom_path(size = 0.75)
p = p + xlim(0, 0.15)
p = p + xlab('true false discovery rate')
p = p + scale_color_manual(values = cbb)
p = p + facet_wrap(~n_cells + effect_size_type, labeller = label_both)
p1 = p

img_name = paste0('p_gene_', p_gene_effects, '_n_genes_', n_total_genes,
  '_mcmc_', n_mcmc_samples, '_', n_mcmc_warmup, '_bin_', bin_configuration,
  '_estype_', 'all', '.pdf')

save_plot(paste0('img/sensitivity_moi_', img_name),
  p, base_height = 7)



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


mageck_results = lapply(1:nrow(metadata_filter),
  function(i) {
    print(i)
    info = metadata_filter[i, ]
    # fname = nimble_fname_to_mageck(info$filename)
    # fname = paste0('results/mageck/', fname)
    # read_mageck(fname)
    process_mageck_results(info)
  })
mageck_results = bind_rows(mageck_results)
mageck_results = mutate(mageck_results, method = 'MAGeCK')

all_results = bind_rows(res, mageck_results)

all_results = dplyr::group_by(all_results, method, seed, moi, effect_size_type, n_cells)

alpha = 0.10

summary_statistics = summarize(all_results,
  sensitivity = sum((value < alpha) & true_effect) / sum(true_effect),
  fpr = sum((value < alpha) & !true_effect) / sum(!true_effect),
  fdr = sum((value < alpha) & !true_effect) / sum(value < alpha))

dplyr::filter(summary_statistics, effect_size_type == 'easy')

p = ggplot(summary_statistics, aes(fpr, sensitivity, color = method, shape = moi))
p = p + geom_point(alpha = 0.5, size = 2)
p = p + facet_wrap(~effect_size_type + n_cells)

save_plot(paste0('img/magic_comparison', img_name),
  p, base_height = 7)
