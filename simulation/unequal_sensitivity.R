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
    sensitivity = tp / sum(true_effect), fdr = fp / rank)

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
fnames = Sys.glob('results/nimble/unequal_*equalmixture*.rds')

parse_unequal_bins_fname = function(x) {
  x = sub('unequal_bins_', '', x)
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

metadata = parse_unequal_bins_fname(basename(fnames))
metadata = mutate(metadata, filename = fnames)

# metadata1000 = dplyr::filter(metadata, n_mcmc_samples == 50000, n_total_genes == 1000,
metadata1000 = dplyr::filter(metadata, n_mcmc_samples == 20000, n_total_genes == 1000,
  seed == 2, p_gene_effects == '0.10', n_cells == '2000000')


p_gene_effects = 0.10
n_total_genes = 1000
sim = readRDS(paste0('results/truth/nimble_', basename(dplyr::filter(metadata1000, p_gene_effects == p_gene_effects)$filename[1])))
oracle_mapping = generate_oracle_mapping(sim$const$guide_to_gene,
  p_gene_effects, n_total_genes)

p_gene_effects = '0.10'
n_total_genes = '1000'
# n_mcmc_samples = '50000'
# n_mcmc_warmup = '25000'
bin_configuration = '1'

tmp_filter = data.frame(
  p_gene_effects = p_gene_effects,
  n_total_genes = n_total_genes,
  # n_mcmc_samples = n_mcmc_samples,
  # n_mcmc_warmup = n_mcmc_warmup,
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

p = ggplot(all_results, aes(fdr, sensitivity, color = method), alpha = 0.2)
p = p + geom_point()
p = p + geom_path(size = 1.5)
p = p + xlim(0, 0.15)
p = p + xlab('true false discovery rate')
p = p + scale_color_manual(values = cbb)
p = p + facet_wrap(~n_cells+effect_size_type, labeller = label_both)
p1 = p
save_plot('img/unequal_sensitivity_p0.10.pdf', p1, base_height = 7)

p = ggplot(res, aes(fdr, sensitivity), alpha = 0.2)
p = p + geom_point()
p = p + geom_path(size = 1.5)
p = p + xlim(0, 0.15)
p = p + xlab('true false discovery rate')
p = p + scale_color_manual(values = cbb)
p = p + facet_wrap(~n_cells+effect_size_type, labeller = label_both)
p1 = p
save_plot('img/waterbear_unequal_sensitivity_p0.10.pdf', p1, base_height = 7)

# waterbear only

metadata1000 = dplyr::filter(metadata, n_mcmc_samples == 20000, n_total_genes == 1000,
  p_gene_effects == '0.10', n_cells == '2000000')

res = lapply(1:nrow(metadata1000),
  function(i) {
    print(i)
    info = metadata1000[i, ]
    read_nimble_results(info)
  })
res = bind_rows(res)
res = mutate(res, method = 'waterbear')

p = ggplot(
  dplyr::filter(res, effect_size_type == 'moderate'),
  aes(fdr, sensitivity, color = factor(bin_configuration)), alpha = 0.2)
# p = p + geom_point()
p = p + geom_path(size = 1.5)
p = p + xlim(0, 0.15)
p = p + xlab('true false discovery rate')
p = p + scale_color_manual(values = cbb)
p1 = p
save_plot('img/waterbear_unequal_sensitivity_p0.10_moderate.pdf', p1, base_height = 7)

bins_1 = list(
  list(
    c(0, 0.10), c(0.10, 0.5),
    c(0.5, 0.75), c(0.75, 1)),
  list(
    c(0, 0.2), c(0.2, 0.5),
    c(0.5, 0.8), c(0.8, 1)),
  list(
    c(0, 0.25), c(0.25, 0.5),
    c(0.5, 0.90), c(0.90, 1))
)

bin_sizes_1 = lapply(seq_along(bins_1),
  function(j) {
    r = bins_1[[j]]
    rep1 = lapply(seq_along(r),
      function(i) {
        data.frame(bin = i, mass = diff(r[[i]]), replicate = j, configuration = 1)
      })
    rep1 = bind_rows(rep1)
})
bin_sizes_1 = bind_rows(bin_sizes_1)
p = ggplot(bin_sizes_1, aes(x = replicate, y = mass, fill = factor(bin)))
p = p + geom_bar(stat = 'identity')
p = p + scale_fill_manual(values = cbb)
save_plot('bin_sizes_configuration1.pdf', p)

bins_2 = list(
  list(
    c(0, 0.25), c(0.25, 0.5),
    c(0.5, 0.90), c(0.90, 1)),
  list(
    c(0, 0.25), c(0.25, 0.5),
    c(0.5, 0.90), c(0.90, 1)),
  list(
    c(0, 0.25), c(0.25, 0.5),
    c(0.5, 0.90), c(0.90, 1))
)

bin_sizes_2 = lapply(seq_along(bins_2),
  function(j) {
    r = bins_2[[j]]
    rep2 = lapply(seq_along(r),
      function(i) {
        data.frame(bin = i, mass = diff(r[[i]]), replicate = j, configuration = 2)
      })
    rep2 = bind_rows(rep2)
})
bin_sizes_2 = bind_rows(bin_sizes_2)

p = ggplot(bin_sizes_2, aes(x = replicate, y = mass, fill = factor(bin)))
p = p + geom_bar(stat = 'identity')
p = p + scale_fill_manual(values = cbb)
save_plot('bin_sizes_configuration2.pdf', p)

bins_3 = list(
  list(
    c(0, 0.25), c(0.25, 0.5),
    c(0.5, 0.75), c(0.75, 1)),
  list(
    c(0, 0.25), c(0.25, 0.5),
    c(0.5, 0.75), c(0.75, 1)),
  list(
    c(0, 0.25), c(0.25, 0.5),
    c(0.5, 0.75), c(0.75, 1))
)

bin_sizes_3 = lapply(seq_along(bins_3),
  function(j) {
    r = bins_3[[j]]
    rep3 = lapply(seq_along(r),
      function(i) {
        data.frame(bin = i, mass = diff(r[[i]]), replicate = j, configuration = 3)
      })
    rep3 = bind_rows(rep3)
})
bin_sizes_3 = bind_rows(bin_sizes_3)
p = ggplot(bin_sizes_3, aes(x = replicate, y = mass, fill = factor(bin)))
p = p + geom_bar(stat = 'identity')
p = p + scale_fill_manual(values = cbb)
save_plot('bin_sizes_configuration3.pdf', p)

all_bin_sizes = bind_rows(bin_sizes_1, bin_sizes_2, bin_sizes_3)

p = ggplot(all_bin_sizes, aes(x = replicate, y = mass, fill = factor(bin)))
p = p + geom_bar(stat = 'identity')
p = p + scale_fill_manual(values = cbb)
p = p + facet_wrap(~configuration, nrow = 3, labeller = label_both)
save_plot('bin_sizes_configurations.pdf', p, base_height = 10, base_width = 4)

pdf('bin_cut_1.pdf')
x = seq(-3, 3, length.out = 1000)
plot(x, dnorm(x))
rep1 = sapply(bins_1[[1]],
  function(x) {
    qnorm(x[2])
  })
rep1 = rep1[1:3]
abline(h = rep, col = 'red')
dev.off()

# end waterbear only

filter(metadata1000, effect_size_type == 'moderate', p_gene_effects == '0.10')

all_results = dplyr::group_by(all_results, method, seed, effect_size_type, n_cells)
alpha = 0.10

summary_statistics = summarize(all_results,
  sensitivity = sum((value < alpha) & true_effect) / sum(true_effect),
  fpr = sum((value < alpha) & !true_effect) / sum(!true_effect),
  fdr = sum((value < alpha) & !true_effect) / sum(value < alpha))
arrange(summary_statistics, effect_size_type, n_cells, method)

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

p = ggplot(all_results, aes(fdr, sensitivity, color = method), alpha = 0.2)
p = p + geom_point()
p = p + geom_path(size = 1.5)
p = p + xlim(0, 0.15)
p = p + xlab('true false discovery rate')
p = p + scale_color_manual(values = cbb)
p = p + facet_wrap(~n_cells+effect_size_type, labeller = label_both)
p1 = p

save_plot('img/unequal_sensitivity_p0.50.pdf', p1, base_height = 7)

all_results = dplyr::group_by(all_results, method, seed, effect_size_type, n_cells)
alpha = 0.10

summary_statistics = summarize(all_results,
  sensitivity = sum((value < alpha) & true_effect) / sum(true_effect),
  fpr = sum((value < alpha) & !true_effect) / sum(!true_effect),
  fdr = sum((value < alpha) & !true_effect) / sum(value < alpha))
arrange(summary_statistics, effect_size_type, n_cells, method)


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

all_results = dplyr::group_by(all_results, method, seed, effect_size_type, n_cells)

alpha = 0.10

summary_statistics = summarize(all_results,
  sensitivity = sum((value < alpha) & true_effect) / sum(true_effect),
  fpr = sum((value < alpha) & !true_effect) / sum(!true_effect),
  fdr = sum((value < alpha) & !true_effect) / sum(value < alpha),
  tnr = sum(value > alpha & !true_effect) / sum(!true_effect))

dplyr::filter(summary_statistics, effect_size_type == 'easy')

p = ggplot(summary_statistics, aes(fdr, sensitivity, color = method, shape = n_cells))
p = p + geom_point(alpha = 0.5, size = 2)
p = p + facet_wrap(~effect_size_type + n_cells)

save_plot(paste0('img/magic_comparison', img_name),
  p, base_height = 7)


# equalmixture


process_subgroup = function(metadata, p_gene_effects, n_total_genes, bin_configuration,
  n_replicates, debug = FALSE, waterbear_only = FALSE) {

  cur_metadata = dplyr::filter(metadata,
    n_total_genes == n_total_genes,
    bin_configuration == bin_configuration,
    p_gene_effects == as.character(p_gene_effects))
  sim = readRDS(paste0('results/truth/nimble_', basename(cur_metadata$filename)[1]))

  oracle_mapping = generate_oracle_mapping(sim$const$guide_to_gene,
    p_gene_effects, n_total_genes)


  p_gene_effects = as.character(p_gene_effects)
  n_total_genes = as.character(n_total_genes)
  bin_configuration = as.character(bin_configuration)
  n_replicates = as.character(n_replicates)

  tmp_filter = data.frame(
    p_gene_effects = p_gene_effects,
    n_total_genes = n_total_genes,
    bin_configuration = bin_configuration,
    n_replicates = n_replicates
  )

  metadata_filter = dplyr::semi_join(cur_metadata, tmp_filter)

  res = lapply(1:nrow(metadata_filter),
    function(i) {
      print(i)
      info = metadata_filter[i, ]

      # TODO: classify effect sizes
      # browser()
      om = get_effect_size_information(info, oracle_mapping)
      read_waterbear_results(info, om)
    })
  res = bind_rows(res)
  # res = mutate(res, method = 'waterbear')
  # res = mutate(res,
  #   coverage = as.numeric(n_cells) / (
  #     as.integer(n_total_genes) * as.integer(n_guides_per_target) + as.integer(n_control_guides)))
  n_control_guides = metadata_filter$n_control_guides
  n_guides_per_target = metadata_filter$n_guides_per_target

  s_res = group_by(res, effect_size_type, n_cells, rank)
  s_res = summarize(s_res, mean_fdr = mean(fdr), sd_fdr = sd(fdr),
    mean_sensitivity = mean(sensitivity), sd_sensitivity = sd(sensitivity),
    mean_estimated_fdr = mean(value), sd_estimated_fdr = sd(value))
  s_res = mutate(s_res,
    coverage = as.numeric(n_cells) / (
      as.integer(n_total_genes[1]) * as.integer(n_guides_per_target[1]) + as.integer(n_control_guides[1])))
  s_res = mutate(s_res, method = 'waterbear')

  if (!waterbear_only) {
    mageck_results = lapply(1:nrow(metadata_filter),
      function(i) {
        print(i)
        info = metadata_filter[i, ]
        om = get_effect_size_information(info, oracle_mapping)
        process_mageck_results(info, om)
      })
    mageck_results = bind_rows(mageck_results)

    s_mageck_res = group_by(mageck_results, effect_size_type, n_cells, rank)
    s_mageck_res = summarize(s_mageck_res, mean_fdr = mean(fdr), sd_fdr = sd(fdr),
      mean_sensitivity = mean(sensitivity), sd_sensitivity = sd(sensitivity),
      mean_estimated_fdr = mean(value), sd_estimated_fdr = sd(value))
    s_mageck_res = mutate(s_mageck_res,
      coverage = as.numeric(n_cells) / (
        as.integer(n_total_genes) * as.integer(n_guides_per_target) + as.integer(n_control_guides)))
    s_mageck_res = mutate(s_mageck_res, method = 'MAGeCK')

    # maude_results = lapply(1:nrow(metadata_filter),
    #   function(i) {
    #     print(i)
    #     info = metadata_filter[i, ]
    #     om = get_effect_size_information(info, oracle_mapping)
    #     process_maude_results(info, om)
    #   })
    # maude_results = bind_rows(maude_results)

    # s_maude_res = group_by(maude_results, effect_size_type, n_cells, rank)
    # s_maude_res = summarize(s_maude_res, mean_fdr = mean(fdr), sd_fdr = sd(fdr),
    #   mean_sensitivity = mean(sensitivity), sd_sensitivity = sd(sensitivity),
    #   mean_estimated_fdr = mean(value), sd_estimated_fdr = sd(value))
    # s_maude_res = mutate(s_maude_res,
    #   coverage = as.numeric(n_cells) / (
    #     as.integer(n_total_genes) * as.integer(n_guides_per_target) + as.integer(n_control_guides)))
    # s_maude_res = mutate(s_maude_res, method = 'MAUDE')

    # all_res = bind_rows(s_res, s_mageck_res, s_maude_res)
    all_res = bind_rows(s_res, s_mageck_res)
    if (debug) {
      # return(list(all_res = all_res, waterbear = res, maude = maude_results,
      #     mageck = mageck_results))
      return(list(all_res = all_res, waterbear = res, mageck = mageck_results))
    }
  } else {
    all_res = bind_rows(s_res)
  }

  all_res
}

process_subgroup_by_method = function(metadata, p_gene_effects, n_total_genes, bin_configuration,
  n_replicates, method = 'waterbear', summary = TRUE, debug = FALSE) {
  cur_metadata = dplyr::filter(metadata,
    n_total_genes == n_total_genes,
    bin_configuration == bin_configuration,
    p_gene_effects == as.character(p_gene_effects))
  sim = readRDS(paste0('results/truth/nimble_', basename(cur_metadata$filename)[1]))

  oracle_mapping = generate_oracle_mapping(sim$const$guide_to_gene,
    p_gene_effects, n_total_genes)

  p_gene_effects = as.character(p_gene_effects)
  n_total_genes = as.character(n_total_genes)
  bin_configuration = as.character(bin_configuration)
  n_replicates = as.character(n_replicates)

  tmp_filter = data.frame(
    p_gene_effects = p_gene_effects,
    n_total_genes = n_total_genes,
    bin_configuration = bin_configuration,
    n_replicates = n_replicates
  )

  metadata_filter = dplyr::semi_join(cur_metadata, tmp_filter)
  # WARNING: assumes the same control guides everywhere
  # n_control_guides = metadata_filter$n_control_guides[1]
  # n_guides_per_target = metadata_filter$n_guides_per_target

  all_res = list()
  if (method == 'waterbear') {
    print('waterbear')
    res = lapply(1:nrow(metadata_filter),
      function(i) {
        print(i)
        info = metadata_filter[i, ]
        om = get_effect_size_information(info, oracle_mapping)
        res = read_waterbear_results(info, om)
        mutate(res, n_control_guides = info$n_control_guides,
          n_guides_per_target = info$n_guides_per_target)
      })
    res = bind_rows(res)
    if (summary) {
      s_res = group_by(res, n_control_guides, effect_size_type, n_cells, n_guides_per_target, rank)
      s_res = summarize(s_res, mean_fdr = mean(fdr), sd_fdr = sd(fdr),
        mean_sensitivity = mean(sensitivity), sd_sensitivity = sd(sensitivity),
        mean_estimated_fdr = mean(value), sd_estimated_fdr = sd(value))
      s_res = mutate(s_res,
        coverage = as.numeric(n_cells) / (
          as.integer(n_total_genes) * as.integer(n_guides_per_target) + as.integer(n_control_guides)))
      s_res = mutate(s_res, method = 'waterbear')
    } else {
      s_res = mutate(res, method = 'waterbear')
    }
    all_res$waterbear = s_res
  } else if (method == 'mageck') {
    print('mageck')
    mageck_results = lapply(1:nrow(metadata_filter),
      function(i) {
        print(i)
        info = metadata_filter[i, ]
        om = get_effect_size_information(info, oracle_mapping)
        res = process_mageck_results(info, om)
        mutate(res, n_control_guides = info$n_control_guides,
          n_guides_per_target = info$n_guides_per_target)
      })
    mageck_results = bind_rows(mageck_results)

    if (summary) {
      s_mageck_res = group_by(mageck_results, n_control_guides, effect_size_type, n_cells, n_guides_per_target, rank)
      s_mageck_res = summarize(s_mageck_res, mean_fdr = mean(fdr), sd_fdr = sd(fdr),
        mean_sensitivity = mean(sensitivity), sd_sensitivity = sd(sensitivity),
        mean_estimated_fdr = mean(value), sd_estimated_fdr = sd(value))
      s_mageck_res = mutate(s_mageck_res,
        coverage = as.numeric(n_cells) / (
          as.integer(n_total_genes) * as.integer(n_guides_per_target) + as.integer(n_control_guides)))
      s_mageck_res = mutate(s_mageck_res, method = 'MAGeCK')
    } else {
      s_mageck_res = mutate(mageck_results, method = 'MAGeCK')
    }
    all_res$mageck = s_mageck_res
  } else if (method == 'maude') {
    maude_results = lapply(1:nrow(metadata_filter),
      function(i) {
        print(i)
        info = metadata_filter[i, ]
        om = get_effect_size_information(info, oracle_mapping)
        res = process_maude_results(info, om)
        mutate(res, n_control_guides = info$n_control_guides,
          n_guides_per_target = info$n_guides_per_target)
      })
    maude_results = bind_rows(maude_results)

    if (summary) {
      s_maude_res = group_by(maude_results, n_control_guides, effect_size_type, n_cells, n_guides_per_target, rank)
      s_maude_res = summarize(s_maude_res, mean_fdr = mean(fdr), sd_fdr = sd(fdr),
        mean_sensitivity = mean(sensitivity), sd_sensitivity = sd(sensitivity),
        mean_estimated_fdr = mean(value), sd_estimated_fdr = sd(value))
      s_maude_res = mutate(s_maude_res,
        coverage = as.numeric(n_cells) / (
          as.integer(n_total_genes) * as.integer(n_guides_per_target) + as.integer(n_control_guides)))
      s_maude_res = mutate(s_maude_res, method = 'MAUDE')
    } else {
      s_maude_res = mutate(maude_results, method = 'MAUDE')
    }
    all_res$maude = s_maude_res
  }
  all_res
}

source('benchmark_common.R')
source('common.R')

fnames = Sys.glob('results/nimble/unequal_*equalmixture*.rds')

metadata = parse_unequal_bins_fname(basename(fnames))
metadata = mutate(metadata, filename = fnames)
metadata = filter(metadata, !grepl('summary', filename))

metadata = dplyr::mutate(metadata,
  summary_filename = paste0(dirname(filename), '/', sub('.rds', '_summary.rds', basename(filename))))
metadata = dplyr::filter(metadata, file.exists(summary_filename))

metadata = dplyr::mutate(metadata,
  posterior_summary = paste0(dirname(filename), '/', sub('.rds', '_posterior_summary.rds', basename(filename))))
metadata = dplyr::filter(metadata, file.exists(posterior_summary))

metadata = dplyr::mutate(metadata,
  maude_filename = paste0(sub('nimble', 'maude', dirname(filename)), '/',
    sub('.rds', '.csv', basename(filename))))
metadata = filter(metadata, file.exists(maude_filename))

p_gene_effects = '0.10'
n_total_genes = '1000'
bin_configuration = '1'
n_replicates = '3'


p_10_maude = process_subgroup_by_method(metadata, p_gene_effects, n_total_genes, bin_configuration, n_replicates,
  method = 'maude',
  summary = FALSE,
  debug = TRUE)

p_10_waterbear = process_subgroup_by_method(metadata, p_gene_effects, n_total_genes, bin_configuration, n_replicates,
  method = 'waterbear',
  summary = FALSE,
  debug = TRUE)

p_10_mageck = process_subgroup_by_method(metadata, p_gene_effects, n_total_genes, bin_configuration, n_replicates,
  method = 'mageck',
  summary = FALSE,
  debug = TRUE)

all_res = bind_rows(p_10_maude$maude, p_10_waterbear$waterbear, p_10_mageck$mageck)

alpha = 0.10
tmp = group_by(all_res,
  method, n_control_guides, seed, effect_size_type, n_cells, n_guides_per_target)
get_est_fdr = function(df, alpha, eps = 0.05) {
  df = mutate(df, d = abs(value - alpha))
  min_distance = min(df$d)
  fdr_summary = filter(df, d == min_distance)
  if (nrow(fdr_summary) > 1) {
    fdr_summary = filter(fdr_summary, rank == max(rank))
  } else if (nrow(fdr_summary) == 0) {
    warning('ERROR: no valid fdr')
    fdr_summary = df
  }
  fdr_summary
}
tmp = group_map(tmp, get_est_fdr, alpha = 0.10, keep = TRUE)
tmp = bind_rows(tmp)
tmp = mutate(tmp, coverage = as.integer(n_cells) /
  (as.integer(n_total_genes) * as.integer(n_guides_per_target) + as.integer(n_control_guides)))

tmp = filter(tmp, n_guides_per_target == '4', coverage <= 1000)

z_95 = abs(qnorm(0.025))
tmp_summary = group_by(tmp, method, coverage) %>%
  summarize(
    m_sensitivity = mean(sensitivity),
    m_sensitivity_sd = sd(sensitivity),
    m_sensitivity_sem = m_sensitivity_sd / sqrt(length(sensitivity)),
    m_95_sensitivity = z_95 * m_sensitivity_sem,
    m_fdr = mean(fdr),
    m_fdr_sd = sd(fdr),
    m_fdr_sem = m_fdr_sd / sqrt(length(fdr)),
    m_95_fdr = z_95 * m_fdr_sem
  )
tmp_summary = filter(tmp_summary, coverage <= 1000)

tmp_summary = mutate(tmp_summary, f_coverage = factor(coverage), n_coverage = as.numeric(f_coverage))
color_mapping = distinct(
  select(ungroup(tmp_summary), f_coverage))
color_mapping = mutate(color_mapping, coverage_color = rep(c('gray', 'white'), length.out = nrow(color_mapping)))
tmp_summary = left_join(tmp_summary, color_mapping, by = 'f_coverage')


# these are the current plots
x_dodge = 0.5
p = ggplot(tmp_summary,
  aes(f_coverage, m_sensitivity, color = method))
p = p + geom_point(position = position_dodge2(width = x_dodge, padding = x_dodge), alpha = 0)
p = p + geom_rect(
  aes(xmin = n_coverage - x_dodge, xmax = n_coverage + x_dodge, ymin = 0, ymax = 1, fill = f_coverage), color = NA, alpha = 0.05)
p = p + scale_fill_manual(values = color_mapping$coverage_color)
p = p + geom_point(position = position_dodge2(width = x_dodge, padding = x_dodge))
p = p + geom_pointrange(aes(ymin = m_sensitivity - m_95_sensitivity,
    ymax = m_sensitivity + m_95_sensitivity),
  position = position_dodge2(width = x_dodge, padding = x_dodge))
p = p + scale_color_brewer(palette = 'Dark2')
p = p + xlab('coverage')
p = p + ylab('sensitivity')
p = p + ylim(0, 1)
save_plot('img/unequal_average_sensitivity_compare.pdf', p)

# these are the current plots
x_dodge = 0.0
p = ggplot(tmp_summary,
  aes(method, m_sensitivity, color = method))
# p = p + geom_point(position = position_dodge2(width = 0, padding = x_dodge), alpha = 0)
# p = p + geom_rect(
#   aes(xmin = n_coverage - x_dodge, xmax = n_coverage + x_dodge, ymin = 0, ymax = 1), color = NA, alpha = 0.05)
# p = p + scale_fill_manual(values = color_mapping$coverage_color)
p = p + geom_point()
p = p + geom_pointrange(aes(ymin = m_sensitivity - m_95_sensitivity,
    ymax = m_sensitivity + m_95_sensitivity),
  )
p = p + scale_x_discrete(expand=c(1, 1))
p = p + scale_color_brewer(palette = 'Dark2')
p = p + xlab('coverage')
p = p + ylab('sensitivity')
p = p + ylim(0, 1)
p = p + facet_wrap(~f_coverage, nrow = 1)
p = p + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(),
  panel.spacing = unit(0.1, "lines"),
  panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.1)
)
save_plot('img/grid_unequal_average_sensitivity_compare.pdf', p)

p = ggplot(tmp_summary,
  aes(f_coverage, m_fdr, color = method))
p = p + geom_point(position = position_dodge2(width = x_dodge, padding = x_dodge), alpha = 0)
p = p + geom_rect(
  aes(xmin = n_coverage - x_dodge, xmax = n_coverage + x_dodge, ymin = 0, ymax = 1, fill = f_coverage), color = NA, alpha = 0.05)
p = p + scale_fill_manual(values = color_mapping$coverage_color)
p = p + geom_point(position = position_dodge2(width = 0.5, padding = 0.5))
p = p + geom_pointrange(aes(ymin = m_fdr - m_95_fdr,
    ymax = m_fdr + m_95_fdr),
  position = position_dodge2(width = 0.5, padding = 0.5))
p = p + scale_color_brewer(palette = 'Dark2')
p = p + xlab('coverage')
p = p + ylab('fdr')
p = p + ylim(0, 1)
p = p + geom_hline(yintercept = alpha, linetype = 3)
save_plot('img/unequal_average_fdr_compare.pdf', p)

p = ggplot(tmp_summary,
  aes(method, m_fdr, color = method))
p = p + geom_point()
p = p + geom_pointrange(aes(ymin = m_fdr - m_95_fdr,
    ymax = m_fdr + m_95_fdr))
p = p + scale_x_discrete(expand=c(1, 1))
p = p + scale_color_brewer(palette = 'Dark2')
p = p + xlab('coverage')
p = p + ylab('fdr')
p = p + ylim(0, 1)
p = p + geom_hline(yintercept = alpha, linetype = 3)
p = p + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(),
  panel.spacing = unit(0.1, "lines"),
  panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.1)
)
p = p + facet_wrap(~f_coverage, nrow = 1)
save_plot('img/grid_unequal_average_fdr_compare.pdf', p)
# bin 2

p_gene_effects = '0.10'
n_total_genes = '1000'
bin_configuration = '2'
n_replicates = '3'


p_10_maude = process_subgroup_by_method(metadata, p_gene_effects, n_total_genes, bin_configuration, n_replicates,
  method = 'maude',
  summary = FALSE,
  debug = TRUE)

p_10_waterbear = process_subgroup_by_method(metadata, p_gene_effects, n_total_genes, bin_configuration, n_replicates,
  method = 'waterbear',
  summary = FALSE,
  debug = TRUE)

p_10_mageck = process_subgroup_by_method(metadata, p_gene_effects, n_total_genes, bin_configuration, n_replicates,
  method = 'mageck',
  summary = FALSE,
  debug = TRUE)

all_res = bind_rows(p_10_maude$maude, p_10_waterbear$waterbear, p_10_mageck$mageck)

alpha = 0.10
tmp = group_by(all_res,
  method, n_control_guides, seed, effect_size_type, n_cells, n_guides_per_target)
get_est_fdr = function(df, alpha, eps = 0.05) {
  df = mutate(df, d = abs(value - alpha))
  min_distance = min(df$d)
  fdr_summary = filter(df, d == min_distance)
  if (nrow(fdr_summary) > 1) {
    fdr_summary = filter(fdr_summary, rank == max(rank))
  } else if (nrow(fdr_summary) == 0) {
    warning('ERROR: no valid fdr')
    fdr_summary = df
  }
  fdr_summary
}
tmp = group_map(tmp, get_est_fdr, alpha = 0.10, keep = TRUE)
tmp = bind_rows(tmp)
tmp = mutate(tmp, coverage = as.integer(n_cells) /
  (as.integer(n_total_genes) * as.integer(n_guides_per_target) + as.integer(n_control_guides)))

tmp = filter(tmp, n_guides_per_target == '4', coverage <= 1000)

z_95 = abs(qnorm(0.025))
tmp_summary = group_by(tmp, method, coverage) %>%
  summarize(
    m_sensitivity = mean(sensitivity),
    m_sensitivity_sd = sd(sensitivity),
    m_sensitivity_sem = m_sensitivity_sd / sqrt(length(sensitivity)),
    m_95_sensitivity = z_95 * m_sensitivity_sem,
    m_fdr = mean(fdr),
    m_fdr_sd = sd(fdr),
    m_fdr_sem = m_fdr_sd / sqrt(length(fdr)),
    m_95_fdr = z_95 * m_fdr_sem
  )
tmp_summary = filter(tmp_summary, coverage <= 1000)

tmp_summary = mutate(tmp_summary, f_coverage = factor(coverage), n_coverage = as.numeric(f_coverage))
color_mapping = distinct(
  select(ungroup(tmp_summary), f_coverage))
color_mapping = mutate(color_mapping, coverage_color = rep(c('gray', 'white'), length.out = nrow(color_mapping)))
tmp_summary = left_join(tmp_summary, color_mapping, by = 'f_coverage')


# these are the current plots
x_dodge = 0.5
p = ggplot(tmp_summary,
  aes(f_coverage, m_sensitivity, color = method))
p = p + geom_point(position = position_dodge2(width = x_dodge, padding = x_dodge), alpha = 0)
p = p + geom_rect(
  aes(xmin = n_coverage - x_dodge, xmax = n_coverage + x_dodge, ymin = 0, ymax = 1, fill = f_coverage), color = NA, alpha = 0.05)
p = p + scale_fill_manual(values = color_mapping$coverage_color)
p = p + geom_point(position = position_dodge2(width = x_dodge, padding = x_dodge))
p = p + geom_pointrange(aes(ymin = m_sensitivity - m_95_sensitivity,
    ymax = m_sensitivity + m_95_sensitivity),
  position = position_dodge2(width = x_dodge, padding = x_dodge))
p = p + scale_color_brewer(palette = 'Dark2')
p = p + xlab('coverage')
p = p + ylab('sensitivity')
p = p + ylim(0, 1)
save_plot('img/unequal2_average_sensitivity_compare.pdf', p)

p = ggplot(tmp_summary,
  aes(f_coverage, m_fdr, color = method))
p = p + geom_point(position = position_dodge2(width = x_dodge, padding = x_dodge), alpha = 0)
p = p + geom_rect(
  aes(xmin = n_coverage - x_dodge, xmax = n_coverage + x_dodge, ymin = 0, ymax = 1, fill = f_coverage), color = NA, alpha = 0.05)
p = p + scale_fill_manual(values = color_mapping$coverage_color)
p = p + geom_point(position = position_dodge2(width = 0.5, padding = 0.5))
p = p + geom_pointrange(aes(ymin = m_fdr - m_95_fdr,
    ymax = m_fdr + m_95_fdr),
  position = position_dodge2(width = 0.5, padding = 0.5))
p = p + scale_color_brewer(palette = 'Dark2')
p = p + xlab('coverage')
p = p + ylab('fdr')
p = p + ylim(0, 1)
p = p + geom_hline(yintercept = alpha, linetype = 3)
save_plot('img/unequal2_average_fdr_compare.pdf', p)

# bin 3

p_gene_effects = '0.10'
n_total_genes = '1000'
bin_configuration = '3'
n_replicates = '3'


p_10_maude = process_subgroup_by_method(metadata, p_gene_effects, n_total_genes, bin_configuration, n_replicates,
  method = 'maude',
  summary = FALSE,
  debug = TRUE)

p_10_waterbear = process_subgroup_by_method(metadata, p_gene_effects, n_total_genes, bin_configuration, n_replicates,
  method = 'waterbear',
  summary = FALSE,
  debug = TRUE)

p_10_mageck = process_subgroup_by_method(metadata, p_gene_effects, n_total_genes, bin_configuration, n_replicates,
  method = 'mageck',
  summary = FALSE,
  debug = TRUE)

all_res = bind_rows(p_10_maude$maude, p_10_waterbear$waterbear, p_10_mageck$mageck)

alpha = 0.10
tmp = group_by(all_res,
  method, n_control_guides, seed, effect_size_type, n_cells, n_guides_per_target)
get_est_fdr = function(df, alpha, eps = 0.05) {
  df = mutate(df, d = abs(value - alpha))
  min_distance = min(df$d)
  fdr_summary = filter(df, d == min_distance)
  if (nrow(fdr_summary) > 1) {
    fdr_summary = filter(fdr_summary, rank == max(rank))
  } else if (nrow(fdr_summary) == 0) {
    warning('ERROR: no valid fdr')
    fdr_summary = df
  }
  fdr_summary
}
tmp = group_map(tmp, get_est_fdr, alpha = 0.10, keep = TRUE)
tmp = bind_rows(tmp)
tmp = mutate(tmp, coverage = as.integer(n_cells) /
  (as.integer(n_total_genes) * as.integer(n_guides_per_target) + as.integer(n_control_guides)))

tmp = filter(tmp, n_guides_per_target == '4', coverage <= 1000)

z_95 = abs(qnorm(0.025))
tmp_summary = group_by(tmp, method, coverage) %>%
  summarize(
    m_sensitivity = mean(sensitivity),
    m_sensitivity_sd = sd(sensitivity),
    m_sensitivity_sem = m_sensitivity_sd / sqrt(length(sensitivity)),
    m_95_sensitivity = z_95 * m_sensitivity_sem,
    m_fdr = mean(fdr),
    m_fdr_sd = sd(fdr),
    m_fdr_sem = m_fdr_sd / sqrt(length(fdr)),
    m_95_fdr = z_95 * m_fdr_sem
  )
tmp_summary = filter(tmp_summary, coverage <= 1000)

tmp_summary = mutate(tmp_summary, f_coverage = factor(coverage), n_coverage = as.numeric(f_coverage))
color_mapping = distinct(
  select(ungroup(tmp_summary), f_coverage))
color_mapping = mutate(color_mapping, coverage_color = rep(c('gray', 'white'), length.out = nrow(color_mapping)))
tmp_summary = left_join(tmp_summary, color_mapping, by = 'f_coverage')


# these are the current plots
x_dodge = 0.5
p = ggplot(tmp_summary,
  aes(f_coverage, m_sensitivity, color = method))
p = p + geom_point(position = position_dodge2(width = x_dodge, padding = x_dodge), alpha = 0)
p = p + geom_rect(
  aes(xmin = n_coverage - x_dodge, xmax = n_coverage + x_dodge, ymin = 0, ymax = 1, fill = f_coverage), color = NA, alpha = 0.05)
p = p + scale_fill_manual(values = color_mapping$coverage_color)
p = p + geom_point(position = position_dodge2(width = x_dodge, padding = x_dodge))
p = p + geom_pointrange(aes(ymin = m_sensitivity - m_95_sensitivity,
    ymax = m_sensitivity + m_95_sensitivity),
  position = position_dodge2(width = x_dodge, padding = x_dodge))
p = p + scale_color_brewer(palette = 'Dark2')
p = p + xlab('coverage')
p = p + ylab('sensitivity')
p = p + ylim(0, 1)
save_plot('img/unequal3_average_sensitivity_compare.pdf', p)

p = ggplot(tmp_summary,
  aes(f_coverage, m_fdr, color = method))
p = p + geom_point(position = position_dodge2(width = x_dodge, padding = x_dodge), alpha = 0)
p = p + geom_rect(
  aes(xmin = n_coverage - x_dodge, xmax = n_coverage + x_dodge, ymin = 0, ymax = 1, fill = f_coverage), color = NA, alpha = 0.05)
p = p + scale_fill_manual(values = color_mapping$coverage_color)
p = p + geom_point(position = position_dodge2(width = 0.5, padding = 0.5))
p = p + geom_pointrange(aes(ymin = m_fdr - m_95_fdr,
    ymax = m_fdr + m_95_fdr),
  position = position_dodge2(width = 0.5, padding = 0.5))
p = p + scale_color_brewer(palette = 'Dark2')
p = p + xlab('coverage')
p = p + ylab('fdr')
p = p + ylim(0, 1)
p = p + geom_hline(yintercept = alpha, linetype = 3)
save_plot('img/unequal3_average_fdr_compare.pdf', p)



p_10 = process_subgroup(metadata,
  p_gene_effects, n_total_genes, bin_configuration, n_replicates, waterbear_only = FALSE)


debug_p_10 = process_subgroup(metadata,
  p_gene_effects, n_total_genes, bin_configuration, n_replicates, waterbear_only = FALSE,
  debug = TRUE)

all_methods = bind_rows(
  dplyr::mutate(debug_p_10$waterbear, method = 'waterbear'),
  # dplyr::mutate(p_10$maude, method = 'MAUDE'),
  dplyr::mutate(debug_p_10$mageck, method = 'MAGeCK'))

all_methods = group_by(all_methods, method, n_cells, seed)
all_methods = dplyr::mutate(all_methods, coverage = as.integer(n_cells) / 4100)
epsilon = 0.001




efdr_10 = dplyr::filter(group_by(all_methods, method, n_cells, seed), value < 0.10)
efdr_10 = dplyr::filter(efdr_10, value == max(value))

filter(efdr_10, method == 'waterbear', coverage == 50)

# get the largest sensitivity value
fdr_10 = dplyr::filter(all_methods, fdr < 0.10 + epsilon, 0.10 - epsilon < fdr)


# TODO: DEBUG what is happening at the low coverage situations
tmp = dplyr::filter(all_methods, coverage == 50, method == 'waterbear')
tmp

p = ggplot(dplyr::filter(efdr_10, coverage <= 5000),
  aes(as.factor(coverage), sensitivity, fill = method))
p = p + geom_boxplot()
# p = p + facet_wrap(~n_cells)
save_plot('img/unequal_bins_all_methods_fdr_10.pdf', p)

# bin configuration 3

p_gene_effects = '0.10'
n_total_genes = '1000'
bin_configuration = '3'
n_replicates = '3'


p_10 = process_subgroup(metadata,
  p_gene_effects, n_total_genes, bin_configuration, n_replicates, waterbear_only = FALSE)

debug_p_10 = process_subgroup(metadata,
  p_gene_effects, n_total_genes, bin_configuration, n_replicates, waterbear_only = FALSE,
  debug = TRUE)

all_methods = bind_rows(
  dplyr::mutate(debug_p_10$waterbear, method = 'waterbear'),
  # dplyr::mutate(p_10$maude, method = 'MAUDE'),
  dplyr::mutate(debug_p_10$mageck, method = 'MAGeCK'))

all_methods = group_by(all_methods, method, n_cells, seed)
all_methods = dplyr::mutate(all_methods, coverage = as.integer(n_cells) / 4100)
epsilon = 0.001


efdr_10 = dplyr::filter(group_by(all_methods, method, n_cells, seed), value < 0.10)
efdr_10 = dplyr::filter(efdr_10, value == max(value))

p = ggplot(dplyr::filter(efdr_10, coverage <= 5000),
  aes(as.factor(coverage), sensitivity, fill = method))
p = p + geom_boxplot()
# p = p + facet_wrap(~n_cells)
save_plot('img/bin_3_unequal_bins_all_methods_fdr_10.pdf', p)

# bin configuration 4

p_gene_effects = '0.10'
n_total_genes = '1000'
bin_configuration = '4'
n_replicates = '3'


# debugonce(read_waterbear_results)
debug_p_10 = process_subgroup(metadata,
  p_gene_effects, n_total_genes, bin_configuration, n_replicates, waterbear_only = FALSE,
  debug = TRUE)

all_methods = bind_rows(
  dplyr::mutate(debug_p_10$waterbear, method = 'waterbear'),
  # dplyr::mutate(p_10$maude, method = 'MAUDE'),
  dplyr::mutate(debug_p_10$mageck, method = 'MAGeCK'))

all_methods = group_by(all_methods, method, n_cells, seed)
all_methods = dplyr::mutate(all_methods, coverage = as.integer(n_cells) / 4100)
epsilon = 0.001


efdr_10 = dplyr::filter(group_by(all_methods, method, n_cells, seed), value < 0.10)
efdr_10 = dplyr::filter(efdr_10, value == max(value))

p = ggplot(dplyr::filter(efdr_10, coverage <= 5000),
  aes(as.factor(coverage), sensitivity, fill = method))
p = p + geom_boxplot()
# p = p + facet_wrap(~n_cells)
save_plot('img/bin_4_unequal_bins_all_methods_fdr_10.pdf', p)

p = ggplot(dplyr::filter(efdr_10, coverage <= 5000),
  aes(as.factor(coverage), fdr, fill = method))
p = p + geom_boxplot()
# p = p + facet_wrap(~n_cells)
save_plot('img/bin_4_unequal_bins_all_methods_fdr_10_fdr.pdf', p)


# testing to verify FDR

dplyr::filter(metadata, n_cells == 205000, bin_configuration == 4)
