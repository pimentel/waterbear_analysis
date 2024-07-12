library('dplyr')
library('tidyr')
library('ggplot2')
library('cowplot')

cbb <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

theme_set(theme_cowplot())

source('benchmark_common.R')


########################################################################
# looking at coverage
fnames = Sys.glob('results/nimble/vary_coverage_*_20000_*.rds')

parse_coverage_fname = function(x) {
  x = sub('vary_coverage_', '', x)
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

metadata = parse_coverage_fname(basename(fnames))
metadata = mutate(metadata, filename = fnames)
metadata = dplyr::filter(metadata, grepl('vary_coverage_[0-9]+_[0-9]+_.*', filename))

metadata = filter(metadata, n_total_genes == 20000)

metadata = dplyr::filter(metadata, !grepl('_summary.rds', filename))
metadata = dplyr::mutate(metadata,
  summary_filename = paste0(dirname(filename), '/', sub('.rds', '_summary.rds', basename(filename))))
metadata = dplyr::filter(metadata, file.exists(summary_filename))

# metadata = dplyr::mutate(metadata,
#   maude_filename = paste0(sub('nimble', 'maude', dirname(filename)), '/',
#     sub('.rds', '.csv', basename(filename))))

metadata = dplyr::filter(metadata, n_mcmc_samples == 20000)

# load the truth
tmp_fname = sub('/nimble/', '/truth/sim_', metadata$filename)
tmp = readRDS(tmp_fname[3])

tmp_es = compute_gene_effect_sizes(tmp$true_parameters, 4, 100)
res = classify_effect_sizes(tmp_es$gene_effect_size, 'coverage_effect_sizes.RData')

tmp = classify_effect_sizes(tmp_es$effect_left, 'coverage_effect_sizes.RData')

x = abs(tmp_es$gene_effect_size[1])
res$lower <= x & x < res$upper

tmp_es$gene_effect_size


# end: load truth

nrow(dplyr::filter(metadata, !file.exists(summary_filename)))

tmp = dplyr::filter(metadata, !file.exists(maude_filename))
nrow(tmp)

# currently, the only thing that is changing is the p_gene_effects

process_subgroup = function(metadata, p_gene_effects, n_total_genes, bin_configuration,
  n_replicates, debug = FALSE) {

  base = paste0('vary_coverage_', p_gene_effects, '_', n_total_genes, '_',
    bin_configuration, '_', n_replicates)

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
      as.integer(n_total_genes) * as.integer(n_guides_per_target) + as.integer(n_control_guides)))
  s_res = mutate(s_res, method = 'waterbear')

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

  maude_results = lapply(1:nrow(metadata_filter),
    function(i) {
      print(i)
      info = metadata_filter[i, ]
      om = get_effect_size_information(info, oracle_mapping)
      process_maude_results(info, om)
    })
  maude_results = bind_rows(maude_results)

  s_maude_res = group_by(maude_results, effect_size_type, n_cells, rank)
  s_maude_res = summarize(s_maude_res, mean_fdr = mean(fdr), sd_fdr = sd(fdr),
    mean_sensitivity = mean(sensitivity), sd_sensitivity = sd(sensitivity),
    mean_estimated_fdr = mean(value), sd_estimated_fdr = sd(value))
  s_maude_res = mutate(s_maude_res,
    coverage = as.numeric(n_cells) / (
      as.integer(n_total_genes) * as.integer(n_guides_per_target) + as.integer(n_control_guides)))
  s_maude_res = mutate(s_maude_res, method = 'MAUDE')

  all_res = bind_rows(s_res, s_mageck_res, s_maude_res)
  if (debug) {
    return(list(all_res = all_res, waterbear = res, maude = maude_results,
        mageck = mageck_results))
  }
  all_res
}

get_estimated_fdr = function(s_res) {
  all_estimated_fdr = lapply(c(0.01, 0.05, 0.10),
    function(l) {
      e_fdr = filter(
        dplyr::filter(
          group_by(s_res, effect_size_type, coverage, method),
          mean_estimated_fdr <= l),
        rank == max(rank)
      )
      mutate(e_fdr, level = l)
    })
  all_estimated_fdr = bind_rows(all_estimated_fdr)
  all_estimated_fdr
}

p_gene_effects = '0.10'
n_total_genes = '1000'
bin_configuration = '1'
n_replicates = '3'

debugonce(compute_benchmark_statistics)

debugonce(get_effect_size_information)

p_10 = process_subgroup(metadata, p_gene_effects, n_total_genes, bin_configuration, n_replicates)

tmp = metadata %>%
  filter(p_gene_effects == '0.10', n_replicates == '3')

process_mageck_results(tmp, TRUE)

file.exists(nimble_fname_to_mageck(tmp$filename[1]))

p_10_efdr = get_estimated_fdr(p_10)

base = paste0('vary_coverage_', p_gene_effects, '_', n_total_genes, '_',
  bin_configuration, '_', n_replicates)


p = ggplot(
  p_10,
  aes(mean_fdr, mean_sensitivity, color = factor(coverage)), alpha = 0.2)
p = p + geom_path(size = 1.5)
p = p + geom_point(aes(shape = factor(level)), data = p_10_efdr, size = 5)
p = p + geom_errorbar(
  aes(
    ymin = mean_sensitivity - sd_sensitivity,
    ymax = mean_sensitivity + sd_sensitivity), data = p_10_efdr, size = 2)
p = p + geom_errorbarh(
  aes(
    xmin = mean_fdr - sd_fdr,
    xmax = mean_fdr + sd_fdr), data = p_10_efdr, size = 2)
# p = p + xlim(0, 0.25)
p = p + xlim(0, 0.15)
p = p + xlab('true false discovery rate')
p = p + scale_color_manual(values = cbb)
p = p + facet_wrap(~method + effect_size_type)
p1 = p

save_plot(paste0('img/summary_', base, '.pdf'), p1, base_height = 7)


tmp_10 = dplyr::filter(p_10, effect_size_type == 'equalmixture')
tmp_10_efdr = dplyr::filter(p_10_efdr, effect_size_type == 'equalmixture')

p = ggplot(
  tmp_10,
  aes(mean_fdr, mean_sensitivity, color = factor(method)), alpha = 0.2)
p = p + geom_path(size = 1.5)
p = p + geom_point(aes(shape = factor(level)), data = tmp_10_efdr, size = 5)
p = p + geom_errorbar(
  aes(
    ymin = mean_sensitivity - sd_sensitivity,
    ymax = mean_sensitivity + sd_sensitivity), data = tmp_10_efdr, size = 2)
p = p + geom_errorbarh(
  aes(
    xmin = mean_fdr - sd_fdr,
    xmax = mean_fdr + sd_fdr), data = tmp_10_efdr, size = 2)
# p = p + xlim(0, 0.25)
p = p + xlim(0, 0.15)
p = p + xlab('true false discovery rate')
p = p + scale_color_manual(values = cbb)
p = p + facet_wrap(~coverage, labeller = 'label_both')
p1 = p
save_plot(paste0('img/equal_mixture_summary_', base, '.pdf'), p1, base_height = 7)

p = ggplot(
  dplyr::filter(tmp_10, method == 'waterbear', coverage < 50000),
  aes(mean_fdr, mean_sensitivity, color = factor(coverage)), alpha = 0.2)
p = p + geom_path(size = 1.5)
# p = p + geom_point(aes(shape = factor(level)), data = dplyr::filter(tmp_10_efdr), size = 5)
# p = p + geom_errorbar(
#   aes(
#     ymin = mean_sensitivity - sd_sensitivity,
#     ymax = mean_sensitivity + sd_sensitivity), data = tmp_10_efdr, size = 2)
# p = p + geom_errorbarh(
#   aes(
#     xmin = mean_fdr - sd_fdr,
#     xmax = mean_fdr + sd_fdr), data = tmp_10_efdr, size = 2)
# p = p + xlim(0, 0.25)
p = p + xlim(0, 0.15)
p = p + xlab('true false discovery rate')
p = p + scale_color_manual(values = cbb)
p1 = p
save_plot(paste0('img/waterbear_equal_mixture_summary_', base, '.pdf'), p1, base_height = 7)


# TODO: attempting to see what is wrong with p = 0.1, difficult

metadata_em = dplyr::filter(metadata, effect_size_type == 'equalmixture')

p_10 = process_subgroup(metadata_em, p_gene_effects, n_total_genes, bin_configuration, n_replicates,
  TRUE)

wb_results = p_10$waterbear
wb_results = mutate(wb_results, pip = 1 - value)

pip_group = function(pip) {
  pip_mapping = vector('character', length(pip))
  pip_mapping[pip < 0.90] = 'none'
  pip_mapping[pip >= 0.99] = 'PIP >= 0.99'
  pip_mapping[pip >= 0.95 & pip < 0.99] = '0.95 <= PIP < 0.99'
  pip_mapping[pip >= 0.90 & pip < 0.95] = '0.90 <= PIP < 0.95'
  pip_mapping = factor(pip_mapping, c('none', '0.90 <= PIP < 0.95', '0.95 <= PIP < 0.99',
      'PIP >= 0.99')
  )
  pip_mapping
}

pip_summary = function(res) {
  res = mutate(res, pip = 1 - value)
  res = group_by(res, effect_size_type, n_cells, seed)
  res = mutate(res, pip_mapping = pip_group(pip))

  res = group_by(res, effect_size_type, n_cells, seed, pip_mapping)

  res = summarize(res, sensitivity = sum(true_effect) / length(true_effect),
    tp = sum(true_effect), fp = sum(!true_effect), fdr = fp / length(true_effect))
  res = dplyr::mutate(res, coverage = as.numeric(n_cells) / 4100)
  res
}

maude_results = pip_summary(p_10$maude)
mageck_results = pip_summary(p_10$mageck)
wb_results = pip_summary(p_10$waterbear)

wb_results = group_by(wb_results, effect_size_type, n_cells, seed)
wb_results = mutate(wb_results, pip_mapping = pip_group(pip))

wb_results = group_by(wb_results, effect_size_type, n_cells, seed, pip_mapping)

wb_results = summarize(wb_results, sensitivity = sum(true_effect) / length(true_effect),
  tp = sum(true_effect), fp = sum(!true_effect), fdr = fp / length(true_effect))
wb_results = dplyr::mutate(wb_results, coverage = as.numeric(n_cells) / 4100)


p = ggplot(wb_results, aes(pip_mapping, tp, fill = factor(coverage)))
p = p + geom_boxplot()
p = p + facet_wrap(~effect_size_type)
p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p = p + theme_light()
p = p + scale_y_continuous(breaks = seq(0, 100, length.out = 21))
p_wb_tp = p

p = ggplot(mageck_results, aes(pip_mapping, tp, fill = factor(coverage)))
p = p + geom_boxplot()
p = p + facet_wrap(~effect_size_type)
p = p + theme_light()
p = p + scale_y_continuous(breaks = seq(0, 100, length.out = 21))
p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p_mageck_tp = p

p = ggplot(maude_results, aes(pip_mapping, tp, fill = factor(coverage)))
p = p + geom_boxplot()
p = p + facet_wrap(~effect_size_type)
p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p = p + theme_light()
p = p + scale_y_continuous(breaks = seq(0, 100, length.out = 21))
p_maude_tp = p

save_plot(
  'img/pip_tp.pdf',
  plot_grid(p_wb_tp, p_mageck_tp, p_maude_tp,
    labels = c('waterbear', 'MAGeCK', 'MAUDE'),
    nrow = 3),
  base_height = 21
)

p = ggplot(wb_results, aes(pip_mapping, tp / (fp + tp), fill = factor(coverage)))
p = p + geom_boxplot()
p = p + facet_wrap(~effect_size_type)
p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p = p + theme_light()
p = p + scale_y_continuous(breaks = seq(0, 100, length.out = 21))
p_wb_tp = p

p = ggplot(mageck_results, aes(pip_mapping, tp / (fp + tp), fill = factor(coverage)))
p = p + geom_boxplot()
p = p + facet_wrap(~effect_size_type)
p = p + theme_light()
p = p + scale_y_continuous(breaks = seq(0, 100, length.out = 21))
p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p_mageck_tp = p

p = ggplot(maude_results, aes(pip_mapping, tp / (fp + tp), fill = factor(coverage)))
p = p + geom_boxplot()
p = p + facet_wrap(~effect_size_type)
p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p = p + theme_light()
p = p + scale_y_continuous(breaks = seq(0, 100, length.out = 21))
p_maude_tp = p

save_plot(
  'img/pip_sensitivity.pdf',
  plot_grid(p_wb_tp, p_mageck_tp, p_maude_tp,
    labels = c('waterbear', 'MAGeCK', 'MAUDE'),
    nrow = 3),
  base_height = 21
)

p = ggplot(
  dplyr::filter(wb_results, pip_mapping != 'none'),
  aes(pip_mapping, fp, fill = factor(coverage)))
p = p + geom_boxplot()
p = p + facet_wrap(~effect_size_type)
p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p = p + ylim(0, 200)
p_wb_fp = p

p = ggplot(
  dplyr::filter(mageck_results, pip_mapping != 'none'),
  aes(pip_mapping, fp, fill = factor(coverage)))
p = p + geom_boxplot()
p = p + facet_wrap(~effect_size_type)
p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p = p + ylim(0, 200)
p_mageck_fp = p

p = ggplot(
  dplyr::filter(maude_results, pip_mapping != 'none'),
  aes(pip_mapping, fp, fill = factor(coverage)))
p = p + geom_boxplot()
p = p + facet_wrap(~effect_size_type)
p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p = p + ylim(0, 200)
p_maude_fp = p

save_plot(
  'img/pip_fp.pdf',
  plot_grid(p_wb_fp, p_mageck_fp, p_maude_fp,
    labels = c('waterbear', 'MAGeCK', 'MAUDE'),
    nrow = 3),
  base_height = 21
)

wb_results_summary = group_by(wb_results, effect_size_type, n_cells)
wb_results_summary = summaries(wb_results)

tmp = dplyr::filter(metadata, p_gene_effects == '0.10', effect_size_type == 'difficult',
  n_cells == '1025000')

tmp_summary = readRDS(tmp$summary_filename[4])
tmp_summary$all.chains[grepl('dispersion', rownames(tmp_summary$all.chains)), ]

tmp_summary$all.chains[grepl('sigma', rownames(tmp_summary$all.chains)), ]

tmp_summary$all.chains[grepl('psi', rownames(tmp_summary$all.chains)), ]

tmp_summary
# dealing with setting psi to the correct value

sn = samples_n$summary$all.chains
mean(sn[grepl('gene_inclusion', rownames(sn)), 'Mean'] > 0.90)
summary(sn[grepl('gene_inclusion', rownames(sn)), 'Mean'] )
sn[grepl('dispersion', rownames(sn)), ]
# end: dealing with setting psi to the correct value


wb_res = read_waterbear_results(tmp[4, ], oracle_mapping)
head(wb_res, 30)

p = ggplot(wb_res, aes(fdr, sensitivity))
p = p + geom_path()
p = p + xlim(0, 1)
save_plot('img/debug_250x_difficult_2.pdf', p)

readRDS(tmp$)
1025000 / (1000 * 4 + 100)

debug_fnames = Sys.glob('results/nimble/vary_coverage_*.rds')
debug_metadata = parse_coverage_fname(basename(debug_fnames))
debug_metadata = mutate(debug_metadata, filename = debug_fnames)
debug_metadata = dplyr::filter(debug_metadata, grepl('vary_coverage_[0-9]+_[0-9]+_.*', filename))
debug_metadata = dplyr::filter(debug_metadata, !grepl('_summary.rds', filename))
debug_metadata = dplyr::filter(debug_metadata, grepl('30000_15000', filename))
debug_metadata = dplyr::mutate(debug_metadata,
  summary_filename = paste0(dirname(filename), '/', sub('.rds', '_summary.rds', basename(filename))))

sim = readRDS(paste0('results/truth/nimble_', basename(debug_metadata$filename)[1]))

oracle_mapping = generate_oracle_mapping(sim$const$guide_to_gene,
  p_gene_effects, n_total_genes)

debug_res = lapply(1:nrow(debug_metadata),
  function(i) {
    print(i)
    info = debug_metadata[i, ]
    tmp = read_waterbear_results(info, oracle_mapping)
    tmp = dplyr::mutate(tmp, p_gene_effects = info$p_gene_effects[1])
    tmp
  })
debug_res = bind_rows(debug_res)
debug_res = mutate(debug_res, method = 'waterbear')
debug_res = mutate(debug_res, coverage = as.numeric(n_cells) / (as.integer(n_total_genes) * 4 + 100))

s_debug_res = group_by(debug_res, effect_size_type, n_cells, rank, p_gene_effects)
s_debug_res = summarize(s_debug_res, mean_fdr = mean(fdr), sd_fdr = sd(fdr),
  mean_sensitivity = mean(sensitivity), sd_sensitivity = sd(sensitivity))
s_debug_res = mutate(s_debug_res, coverage = as.numeric(n_cells) / (as.integer(n_total_genes) * 4 + 100))
s_debug_res = mutate(s_debug_res, method = 'waterbear')

# tmp = dplyr::filter(res, bin_configuration == 1, effect_size_type == 'easy')
# p = ggplot(tmp, aes(fdr, tp, color = moi), alpha = 0.2)
p = ggplot(s_debug_res, aes(mean_fdr, mean_sensitivity, color = factor(coverage)), alpha = 0.2)
p = p + geom_point()
# p = p + geom_line(size = 1.5)
p = p + geom_path(size = 1.5)
p = p + xlim(0, 1)
p = p + xlab('true false discovery rate')
p = p + scale_color_manual(values = cbb)
p = p + facet_wrap(~effect_size_type + p_gene_effects)
p1 = p
save_plot(paste0('img/debug_summary_', base, '.pdf'), p1, base_height = 7)

# TODO: remove this once MAUDE is working on everything
metadata = dplyr::filter(metadata, file.exists(maude_filename))


p_gene_effects = 0.10
n_total_genes = 1000

metadata1000 = dplyr::filter(metadata, n_total_genes == 1000,
  bin_configuration == 1, p_gene_effects == '0.10')

sim = readRDS(paste0('results/truth/nimble_', basename(metadata1000$filename)[1]))

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

res = lapply(1:nrow(metadata_filter),
  function(i) {
    print(i)
    info = metadata_filter[i, ]
    read_waterbear_results(info, oracle_mapping)
  })
res = bind_rows(res)
res = mutate(res, method = 'waterbear')
res = mutate(res, coverage = as.numeric(n_cells) / (as.integer(n_total_genes) * 4 + 100))

s_res = group_by(res, effect_size_type, n_cells, rank)
s_res = summarize(s_res, mean_fdr = mean(fdr), sd_fdr = sd(fdr),
  mean_sensitivity = mean(sensitivity), sd_sensitivity = sd(sensitivity))
s_res = mutate(s_res, coverage = as.numeric(n_cells) / (as.integer(n_total_genes) * 4 + 100))
s_res = mutate(s_res, method = 'waterbear')

base = paste0('vary_coverage_', p_gene_effects, '_', n_total_genes, '_',
  bin_configuration, '_', n_replicates)

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
save_plot(paste0('img/summary_', base, '.pdf'), p1, base_height = 7)

# tmp = dplyr::filter(res, bin_configuration == 1, effect_size_type == 'easy')
# p = ggplot(tmp, aes(fdr, tp, color = moi), alpha = 0.2)
tmp = all_res
# tmp = dplyr::filter(all_res, grepl('(waterbear|maude)', method))
# p = ggplot(all_res, aes(mean_fdr, mean_sensitivity,
p = ggplot(tmp, aes(mean_fdr, mean_sensitivity,
    color = factor(coverage), linetype = method), alpha = 0.2)
# p = p + geom_point()
# p = p + geom_line(size = 1.5)
p = p + geom_path(size = 1.5)
# p = p + xlim(0, 0.15)
p = p + xlim(0, 1)
p = p + xlab('true false discovery rate')
p = p + scale_color_manual(values = cbb)
# p = p + facet_wrap(~effect_size_type)
p = p + facet_wrap(~effect_size_type + method)
p1 = p
save_plot(paste0('img/summary_', base, '.pdf'), p1, base_height = 7)


mageck_results = lapply(1:nrow(metadata_filter),
  function(i) {
    print(i)
    info = metadata_filter[i, ]
    process_mageck_results(info, oracle_mapping)
  })
mageck_results = bind_rows(mageck_results)
mageck_results = mutate(mageck_results, method = 'MAGeCK')
mageck_results = mutate(mageck_results, coverage = as.numeric(n_cells) / (as.integer(n_total_genes) * 4 + 100))

s_mageck_res = group_by(mageck_results, effect_size_type, n_cells, rank)
s_mageck_res = summarize(s_mageck_res, mean_fdr = mean(fdr), sd_fdr = sd(fdr),
  mean_sensitivity = mean(sensitivity), sd_sensitivity = sd(sensitivity))
s_mageck_res = mutate(s_mageck_res, coverage = as.numeric(n_cells) / (as.integer(n_total_genes) * 4 + 100))
s_mageck_res = mutate(s_mageck_res, method = 'MAGeCK')


debugonce(process_maude_results)
tmp = process_maude_results(metadata_filter[1, ])

maude_results = lapply(1:nrow(metadata_filter),
  function(i) {
    print(i)
    info = metadata_filter[i, ]
    process_maude_results(info, oracle_mapping)
  })
maude_results = bind_rows(maude_results)
maude_results = mutate(maude_results, method = 'maude')
maude_results = mutate(maude_results, coverage = as.numeric(n_cells) / (as.integer(n_total_genes) * 4 + 100))

s_maude_res = group_by(maude_results, effect_size_type, n_cells, rank)
s_maude_res = summarize(s_maude_res, mean_fdr = mean(fdr), sd_fdr = sd(fdr),
  mean_sensitivity = mean(sensitivity), sd_sensitivity = sd(sensitivity))
s_maude_res = mutate(s_maude_res, coverage = as.numeric(n_cells) / (as.integer(n_total_genes) * 4 + 100))
s_maude_res = mutate(s_maude_res, method = 'maude')

all_res = bind_rows(s_res, s_mageck_res, s_maude_res)


########################################################################
# looking at higher proportion
########################################################################

p_gene_effects = 0.25
pge = as.character(p_gene_effects)
n_total_genes = 1000


metadata1000 = dplyr::filter(metadata, n_total_genes == 1000,
  bin_configuration == 1, p_gene_effects == pge)

sim = readRDS(paste0('results/truth/nimble_', basename(metadata1000$filename)[1]))
# sim = readRDS(paste0('results/truth/nimble_', basename(metadata_large$filename)[1]))

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

# debugonce(read_waterbear_results)
res = lapply(1:nrow(metadata_filter),
  function(i) {
    print(i)
    info = metadata_filter[i, ]
    read_waterbear_results(info)
  })
res = bind_rows(res)
res = mutate(res, method = 'waterbear')
res = mutate(res, coverage = as.numeric(n_cells) / (as.integer(n_total_genes) * 4 + 100))

s_res = group_by(res, effect_size_type, n_cells, rank)
s_res = summarize(s_res, mean_fdr = mean(fdr), sd_fdr = sd(fdr),
  mean_sensitivity = mean(sensitivity), sd_sensitivity = sd(sensitivity),
  mean_estimated_fdr = mean(value), sd_estimated_fdr = sd(value))
s_res = mutate(s_res, coverage = as.numeric(n_cells) / (as.integer(n_total_genes) * 4 + 100))
s_res = mutate(s_res, method = 'waterbear')

all_estimated_fdr = lapply(c(0.01, 0.05, 0.10),
  function(l) {
    e_fdr = filter(
      dplyr::filter(
        group_by(s_res, effect_size_type, coverage, method),
        mean_estimated_fdr <= l),
      rank == max(rank)
    )
    mutate(e_fdr, level = l)
  })
all_estimated_fdr = bind_rows(all_estimated_fdr)

tmp = average_fdr_cutoff(res)

base = paste0('vary_coverage_', p_gene_effects, '_', n_total_genes, '_',
  bin_configuration, '_', n_replicates)

mageck_results = lapply(1:nrow(metadata_filter),
  function(i) {
    print(i)
    info = metadata_filter[i, ]
    process_mageck_results(info)
  })
mageck_results = bind_rows(mageck_results)
mageck_results = mutate(mageck_results, method = 'MAGeCK')
mageck_results = mutate(mageck_results, coverage = as.numeric(n_cells) / (as.integer(n_total_genes) * 4 + 100))

all_res = bind_rows(res, mageck_results)

s_res = group_by(all_res, method, effect_size_type, n_cells, rank)
s_res = summarize(s_res, mean_fdr = mean(fdr), sd_fdr = sd(fdr),
  mean_sensitivity = mean(sensitivity), sd_sensitivity = sd(sensitivity),
  mean_estimated_fdr = mean(value), sd_estimated_fdr = sd(value))
s_res = mutate(s_res, coverage = as.numeric(n_cells) / (as.integer(n_total_genes) * 4 + 100))

# tmp = dplyr::filter(res, bin_configuration == 1, effect_size_type == 'easy')
# p = ggplot(tmp, aes(fdr, tp, color = moi), alpha = 0.2)
all_estimated_fdr = dplyr::filter(all_estimated_fdr, coverage <= 5000)
p = ggplot(
  dplyr::filter(s_res, coverage <= 5000),
  aes(mean_fdr, mean_sensitivity, color = factor(coverage)), alpha = 0.2)
# p = p + geom_point()
# p = p + geom_line(size = 1.5)
p = p + geom_path(aes(linetype = method), size = 1.5)
p = p + geom_point(aes(shape = factor(level)), data = all_estimated_fdr, size = 5)
p = p + geom_errorbar(
  aes(
    ymin = mean_sensitivity - sd_sensitivity,
    ymax = mean_sensitivity + sd_sensitivity), data = all_estimated_fdr, size = 2, width = 0.3)
p = p + geom_errorbarh(
  aes(
    xmin = mean_fdr - sd_fdr,
    xmax = mean_fdr + sd_fdr), data = all_estimated_fdr, size = 2, width = 0.3)
# p = p + xlim(0, 0.25)
p = p + xlim(0, 0.15)
p = p + xlab('true false discovery rate')
p = p + scale_color_manual(values = cbb)
p = p + facet_wrap(~method + effect_size_type)
p1 = p
save_plot(paste0('img/summary_', base, '.pdf'), p1, base_height = 7)

all_res = bind_rows(s_res, s_mageck_res)

# investigating the issue with coverage == 50

metadata_debug = dplyr::filter(metadata1000, n_cells == 205000,
  effect_size_type == 'easy')

res_debug = dplyr::filter(res, n_cells == 205000, effect_size_type == 'easy')

get_estimated_fdr_cutoff(res_debug, 0.10)



# seems to be an issue with all of the samples
p = ggplot(res_debug, aes(fdr, sensitivity, color = factor(seed)))
p = p + geom_point()
p = p + geom_path(size = 1.5)
p = p + xlim(0, 0.15)
p = p + xlab('true false discovery rate')
# p = p + scale_color_manual(values = cbb)
p1 = p
save_plot(paste0('img/debug_summary_', base, '.pdf'), p1, base_height = 7)

res_debug = dplyr::filter(res_debug, seed == 4)
res_debug_mageck = dplyr::filter(mageck_results, seed == 4, n_cells == 205000,
  effect_size_type == 'easy')

get_estimated_fdr_cutoff(res_debug, 0.10)

tmp = nest_by(res, seed) %>%
  mutate(estimated_fdr = list(get_estimated_fdr_cutoff(data, 0.10)))

mutate(nest_by(res, seed), estimated_fdr = get_estimated_fdr_cutoff(data, 0.10))



summarize(res_debug, get_estimated_fdr_cutoff(data, 010))

mutate(get_estimated_fdr_cutoff(across(res_debug_mageck), 0.10))

get_estimated_fdr_cutoff(res_debug_mageck, 0.10)

metadata_debug = dplyr::filter(metadata1000, n_cells == 205000,
  effect_size_type == 'easy', seed == 4)


truth = readRDS(paste0('results/truth/sim_', basename(metadata_debug$filename)[1]))



fnames_large = Sys.glob('results/nimble/vary_coverage_*50000_25000*.rds')

metadata_large = parse_coverage_fname(basename(fnames_large))
metadata_large = mutate(metadata_large, filename = fnames_large)
metadata_large = dplyr::filter(metadata_large, grepl('vary_coverage_[0-9]+_[0-9]+_.*', filename))
metadata_large = dplyr::filter(metadata_large, !grepl('_summary.rds', filename))
metadata_large = dplyr::mutate(metadata_large,
  summary_filename = paste0(dirname(filename), '/', sub('.rds', '_summary.rds', basename(filename))))
metadata_large = dplyr::filter(metadata_large, file.exists(summary_filename))


res_large = lapply(1:nrow(metadata_large),
  function(i) {
    print(i)
    info = metadata_large[i, ]
    read_waterbear_results(info)
  })
res_large = bind_rows(res_large)
res_large = mutate(res_large, method = 'waterbear')
res_large = mutate(res_large, coverage = as.numeric(n_cells) / (as.integer(n_total_genes) * 4 + 100))

s_res = group_by(res_large, method, effect_size_type, n_cells, rank)
s_res = summarize(s_res, mean_fdr = mean(fdr), sd_fdr = sd(fdr),
  mean_sensitivity = mean(sensitivity), sd_sensitivity = sd(sensitivity),
  mean_estimated_fdr = mean(value), sd_estimated_fdr = sd(value))
s_res = mutate(s_res, coverage = as.numeric(n_cells) / (as.integer(n_total_genes) * 4 + 100))

# seems to be an issue with all of the samples
p = ggplot(res_large, aes(fdr, sensitivity, color = factor(seed)))
p = p + geom_point()
p = p + geom_path(size = 1.5)
p = p + xlim(0, 0.15)
p = p + xlab('true false discovery rate')
# p = p + scale_color_manual(values = cbb)
p1 = p
save_plot(paste0('img/debug_summary_', base, '.pdf'), p1, base_height = 7)

# seems to be an issue with all of the samples
p = ggplot(s_res, aes(mean_fdr, mean_sensitivity))
p = p + geom_point()
p = p + geom_path(size = 1.5)
p = p + xlim(0, 0.15)
p = p + xlab('true false discovery rate')
# p = p + scale_color_manual(values = cbb)
p1 = p
save_plot(paste0('img/5debug_summary_mean_', base, '.pdf'), p1, base_height = 7)

# TODO: look carefully at the true effect sizes and see how accurate the inference is

########################################################################

metadata1000 = dplyr::filter(metadata, n_mcmc_samples == 50000, n_total_genes == 1000,
  effect_size_type == 'easy', bin_configuration == 1)

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
n_replicates = '3'

tmp_filter = data.frame(
  p_gene_effects = p_gene_effects,
  n_total_genes = n_total_genes,
  n_mcmc_samples = n_mcmc_samples,
  n_mcmc_warmup = n_mcmc_warmup,
  bin_configuration = bin_configuration,
  n_replicates = n_replicates
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
mageck_results = mutate(mageck_results, coverage = as.numeric(n_cells) / (as.integer(n_total_genes) * 4 + 100))

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

tmp = bind_rows(mageck_results, res)


# tmp = dplyr::filter(res, bin_configuration == 1, effect_size_type == 'easy')
# p = ggplot(tmp, aes(fdr, tp, color = moi), alpha = 0.2)
p = ggplot(tmp, aes(fdr, sensitivity, color = factor(coverage), linetype = method), alpha = 0.2)
# p = ggplot(res, aes(fdr, sensitivity, color = factor(coverage)), alpha = 0.2)
p = p + geom_point()
# p = p + geom_line(size = 1.5)
p = p + geom_path(size = 1.5)
p = p + xlim(0, 0.15)
p = p + xlab('true false discovery rate')
p = p + scale_color_manual(values = cbb)
p1 = p
save_plot('img/coverage_easy.pdf', p1, base_height = 7)

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
save_plot('img/coverage_difficult.pdf', p1, base_height = 7)

########################################################################

metadata1000 = dplyr::filter(metadata, n_mcmc_samples == 50000, n_total_genes == 1000,
  effect_size_type == 'moderate', bin_configuration == 1)

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
n_replicates = '3'

tmp_filter = data.frame(
  p_gene_effects = p_gene_effects,
  n_total_genes = n_total_genes,
  n_mcmc_samples = n_mcmc_samples,
  n_mcmc_warmup = n_mcmc_warmup,
  bin_configuration = bin_configuration,
  n_replicates = n_replicates
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

tmp = bind_rows(res, mageck_results)
tmp = filter(tmp, as.integer(as.character(coverage)) <= 1000)
tmp = mutate(tmp, method = factor(method, levels = c('waterbear', 'MAGeCK'), labels = c('4 bins', '2 bins')))

# tmp = dplyr::filter(res, bin_configuration == 1, effect_size_type == 'easy')
# p = ggplot(tmp, aes(fdr, tp, color = moi), alpha = 0.2)
# p = ggplot(res, aes(fdr, sensitivity, color = factor(coverage)), alpha = 0.2)
p = ggplot(tmp, aes(fdr, sensitivity, color = factor(coverage), linetype = method), alpha = 0.2)
# p = p + geom_point()
# p = p + geom_line(size = 1.5)
p = p + geom_path(size = 1.5)
p = p + xlim(0, 0.15)
p = p + xlab('false discovery rate')
p = p + scale_color_manual(values = cbb, name = 'Cells per guide')
p = p + scale_linetype_manual(values = c(1, 3), name = 'Number of\nfractions (bins)')
p = p + theme_cowplot(20)
p1 = p
save_plot('img/coverage_moderate.pdf', p1, base_height = 7)
save_plot('img/coverage_moderate.png', p1, base_height = 7)

save.image('coverage_moderate.RData')

load('coverage_moderate.RData')

########################################################################
# prop 50
########################################################################

metadata1000 = dplyr::filter(metadata, n_mcmc_samples == 50000, n_total_genes == 1000,
  effect_size_type == 'easy', bin_configuration == 1)

p_gene_effects = 0.50
n_total_genes = 1000

sim = readRDS(paste0('results/coverage/truth/nimble_vary_coverage_2_4_1000_100_0.50_3_50000_25000_3_410000_easy_200.00_4.rds'))

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
save_plot(paste0('img/coverage_easy_', p_gene_effects, '.pdf'), p1, base_height = 7)

########################################################################

metadata1000 = dplyr::filter(metadata, n_mcmc_samples == 50000, n_total_genes == 1000,
  effect_size_type == 'difficult', bin_configuration == 1)

p_gene_effects = 0.50
n_total_genes = 1000

sim = readRDS(paste0('results/coverage/truth/nimble_vary_coverage_2_4_1000_100_0.50_3_50000_25000_3_410000_easy_200.00_4.rds'))

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
save_plot(paste0('img/coverage_difficult_', p_gene_effects, '.pdf'), p1, base_height = 7)


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

metadata1000 = dplyr::filter(metadata, n_total_genes == 1000,
  bin_configuration == 1, p_gene_effects == '0.10')

metadata1000 = dplyr::filter(metadata, n_mcmc_samples == 50000, n_total_genes == 1000,
  effect_size_type == 'easy', bin_configuration == 1)

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

tmp = dplyr::filter(res, coverage == '500')
tmp$true_effect %>% sum


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

# equalmixture summary
wb_results = dplyr::filter(p_10$waterbear, effect_size_type == 'equalmixture')
wb_results = dplyr::group_by(wb_results, n_cells, seed, type)

wb_sensitivity = summarize(
  dplyr::filter(wb_results, true_effect), sensitivity = mean(true_effect & value <= 0.10))
wb_sensitivity = dplyr::mutate(wb_sensitivity, method = 'waterbear')

mageck_results = dplyr::filter(p_10$mageck, effect_size_type == 'equalmixture')
mageck_results = dplyr::group_by(mageck_results, n_cells, seed, type)

mageck_sensitivity = summarize(
  dplyr::filter(mageck_results, true_effect), sensitivity = mean(true_effect & value <= 0.10))
mageck_sensitivity = dplyr::mutate(mageck_sensitivity, method = 'MAGeCK')

maude_results = dplyr::filter(p_10$maude, effect_size_type == 'equalmixture')
maude_results = dplyr::group_by(maude_results, n_cells, seed, type)

maude_sensitivity = summarize(
  dplyr::filter(maude_results, true_effect), sensitivity = mean(true_effect & value <= 0.10))
maude_sensitivity = dplyr::mutate(maude_sensitivity, method = 'MAUDE')

all_em = bind_rows(
  wb_sensitivity,
  mageck_sensitivity,
  maude_sensitivity
)

equal_mixture_sensitivity = function(df, label, alpha = 0.10, over_everything = FALSE) {
  results = dplyr::filter(df, effect_size_type == 'equalmixture')
  if (over_everything) {
    results = dplyr::group_by(results, n_cells, seed)
  } else {
    results = dplyr::group_by(results, n_cells, seed, type)
  }

  res = dplyr::summarize(
    results,
    sensitivity = sum(true_effect & value < alpha) / sum(true_effect)
    # fdr = sum(!true_effect & (value <= alpha)) / sum(value <= alpha)
    # precision = sum(true_effect & value < alpha) / sum(value < alpha)
  )
  dplyr::mutate(res, method = label)
}

equal_mixture_fdr = function(df, label, alpha = 0.10) {
  results = dplyr::filter(df, effect_size_type == 'equalmixture')
  results = dplyr::group_by(results, n_cells, seed, type)

  res = dplyr::summarize(
    group_by(ungroup(results), n_cells, seed),
    fdr = sum(!true_effect & (value <= alpha)) / sum(value <= alpha),
    precision = sum(true_effect & value < alpha) / sum(value < alpha)
    )
  dplyr::mutate(res, method = label)
}

all_em = bind_rows(
  equal_mixture_sensitivity(p_10$waterbear, 'waterbear'),
  equal_mixture_sensitivity(p_10$mageck, 'MAGeCK'),
  equal_mixture_sensitivity(p_10$maude, 'MAUDE')
  )
# all_em = dplyr::filter(all_em, type != 'none')

all_em = dplyr::mutate(all_em, coverage = as.integer(n_cells) / 4100)
all_em_50x = dplyr::filter(all_em, coverage == 50)

p = ggplot(dplyr::filter(all_em_50x, type != 'none'),
  aes(as.factor(type), sensitivity, fill = type))
p = p + geom_boxplot()
save_plot('img/equalmixture_50x_sensitivity.pdf', p, base_height = 7)

p = ggplot(dplyr::filter(all_em, type == 'difficult'),
  aes(method, sensitivity, fill = method))
p = p + geom_boxplot()
p = p + facet_wrap(~n_cells + type, labeller = 'label_both')
save_plot('difficult_sensitivity.pdf', p, base_height = 7)

# p = ggplot(dplyr::filter(all_em, type == 'difficult'),
p = ggplot(
  dplyr::filter(all_em, type != 'none'),
  aes(as.factor(as.integer(n_cells) / 4100), sensitivity, fill = method))
p = p + geom_boxplot()
p = p + facet_wrap(~type, labeller = 'label_both', nrow = 3)
save_plot('sensitivity_merge.pdf', p, base_height = 9, base_width = 14)

pooled_sensitivity = bind_rows(
  equal_mixture_sensitivity(p_10$waterbear, 'waterbear', over_everything = TRUE),
  equal_mixture_sensitivity(p_10$mageck, 'MAGeCK', over_everything = TRUE),
  equal_mixture_sensitivity(p_10$maude, 'MAUDE', over_everything = TRUE)
)

p = ggplot(
  pooled_sensitivity,
  aes(as.factor(as.integer(n_cells) / 4100), sensitivity, fill = method))
p = p + geom_boxplot()
save_plot('pooled_sensitivity.pdf', p, base_height = 9, base_width = 14)



all_fdr = bind_rows(
  equal_mixture_fdr(p_10$waterbear, 'waterbear'),
  equal_mixture_fdr(p_10$mageck, 'MAGeCK'),
  equal_mixture_fdr(p_10$maude, 'MAUDE')
  )

p = ggplot(all_fdr,
  aes(as.factor(as.integer(n_cells) / 4100), fdr, fill = method))
p = p + geom_boxplot()
p = p + ylim(0, 1)
# p = p + facet_wrap(~n_cells, labeller = 'label_both')
save_plot('difficult_fdr.pdf', p, base_height = 7)

p = ggplot(all_fdr,
  aes(method, precision, fill = method))
p = p + geom_boxplot()
p = p + facet_wrap(~n_cells, labeller = 'label_both')
save_plot('difficult_precision.pdf', p, base_height = 7)
