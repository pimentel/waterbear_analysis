library('dplyr')
library('tidyr')
library('ggplot2')
library('cowplot')

cbb <- c(
  "#F0E442",
  "#E69F00",
  "#009E73",
  "#56B4E9",
  "#0072B2", "#D55E00", "#CC79A7", "#000000")

theme_set(theme_cowplot(25))
# theme_set(theme_cowplot())

strip_text = function(x) {
  y = sub('^.+\\[', '', x)
  sub('\\]', '', y)
}

source('benchmark_common.R')



########################################################################
# 2 bin mode
########################################################################

fnames = Sys.glob('results/nimble/vary_coverage_2_bins*.rds')

fnames = fnames[grepl('equalmixture', fnames)]


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

metadata = dplyr::filter(metadata, !grepl('_summary.rds', filename))
metadata = dplyr::mutate(metadata,
  summary_filename = paste0(dirname(filename), '/', sub('.rds', '_summary.rds', basename(filename))))
metadata = dplyr::filter(metadata, file.exists(summary_filename))

########################################################################
# p 0.10
########################################################################

metadata1000 = dplyr::filter(metadata, n_total_genes == 1000,
  bin_configuration == 1, p_gene_effects == '0.10')

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


res = lapply(1:nrow(metadata_filter),
  function(i) {
    print(i)
    info = metadata_filter[i, ]
    read_waterbear_results(info, oracle_mapping)
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

# tmp = dplyr::filter(res, bin_configuration == 1, effect_size_type == 'easy')
# p = ggplot(tmp, aes(fdr, tp, color = moi), alpha = 0.2)
p = ggplot(s_res, aes(mean_fdr, mean_sensitivity, color = factor(coverage)), alpha = 0.2)
p = p + geom_point()
# p = p + geom_line(size = 1.5)
p = p + geom_path(size = 1.5)
p = p + xlim(0, 0.15p)
p = p + xlab('true false discovery rate')
p = p + scale_color_manual(values = cbb)
p = p + facet_wrap(~effect_size_type)
p1 = p
save_plot(paste0('img/summary_', base, '.pdf'), p1, base_height = 7)

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

# dealing with equal mixture
#
metadata_em = dplyr::filter(metadata, effect_size_type == 'equalmixture')

process_subgroup = function(metadata, p_gene_effects, n_total_genes, bin_configuration,
  n_replicates, debug = FALSE, waterbear_only = FALSE) {

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
  n_control_guides = metadata_filter$n_control_guides[1]
  n_guides_per_target = metadata_filter$n_guides_per_target[1]

  s_res = group_by(res, effect_size_type, n_cells, rank)
  s_res = summarize(s_res, mean_fdr = mean(fdr), sd_fdr = sd(fdr),
    mean_sensitivity = mean(sensitivity), sd_sensitivity = sd(sensitivity),
    mean_estimated_fdr = mean(value), sd_estimated_fdr = sd(value))
  s_res = mutate(s_res,
    coverage = as.numeric(n_cells) / (
      as.integer(n_total_genes) * as.integer(n_guides_per_target) + as.integer(n_control_guides)))
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
  } else {
    all_res = bind_rows(s_res)
  }


  all_res
}

p_gene_effects = '0.10'
n_total_genes = '1000'
bin_configuration = '1'
n_replicates = '3'

p_10_3_replicates = process_subgroup(metadata, p_gene_effects, n_total_genes,
  bin_configuration, n_replicates, waterbear_only = TRUE)

p = ggplot(
  dplyr::filter(p_10_3_replicates, method == 'waterbear', coverage < 50000),
  aes(mean_fdr, mean_sensitivity, color = factor(coverage), linetype = n_replicates), alpha = 0.2)
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
# save_plot(paste0('img/waterbear_equal_mixture_summary_compare_coverage_', base, '.pdf'), p1, base_height = 7)
save_plot(paste0('img/tmp_', base, '.pdf'), p1, base_height = 7)

fnames = Sys.glob('results/nimble/vary_coverage_*.rds')

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

metadata_4b = parse_coverage_fname(basename(fnames))
metadata_4b = mutate(metadata_4b, filename = fnames)
metadata_4b = dplyr::filter(metadata_4b, grepl('vary_coverage_[0-9]+_[0-9]+_.*', filename))
metadata_4b = dplyr::filter(metadata_4b, !grepl('_summary.rds', filename))
metadata_4b = dplyr::mutate(metadata_4b,
  summary_filename = paste0(dirname(filename), '/', sub('.rds', '_summary.rds', basename(filename))))
metadata_4b = dplyr::filter(metadata_4b, file.exists(summary_filename))

metadata_4b = dplyr::mutate(metadata_4b,
  maude_filename = paste0(sub('nimble', 'maude', dirname(filename)), '/',
    sub('.rds', '.csv', basename(filename))))

metadata_4b = dplyr::filter(metadata_4b, n_mcmc_samples == 20000)
metadata_4b = dplyr::filter(metadata_4b, effect_size_type == 'equalmixture')

p_10_3_replicates_4b = process_subgroup(metadata_4b, p_gene_effects, n_total_genes,
  bin_configuration, n_replicates, waterbear_only = TRUE)


tmp_10 = bind_rows(
  dplyr::mutate(p_10_3_replicates, n_bins = 2),
  dplyr::mutate(p_10_3_replicates_4b, n_bins = 4)
)

base = paste0('vary_coverage_2_vs_4_bins_', p_gene_effects, '_', n_total_genes, '_',
  bin_configuration, '_', n_replicates)

p = ggplot(
  dplyr::filter(tmp_10, method == 'waterbear', coverage < 50000),
  aes(mean_fdr, mean_sensitivity, color = factor(coverage), linetype = as.factor(n_bins)), alpha = 0.2)
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
# p = p + xlim(0, 0.15)
p = p + xlab('true false discovery rate')
p = p + scale_color_manual(values = cbb)
p1 = p
save_plot(paste0('img/waterbear_equal_mixture_summary_compare_2_vs_4_bins_', base, '.pdf'), p1, base_height = 7)



# NEW
summary_2_bins = readRDS('2_bins_coverage.rds')
summary_4_bins = readRDS('p_10_coverage.rds')
summary_4_bins_mode4 = readRDS('p_10_coverage_bin4.rds')
summary_2_bins_mode4 = readRDS('2_bins_coverage_bin4.rds')
summary_2_bins_mode5 = readRDS('2_bins_coverage_bin5.rds')

summary_4_bins = filter(summary_4_bins, method == 'waterbear')
summary_4_bins = mutate(summary_4_bins, coverage = as.integer(n_cells) / 4100)
summary_4_bins_mode4 = mutate(summary_4_bins_mode4, coverage = as.integer(n_cells) / 4100)

summary_2_bins = mutate(summary_2_bins, method = '2 bins (0.15)')
summary_4_bins = mutate(summary_4_bins, method = '4 bins (0.15)')
summary_2_bins_mode4 = mutate(summary_2_bins_mode4, method = '2 bins (0.05)')
summary_2_bins_mode5 = mutate(summary_2_bins_mode5, method = '2 bins (0.01)')
summary_4_bins_mode4 = mutate(summary_4_bins_mode4, method = '4 bins (0.05)')

all_summary = bind_rows(
  summary_2_bins,
  summary_2_bins_mode4,
  summary_2_bins_mode5,
  summary_4_bins,
  summary_4_bins_mode4)

tmp = filter(all_summary, coverage <=500, method != 'MAGeCK')

set.seed(20)


# tmp = dplyr::filter(res, bin_configuration == 1, effect_size_type == 'easy')
# p = ggplot(tmp, aes(fdr, tp, color = moi), alpha = 0.2)
p = ggplot(
  filter(all_summary, coverage <=500,
    method != 'MAGeCK'),
  aes(mean_fdr, mean_sensitivity, color = method), alpha = 0.2)
# p = p + geom_point()
# p = p + geom_line(size = 1.5)
p = p + geom_path(size = 1.5)
p = p + xlim(0, 0.15)
# p = p + xlim(0, 1)
p = p + xlab('True false discovery rate')
p = p + scale_color_manual(values = (cbb))
p = p + ylab('')
p = p + facet_wrap(~coverage)
p = p + theme(
  # panel.spacing = unit(2, "lines"),
  legend.position = c(0.2, 0.8),
  axis.text.x = element_text(angle = 30, vjust = 0.5))
p1 = p
save_plot('f5_img_bin_comparison.pdf', p1, base_height = 7)
