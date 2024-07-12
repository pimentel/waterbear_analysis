library('dplyr')
library('tidyr')
library('ggplot2')
library('cowplot')

cbb <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

theme_set(theme_cowplot())

source('benchmark_common.R')


########################################################################
# looking at coverage
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

metadata = parse_coverage_fname(basename(fnames))
metadata = mutate(metadata, filename = fnames)
metadata = dplyr::filter(metadata, grepl('vary_coverage_[0-9]+_[0-9]+_.*', filename))
metadata = dplyr::filter(metadata, !grepl('_summary.rds', filename))
metadata = dplyr::mutate(metadata,
  summary_filename = paste0(dirname(filename), '/', sub('.rds', '_summary.rds', basename(filename))))
metadata = dplyr::mutate(metadata,
  posterior_summary = sub('summary.rds', 'posterior_summary.rds', summary_filename))
metadata = dplyr::filter(metadata, file.exists(summary_filename))
metadata = dplyr::filter(metadata, file.exists(posterior_summary))

metadata = dplyr::mutate(metadata,
  maude_filename = paste0(sub('nimble', 'maude', dirname(filename)), '/',
    sub('.rds', '.csv', basename(filename))))

metadata = dplyr::filter(metadata, n_mcmc_samples == 20000)

# equal mixture
metadata_em = dplyr::filter(metadata, effect_size_type == 'equalmixture')

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

  print('waterbear')
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

  print('mageck')
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

        # TODO: classify effect sizes
        # browser()
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

    # s_mageck_res = group_by(mageck_results, effect_size_type, n_cells, rank)
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
  # all_res = bind_rows(s_res, s_mageck_res, s_maude_res)
  # if (debug) {
  #   return(list(all_res = all_res, waterbear = res, maude = maude_results,
  #       mageck = mageck_results))
  # }
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

# TODO: check on local false sign rate and whether or not need to compute
p_10 = process_subgroup(metadata_em, p_gene_effects, n_total_genes, bin_configuration, n_replicates,
  debug = TRUE)

debugonce(process_subgroup_by_method)

metadata_em_4_guides = dplyr::filter(metadata_em, n_guides_per_target == 4)

p_10_maude = process_subgroup_by_method(metadata_em_4_guides, p_gene_effects, n_total_genes, bin_configuration, n_replicates,
  method = 'maude',
  summary = FALSE,
  debug = TRUE)

p_10_waterbear = process_subgroup_by_method(metadata_em_4_guides, p_gene_effects, n_total_genes, bin_configuration, n_replicates,
  method = 'waterbear',
  summary = FALSE,
  debug = TRUE)

p_10_mageck = process_subgroup_by_method(metadata_em_4_guides, p_gene_effects, n_total_genes, bin_configuration, n_replicates,
  method = 'mageck',
  summary = FALSE,
  debug = TRUE)

# wb_res = p_10_waterbear$waterbear

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


summarize(group_by(tmp, as.factor(coverage), n_guides_per_target), length(coverage)) %>%
  filter(n_guides_per_target == 4)

tmp = filter(tmp, n_guides_per_target == '4', coverage <= 1000)

# this is where I left off

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


p = ggplot(tmp, aes(as.factor(coverage), sensitivity, fill = method))
p = p + geom_boxplot()
save_plot('img/tmp.pdf', p)


saveRDS(tmp_summary, file ='tmp_summary.rds')

tmp_summary = readRDS('tmp_summary.rds')

tmp_summary = mutate(tmp_summary, f_coverage = factor(coverage), n_coverage = as.numeric(f_coverage))
color_mapping = distinct(
  select(ungroup(tmp_summary), f_coverage))
color_mapping = mutate(color_mapping, coverage_color = rep(c('gray', 'white'), length.out = nrow(color_mapping)))
tmp_summary = left_join(tmp_summary, color_mapping, by = 'f_coverage')

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
save_plot('img/average_sensitivity_compare.pdf', p)

color_map = c('MAGeCK' = '#7570B3', 'waterbear' = '#D95F02', 'MAUDE' = '#c0c0c0')

p = ggplot(tmp_summary,
  aes(method, m_fdr, color = method))
p = p + geom_point()
p = p + geom_pointrange(aes(ymin = m_fdr - m_95_fdr,
    ymax = m_fdr + m_95_fdr))
p = p + scale_x_discrete(expand=c(1, 1))
# p = p + scale_color_brewer(palette = 'Dark2')
p = p + scale_color_manual(values = color_map)
p = p + xlab('Coverage')
p = p + ylab('True false discovery rate')
p = p + ylim(0, 1)
p = p + geom_hline(yintercept = alpha, linetype = 3)
p = p + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank())
# p = p + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(),
#   panel.spacing = unit(0.1, "lines"),
#   panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_rect(colour = "black", size=0.1)
# )
p = p + facet_wrap(~f_coverage, nrow = 1)
save_plot('img/grid_average_fdr_compare.pdf', p)

# these are the current plots
x_dodge = 0.0
p = ggplot(
  dplyr::filter(tmp_summary, method != 'MAUDE'),
  aes(method, m_sensitivity, color = method))
p = p + geom_point()
p = p + geom_pointrange(aes(ymin = m_sensitivity - m_95_sensitivity,
    ymax = m_sensitivity + m_95_sensitivity),
  )
p = p + scale_x_discrete(expand=c(1, 1))
p = p + scale_color_manual(values = color_map)
p = p + xlab('Coverage')
p = p + ylab('Sensitivity')
p = p + ylim(0, 1)
p = p + facet_wrap(~f_coverage, nrow = 1)
p = p + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank())
# p = p + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(),
#   panel.spacing = unit(0.1, "lines"),
#   panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_rect(colour = "black", size=0.1)
# )
save_plot('img/grid_average_sensitivity_compare.pdf', p)

# these are the current plots
x_dodge = 0.0
p = ggplot(
  tmp_summary,
  aes(method, m_sensitivity, color = method))
p = p + geom_point()
p = p + geom_pointrange(aes(ymin = m_sensitivity - m_95_sensitivity,
    ymax = m_sensitivity + m_95_sensitivity),
  )
p = p + scale_x_discrete(expand=c(1, 1))
p = p + scale_color_manual(values = color_map)
p = p + xlab('coverage')
p = p + ylab('sensitivity')
p = p + ylim(0, 1)
p = p + facet_wrap(~f_coverage, nrow = 1)
p = p + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank())
# p = p + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(),
#   panel.spacing = unit(0.1, "lines"),
#   panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_rect(colour = "black", size=0.1)
# )
save_plot('img/grid_average_sensitivity_compare_with_maude.pdf', p)
