library('dplyr')
library('tidyr')
library('ggplot2')
library('RColorBrewer')
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

# process_mageck_results = function(metadata_row) {
#   info = metadata_row
#   fname = nimble_fname_to_mageck(info$filename)
#   fname = paste0('results/mageck/', fname)
#   mageck = read_mageck(fname)
#   mageck = inner_join(mageck, oracle_mapping, by = 'gene_id')
#   benchmark = compute_benchmark_statistics(mageck)
#   dplyr::mutate(benchmark, bin_configuration = info$bin_configuration[1],
#     seed = info$seed,
#     moi = info$moi, effect_size_type = info$effect_size_type,
#     n_cells = info$n_cells)
# }

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
  (as.integer(n_total_genes) * as.integer(n_guides_per_target) + as.(n_control_guides)))

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

color_map = c('MAGeCK' = '#7570B3', 'waterbear' = '#D95F02', 'MAUDE' = '#c0c0c0')


# these are the current plots
x_dodge = 0.0
p = ggplot(tmp_summary,
  aes(method, m_sensitivity, color = method))
p = p + geom_point()
p = p + geom_pointrange(aes(ymin = m_sensitivity - m_95_sensitivity,
    ymax = m_sensitivity + m_95_sensitivity),
  )
p = p + scale_x_discrete(expand=c(1, 1))
# p = p + scale_color_brewer(palette = 'Dark2')
p = p + scale_color_manual(values = color_map)
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
save_plot('img/grid_unequal_average_sensitivity_compare_all_methods.pdf', p)

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
save_plot('img/grid_unequal_average_sensitivity_compare_no_maude.pdf', p)

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
# p = p + scale_color_brewer(palette = 'Dark2')
p = p + scale_color_manual(values = color_map)
p = p + xlab('coverage')
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
# p = p + scale_color_brewer(palette = 'Dark2')
p = p + scale_color_manual(values = color_map)
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

p = ggplot(
  dplyr::filter(tmp_summary,
    method != 'MAUDE'),
  aes(method, m_fdr, color = method))
p = p + geom_point()
p = p + geom_pointrange(aes(ymin = m_fdr - m_95_fdr,
    ymax = m_fdr + m_95_fdr))
p = p + scale_x_discrete(expand=c(1, 1))
# p = p + scale_color_brewer(palette = 'Dark2')
p = p + scale_color_manual(values = color_map)
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
save_plot('img/grid_unequal_average_fdr_compare_no_maude.pdf', p)


