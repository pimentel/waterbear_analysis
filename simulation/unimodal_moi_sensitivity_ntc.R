library('dplyr')
library('tidyr')
library('ggplot2')
library('cowplot')

cbb <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

source('benchmark_common.R')

theme_set(theme_cowplot(25))

strip_text = function(x) {
  y = sub('.+\\[', '', x)
  sub('\\]', '', y)
}

compute_benchmark_statistics = function(fdr_table) {
  # assumed you have grouped on appropriate columns (at least the method)
  stopifnot(!is.null(fdr_table))

  fdr_table = dplyr::arrange(fdr_table, value)
  fdr_table = dplyr::mutate(fdr_table,
    fp = cumsum(!true_effect), tp = cumsum(true_effect),
    fn = sum(true_effect) - tp, rank = 1:length(true_effect),
    sensitivity = tp / (fn + tp), fdr = fp / rank, est_fdr = value)

  fdr_table
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
      s_res = group_by(res, n_control_guides, effect_size_type, n_moi, n_cells, n_guides_per_target, rank)
      s_res = summarize(s_res, mean_fdr = mean(fdr), sd_fdr = sd(fdr),
        mean_sensitivity = mean(sensitivity), sd_sensitivity = sd(sensitivity),
        mean_estimated_fdr = mean(value), sd_estimated_fdr = sd(value))
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
      s_mageck_res = group_by(mageck_results, n_control_guides, effect_size_type, n_moi, n_cells, n_guides_per_target, rank)
      s_mageck_res = summarize(s_mageck_res, mean_fdr = mean(fdr), sd_fdr = sd(fdr),
        mean_sensitivity = mean(sensitivity), sd_sensitivity = sd(sensitivity),
        mean_estimated_fdr = mean(value), sd_estimated_fdr = sd(value))
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
      s_maude_res = group_by(maude_results, n_control_guides, effect_size_type, n_moi, n_cells, n_guides_per_target, rank)
      s_maude_res = summarize(s_maude_res, mean_fdr = mean(fdr), sd_fdr = sd(fdr),
        mean_sensitivity = mean(sensitivity), sd_sensitivity = sd(sensitivity),
        mean_estimated_fdr = mean(value), sd_estimated_fdr = sd(value))
      s_maude_res = mutate(s_maude_res, method = 'MAUDE')
    } else {
      s_maude_res = mutate(maude_results, method = 'MAUDE')
    }
    all_res$maude = s_maude_res
  }
  all_res
}

fnames = Sys.glob('results/nimble/unimodal_moi_*.rds')

parse_unimodal_moi_fname = function(x) {
  x = sub('unimodal_moi_', '', x)
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
        moi = s[10],
        n_cells = s[11],
        effect_size_type = s[12],
        dispersion = s[13],
        n_mcmc_chains = s[14],
        stringsAsFactors = FALSE)
    })
  dplyr::bind_rows(ret)
}

metadata = parse_unimodal_moi_fname(basename(fnames))
metadata = mutate(metadata, filename = fnames)

# first, get only the ones with varying control guides
metadata = filter(metadata, p_gene_effects == '0.10', n_total_genes == '1000',
  n_replicates == '3', n_mcmc_samples == '20000',
  bin_configuration == '1'
)
metadata = filter(metadata, grepl('_summary.rds', filename))
metadata = mutate(metadata, summary_filename = filename)
metadata = mutate(metadata, filename = sub('_summary', '', filename))

metadata = dplyr::filter(metadata, file.exists(summary_filename))
metadata$n_control_guides %>% table

tmp = filter(metadata, n_control_guides == '10')

# metadata = parse_unimodal_moi_fname(basename(fnames))
# metadata = mutate(metadata, filename = fnames)

# metadata = dplyr::filter(metadata, !grepl('_summary.rds', filename))
# metadata = dplyr::filter(metadata, effect_size_type == 'equalmixture')

p_gene_effects = '0.10'
n_total_genes = '1000'
bin_configuration = '1'
n_replicates = '3'

tmp = dplyr::filter(metadata, n_control_guides == '10')

tmp = dplyr::filter(metadata, n_control_guides == '100',
  moi == 0.3 | moi == 1 | moi == 2 | moi == 5 | moi == 10)


tmp = dplyr::filter(metadata, n_control_guides == '1000')

# debugonce(process_subgroup_by_method)
# debugonce(get_effect_size_information)
# debugonce(classify_effect_sizes)
debugonce(read_waterbear_results)

ntc_10_waterbear = process_subgroup_by_method(
  dplyr::filter(metadata,
    n_control_guides == '10'), p_gene_effects, n_total_genes, bin_configuration, n_replicates,
  method = 'waterbear',
  summary = FALSE,
  debug = TRUE)
ntc_10_waterbear$waterbear$method = 'waterbear (10)'

ntc_100_waterbear = process_subgroup_by_method(
  dplyr::filter(metadata,
    n_control_guides == '100',
    moi == 0.3 | moi == 1 | moi == 2 | moi == 5 | moi == 10
    # moi == 0.3
    ), p_gene_effects, n_total_genes, bin_configuration, n_replicates,
  method = 'waterbear',
  summary = FALSE,
  debug = TRUE)
ntc_100_waterbear$waterbear$method = 'waterbear (100)'

ntc_1000_waterbear = process_subgroup_by_method(
  dplyr::filter(metadata,
    n_control_guides == '1000',
    moi == 0.3 | moi == 1 | moi == 2 | moi == 5 | moi == 10
    ), p_gene_effects, n_total_genes, bin_configuration, n_replicates,
  method = 'waterbear',
  summary = FALSE,
  debug = TRUE)
ntc_1000_waterbear$waterbear$method = 'waterbear (1000)'

p_10_mageck = process_subgroup_by_method(metadata, p_gene_effects, n_total_genes, bin_configuration, n_replicates,
  method = 'mageck',
  summary = FALSE,
  debug = TRUE)

all_res = bind_rows(p_10_maude$maude, p_10_waterbear$waterbear, p_10_mageck$mageck)

all_res = bind_rows(ntc_10_waterbear$waterbear, ntc_100_waterbear$waterbear, ntc_1000_waterbear$waterbear)

alpha = 0.10
tmp = group_by(all_res,
  method, n_control_guides, seed, effect_size_type, moi, n_cells, n_guides_per_target)
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

z_95 = abs(qnorm(0.025))
tmp_summary = group_by(tmp, method, moi, n_cells) %>%
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

tmp_summary = ungroup(tmp_summary)
tmp_summary = mutate(tmp_summary,
  f_moi = factor(as.numeric(moi), levels = sort(as.numeric(unique(moi))), ordered = TRUE), n_moi = as.numeric(f_moi))
tmp_summary$f_moi

tmp_summary = mutate(tmp_summary,
  `number of cells` = factor(as.numeric(n_cells), levels = sort(as.numeric(unique(n_cells))), ordered = TRUE))

color_mapping = distinct(
  select(ungroup(tmp_summary), f_moi))
color_mapping = mutate(color_mapping, moi_color = rep(c('gray', 'white'), length.out = nrow(color_mapping)))
tmp_summary = left_join(tmp_summary, color_mapping, by = 'f_moi')

# these are the current plots
x_dodge = 0.5
p = ggplot(tmp_summary,
  aes(f_moi, m_sensitivity, color = method))
p = p + geom_point(position = position_dodge2(width = x_dodge, padding = x_dodge), alpha = 0)
p = p + geom_rect(
  aes(xmin = n_moi - x_dodge, xmax = n_moi + x_dodge, ymin = 0, ymax = 1, fill = f_moi), color = NA, alpha = 0.05, show.legend = FALSE)
p = p + scale_fill_manual('moi', values = color_mapping$moi_color)
p = p + geom_point(position = position_dodge2(width = x_dodge, padding = x_dodge))
p = p + geom_pointrange(aes(ymin = m_sensitivity - m_95_sensitivity,
    ymax = m_sensitivity + m_95_sensitivity),
  position = position_dodge2(width = x_dodge, padding = x_dodge))
p = p + scale_color_brewer(palette = 'Dark2')
p = p + xlab('moi')
p = p + ylab('sensitivity')
p = p + ylim(0, 1)
p = p + facet_wrap(~`number of cells`, nrow = 3, labeller = label_both)
save_plot('img/ntc_10_moi_average_sensitivity_compare.pdf', p,
  base_height = 14)

p = ggplot(tmp_summary,
  aes(f_moi, m_fdr, color = method))
p = p + geom_point(position = position_dodge2(width = x_dodge, padding = x_dodge), alpha = 0)
p = p + geom_rect(
  aes(xmin = n_moi - x_dodge, xmax = n_moi + x_dodge, ymin = 0, ymax = 1, fill = f_moi), color = NA, alpha = 0.05, show.legend = FALSE)
p = p + scale_fill_manual(values = color_mapping$moi_color)
p = p + geom_point(position = position_dodge2(width = 0.5, padding = 0.5))
p = p + geom_pointrange(aes(ymin = m_fdr - m_95_fdr,
    ymax = m_fdr + m_95_fdr),
  position = position_dodge2(width = 0.5, padding = 0.5))
p = p + scale_color_brewer(palette = 'Dark2')
p = p + xlab('moi')
p = p + ylab('false discovery rate')
p = p + ylim(0, 1)
p = p + geom_hline(yintercept = alpha, linetype = 3)
p = p + facet_wrap(~`number of cells`, nrow = 3, labeller = label_both)
save_plot('img/ntc_moi_average_fdr_compare.pdf', p,
  base_height = 14)


x_dodge = 0.5
p = ggplot(filter(tmp_summary, n_cells == 50000),
  aes(f_moi, m_sensitivity, color = method))
p = p + geom_point(position = position_dodge2(width = x_dodge, padding = x_dodge), alpha = 0)
p = p + geom_rect(
  aes(xmin = n_moi - x_dodge, xmax = n_moi + x_dodge, ymin = 0, ymax = 1, fill = f_moi), color = NA, alpha = 0.05)
p = p + scale_fill_manual(values = color_mapping$moi_color)
p = p + geom_point(position = position_dodge2(width = x_dodge, padding = x_dodge))
p = p + geom_pointrange(aes(ymin = m_sensitivity - m_95_sensitivity,
    ymax = m_sensitivity + m_95_sensitivity),
  position = position_dodge2(width = x_dodge, padding = x_dodge))
p = p + scale_color_brewer(palette = 'Dark2')
p = p + xlab('moi')
p = p + ylab('sensitivity')
p = p + ylim(0, 1)
save_plot('img/moi_50000_average_sensitivity_compare2.pdf', p)

# new plots
color_map = c('MAGeCK' = '#7570B3', 'waterbear' = '#D95F02', 'MAUDE' = '#c0c0c0')
p = ggplot(filter(tmp_summary, n_cells == 50000),
  aes(f_moi, m_fdr, color = method))
p = p + geom_point()
p = p + geom_pointrange(aes(ymin = m_fdr - m_95_fdr,
    ymax = m_fdr + m_95_fdr))
p = p + scale_x_discrete(expand=c(1, 1))
# p = p + scale_color_brewer(palette = 'Dark2')
p = p + scale_color_manual(values = color_map)
p = p + xlab('moi')
p = p + ylab('fdr')
p = p + ylim(0, 1)
p = p + geom_hline(yintercept = alpha, linetype = 3)
p = p + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(),
  panel.spacing = unit(0.1, "lines"),
  panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.1)
)
p = p + facet_wrap(~f_moi, nrow = 1)
save_plot('img/grid_moi_50000_average_fdr_compare.pdf', p)

# these are the current plots
x_dodge = 0.0
p = ggplot(
  dplyr::filter(tmp_summary, method != 'MAUDE', n_cells == 50000),
  aes(method, m_sensitivity, color = method))
# p = p + geom_point()
p = p + geom_pointrange(aes(ymin = m_sensitivity - m_95_sensitivity,
    ymax = m_sensitivity + m_95_sensitivity),
  )
p = p + scale_x_discrete(expand=c(1, 1))
p = p + scale_color_manual(values = color_map)
p = p + xlab('coverage')
p = p + ylab('sensitivity')
p = p + ylim(0, 1)
p = p + facet_wrap(~f_moi, nrow = 1)
p = p + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(),
  panel.spacing = unit(0.1, "lines"),
  panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.1)
)
save_plot('img/grid_moi_50000_average_sensitivity_compare2.pdf', p)

# these are the current plots
x_dodge = 0.0
p = ggplot(
  dplyr::filter(tmp_summary, n_cells == 50000),
  aes(method, m_sensitivity, color = method))
# p = p + geom_point()
p = p + geom_pointrange(aes(ymin = m_sensitivity - m_95_sensitivity,
    ymax = m_sensitivity + m_95_sensitivity),
  )
p = p + scale_x_discrete(expand=c(1, 1))
p = p + scale_color_manual(values = color_map)
p = p + xlab('coverage')
p = p + ylab('sensitivity')
p = p + ylim(0, 1)
p = p + facet_wrap(~f_moi, nrow = 1)
p = p + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(),
  panel.spacing = unit(0.1, "lines"),
  panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.1)
)
save_plot('img/grid_moi_50000_average_sensitivity_compare_with_maude.pdf', p)
