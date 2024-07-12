install.packages(c(
    'MCMCpack',
    'dplyr',
    'nimble',
    'extraDistr'
    ))
# assumes it has been run from `run_dm.R`
# assumes the following variables have been set:
# - seed
# - n_chains
# - n_iter
# - n_warmup

n_cells_per_guide = as.integer(args[1])
n_genes = as.integer(args[2])
p_effects = as.numeric(args[3])
n_samples = as.integer(args[4])
bin_configuration = as.integer(args[5])
which_effects = as.character(args[6])

if (!(which_effects %in% c('easy', 'moderate', 'difficult'))) {
  stop(paste0('invalid which_effects: ', which_effects))
}

get_bin_configuration = function(bin_mapping) {
  bin_mapping = as.integer(bin_mapping)
  if (bin_mapping == 1) {
    return(c(0.15, 0.35, 0.35, 0.15))
  } else if (bin_mapping == 2) {
    return(c(0.2, 0.3, 0.3, 0.2))
  } else if (bin_mapping == 3) {
    return(c(0.25, 0.25, 0.25, 0.25))
  }
  stop('unrecognized bin_mapping mode.')
}

bin_sizes = get_bin_configuration(bin_configuration)

# computed by a_posterior * number of guides (5998)
guide_rate_total_mass = 38782.1

parameters = list(
  sigma_guide = 0.08732188,
  sigma_gene = 0.2381384,
  dispersion = 200,
  n_guides_per_gene = 4,
  n_reporter_cells = 1e6,
  guide_rate_dispersion = guide_rate_total_mass,
  n_control_guides = 100
  )

parameters$n_samples = n_samples
parameters$n_genes = n_genes
parameters$p_effects = p_effects
parameters$n_guides_per_gene = 4
parameters$bin_sizes = bin_sizes
parameters$n_cells_per_guide = n_cells_per_guide

load('./coverage_effect_sizes.RData')

if (which_effects == 'easy') {
  parameters$gene_effect_distribution = all_easy_summary$mean
} else if (which_effects == 'moderate') {
  parameters$gene_effect_distribution = all_moderate_summary$mean
} else if (which_effects == 'difficult') {
  parameters$gene_effect_distribution = all_difficult_summary$mean
} else {
  stop(paste0('unrecognized which_effects: ', which_effects))
}

get_filename = function() {
  paste(
    'coverage_observed_effect',
    which_effects,
    seed,
    n_chains,
    n_iter,
    n_warmup,
    n_samples,
    sprintf('%.2f', round(p_effects, 2)),
    n_genes,
    parameters$n_guides_per_gene,
    bin_configuration,
    n_cells_per_guide
    , sep = '_', collapse = '')
}

# TODO: write a parser for the file names
parse_coverage_observed_effect = function(fname) {
  original_fname = fname
  fname = sub('coverage_observed_effect_', '', fname)
  fname = sub('.rds', '', fname)
  split_s = strsplit(fname, '_')
  ret = lapply(split_s,
    function(s) {
      params = data.frame(
        which_effects = as.character(s[1]),
        seed = as.integer(s[2]),
        n_chains = as.integer(s[3]),
        n_iter = as.integer(s[4]),
        n_warmup = as.integer(s[5]),
        n_samples = as.integer(s[6]),
        p_effects = as.numeric(s[7]),
        n_genes = as.integer(s[8]),
        n_guides_per_gene = as.integer(s[9]),
        bin_configuration = as.integer(s[10]),
        n_cells_per_guide = as.integer(s[11]),
        stringsAsFactors = FALSE)
      params
    })
  ret = bind_rows(ret)
  dplyr::mutate(ret, filename = original_fname)
}

coverage_observed_effect_filename_to_params = function(fname) {
  df = parse_coverage_observed_effect(fname)
  parameters = list()
  parameters$n_cells_per_guide = df$cells_per_guide[1]
  parameters$n_genes = df$n_genes
  parameters$p_effects = df$p_effects
  parameters$n_samples = df$n_samples
  parameters$bin_configuration =
}
