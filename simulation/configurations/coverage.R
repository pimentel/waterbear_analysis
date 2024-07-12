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

get_filename = function() {
  paste(
    'coverage',
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
parse_coverage = function(fname) {
  original_fname = fname
  fname = sub('coverage_', '', fname)
  fname = sub('.rds', '', fname)
  split_s = strsplit(fname, '_')
  ret = lapply(split_s,
    function(s) {
      params = data.frame(
        seed = as.integer(s[1]),
        n_chains = as.integer(s[2]),
        n_iter = as.integer(s[3]),
        n_warmup = as.integer(s[4]),
        n_samples = as.integer(s[5]),
        p_effects = as.numeric(s[6]),
        n_genes = as.integer(s[7]),
        n_guides_per_gene = as.integer(s[8]),
        bin_configuration = as.integer(s[9]),
        n_cells_per_guide = as.integer(s[10]),
        stringsAsFactors = FALSE)
      params
    })
  ret = bind_rows(ret)
  dplyr::mutate(ret, filename = original_fname)
}
