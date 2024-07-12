fname_to_simulation_type = function(x) {
  sapply(x,
    function(f) {
      if (grepl('vary_coverage_[0-9]+_[0-9]+_', f)) {
        return('vary_coverage')
      } else if (grepl('vary_coverage_2_bins_', f)) {
        return('vary_coverage_2_bins')
      } else if (grepl('unimodal_moi_', f)){
        return('unimodal_moi')
      } else if (grepl('unequal_bins_', f)) {
        return('unequal_bins')
      } else {
        return('unidentified')
      }
    })
}

parse_filename = function(fname) {
  sim_type = fname_to_simulation_type(fname)
  switch(sim_type,
    vary_coverage = parse_coverage_fname,
    unimodal_moi = parse_unimodal_moi_fname,
    unequal_bins = parse_unequal_bins_fname
    )
}

fname_to_bins = function(fnames) {
  simulation_type = fname_to_simulation_type(fnames)
  lapply(seq_along(simulation_type),
    function(i) {
      sim_type = simulation_type[i]
      f = basename(fnames[i])
      file_metadata = parse_filename(f)(f)
      bins = get_bin_configuration(sim_type, file_metadata$bin_configuration)
      # browser()
      bin_sizes = lapply(bins, function(x) sapply(x, diff))
      if (all(sapply(bin_sizes[2:length(bin_sizes)], identical, bin_sizes[1]))) {
        bin_sizes = bin_sizes[[1]]
      }
      # bin_sizes = sapply(bins[[1]], diff)
      bin_sizes
    })
}

unequal_bins_env = new.env()
unequal_bins_env$bins_1 = list(
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

unequal_bins_env$bins_2 = list(
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

unequal_bins_env$bins_3 = list(
  list(
    c(0, 0.10), c(0.10, 0.40),
    c(0.40, 0.90), c(0.90, 1)),
  list(
    c(0, 0.15), c(0.15, 0.5),
    c(0.5, 0.9), c(0.9, 1)),
  list(
    c(0, 0.20), c(0.20, 0.45),
    c(0.45, 0.9), c(0.90, 1))
)

unequal_bins_env$bins_4 = list(
  list(
    c(0, 0.10), c(0.10, 0.50),
    c(0.50, 0.75), c(0.75, 1)),
  list(
    c(0, 0.25), c(0.25, 0.5),
    c(0.5, 0.75), c(0.75, 1)),
  list(
    c(0, 0.25), c(0.25, 0.5),
    c(0.50, 0.9), c(0.90, 1))
)

get_bin_configuration = function(simulation_type, bin_type) {
  bin_type = as.integer(bin_type)
  switch(simulation_type,
    'vary_coverage' = switch(bin_type,
      list(list(
          c(0, 0.15), c(0.15, 0.5),
          c(0.5, 0.85), c(1 - 0.15, 1))),
      list(list(
          c(0, 0.2), c(0.2, 0.5),
          c(0.5, 0.8), c(0.8, 1))),
      list(list(
          c(0, 0.25), c(0.25, 0.5),
          c(0.5, 0.75), c(0.75, 1)))
      ),
    'unimodal_moi' = switch(bin_type,
      list(list(
          c(0, 0.15), c(0.15, 0.5),
          c(0.5, 0.85), c(1 - 0.15, 1))),
      list(list(
          c(0, 0.2), c(0.2, 0.5),
          c(0.5, 0.8), c(0.8, 1))),
      list(list(
          c(0, 0.25), c(0.25, 0.5),
          c(0.5, 0.75), c(0.75, 1)))
    ),
    'unequal_bins' = switch(bin_type,
      {unequal_bins_env$bins_1},
      {unequal_bins_env$bins_2},
      {unequal_bins_env$bins_3},
      {unequal_bins_env$bins_4}
    )
  )
}

parse_coverage_fname = function(x) {
  x = sub('vary_coverage_', '', x)
  x = sub('(.rds|.tsv)', '', x)
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


parse_unimodal_moi_fname = function(x) {
  x = sub('unimodal_moi_', '', x)
  x = sub('(.rds|.tsv)', '', x)
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


parse_unequal_bins_fname = function(x) {
  x = sub('unequal_bins_', '', x)
  x = sub('(.rds|.tsv)', '', x)
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
