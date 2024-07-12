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

process_mageck_results = function(metadata_row, oracle_mapping) {
  info = metadata_row
  fname = nimble_fname_to_mageck(info$filename)
  # fname = paste0('results/mageck/', fname)
  mageck = read_mageck(fname)
  mageck = inner_join(mageck, oracle_mapping, by = 'gene_id')
  benchmark = compute_benchmark_statistics(mageck)
  dplyr::mutate(benchmark, bin_configuration = info$bin_configuration[1],
    seed = info$seed,
    moi = info$moi, effect_size_type = info$effect_size_type,
    n_cells = info$n_cells)
}


process_maude_results = function(metadata_row, oracle_mapping) {
  info = metadata_row
  maude = read.csv(info$maude_filename, header = TRUE, stringsAsFactors = FALSE)
  maude = dplyr::mutate(maude, Gene = as.character(Gene))
  maude = dplyr::select(maude, gene_id = Gene, value = fdr)
  maude = inner_join(maude, oracle_mapping, by = 'gene_id')
  benchmark = compute_benchmark_statistics(maude)
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
  file.path(
    sub('nimble', 'mageck', dirname(fname)),
    sub('.rds$', '.gene_summary.txt', basename(fname))
  )
}

read_waterbear_results = function(metadata_row, oracle_mapping, recode = FALSE) {
  info = metadata_row
  if (!is.null(info$posterior_summary)) {
    gene_inclusion = readRDS(info$posterior_summary)
    gene_inclusion = inner_join(
      dplyr::mutate(gene_inclusion, gene_id = as.character(gene_id)),
      oracle_mapping)
    gene_inclusion = mutate(gene_inclusion, value = lfsr)
    gene_inclusion = arrange(gene_inclusion, value, desc(abs(gene_shift)))
  } else {
    x = readRDS(info$summary_filename)

    # wo = readRDS(paste0('results/truth/wo_', basename(info$filename)[1]))
    # gene_inclusion = waterbear:::wb_recode(x, wo)
    # gene_inclusion = dplyr::mutate(gene_inclusion, value = 1 - Mean)
    # gene_inclusion = dplyr::rename(gene_inclusion, gene_id = gene)

    # TODO: ideally, add the local false sign rate to WB output
    all_chains = x[['all.chains']]
    var_names = rownames(all_chains)
    gene_inclusion = all_chains[grepl('gene_inclusion', var_names), ]
    gene_inclusion = data.frame(variable = rownames(gene_inclusion), gene_inclusion)
    gene_inclusion = dplyr::mutate(gene_inclusion, gene_id = strip_text(variable),
      value = 1 - Mean)
    gene_inclusion = inner_join(gene_inclusion, oracle_mapping, by = c('gene_id'))

    gene_effect = all_chains[grepl('gene_shift', var_names), ]
    gene_effect = data.frame(variable = rownames(gene_effect), gene_effect)
    gene_effect = dplyr::mutate(gene_effect, gene_id = strip_text(variable),
      waterbear_effect_size = Mean, waterbear_95_lower = `X95.CI_low`,
      waterbear_95_upper = `X95.CI_upp`)
    gene_effect = dplyr::mutate(gene_effect,
      ci_95_significant = (waterbear_95_lower < 0 & waterbear_95_upper < 0) |
        (0 < waterbear_95_lower & 0 < waterbear_95_upper))
    gene_inclusion = left_join(gene_inclusion,
      dplyr::select(gene_effect, gene_id, waterbear_effect_size, ci_95_significant),
      by = c('gene_id'))
    if (any(gene_inclusion$value < 0.10 & !gene_inclusion$ci_95_significant)) {
      warning('CI not significant but value is!')
    }
    gene_inclusion = arrange(gene_inclusion, value, desc(abs(waterbear_effect_size)))
  }

  # gene_inclusion = all_chains[grepl('gene_inclusion', var_names), ]
  # gene_inclusion = data.frame(variable = rownames(gene_inclusion), gene_inclusion)
  # gene_inclusion = dplyr::mutate(gene_inclusion,
  #   gene_id = as.character(as.factor(strip_text(variable))),
  #   value = 1 - Mean)
  # gene_inclusion = inner_join(gene_inclusion, oracle_mapping, by = c('gene_id'))
  # gene_inclusion = arrange(gene_inclusion, value)
  benchmark = compute_benchmark_statistics(gene_inclusion)
  res = dplyr::mutate(benchmark, bin_configuration = info$bin_configuration[1],
    seed = info$seed,
    moi = info$moi, effect_size_type = info$effect_size_type,
    n_cells = info$n_cells, n_guides_per_target = info$n_guides_per_target)
  if ('coverage_per_gene' %in% colnames(info)) {
    res = dplyr::mutate(res, coverage_per_gene = info$coverage_per_gene)
  }
  res
}

get_estimated_fdr_cutoff = function(sensitivity_results, level) {
  res = dplyr::mutate(sensitivity_results, estimated_level = as.character(level))
  res = dplyr::filter(res, value <= level)
  if (nrow(res) > 0) {
    res = dplyr::filter(res, rank == max(rank))
  } else {
    res = sensitivity_results[1, ]
    res = dplyr::mutate(res, estimated_level = as.character(level))
    res = dplyr::mutate(res, value = NA, fdr = NA, sensitivity = NA, rank = NA, Mean = NA)
  }
  res
}

get_estimated_fdr_cutoff = function(sensitivity_results, level) {
  res = dplyr::mutate(sensitivity_results, estimated_level = as.character(level))
  res = dplyr::filter(res, value <= level)
  if (nrow(res) > 0) {
    res = dplyr::filter(res, rank == max(rank))
  } else {
    res = sensitivity_results[1, ]
    res = dplyr::mutate(res, estimated_level = as.character(level))
    res = dplyr::mutate(res, value = NA, fdr = NA, sensitivity = NA, rank = NA, Mean = NA)
  }
  res
}

average_fdr_cutoff = function(res) {
  r = lapply(c(0.01, 0.05, 0.10),
    function(alpha) {
      est_fdr = mutate(
        nest_by(res, effect_size_type, coverage, seed),
        estimated_fdr = (get_estimated_fdr_cutoff(data, alpha)))
      summarize(group_by(est_fdr, effect_size_type, coverage, estimated_fdr$estimated_level),
        mean_fdr = mean(estimated_fdr$fdr, na.rm = TRUE),
        sd_fdr = sd(estimated_fdr$fdr, na.rm = TRUE),
        mean_sensitivity = mean(estimated_fdr$sensitivity, na.rm = TRUE))
    })
  bind_rows(r)
}

compute_gene_effect_sizes = function(true_parameters, n_guides_per_gene, n_control_guides) {
  n_guides = nrow(true_parameters)
  n_test_guides = n_guides - n_control_guides
  n_genes = n_test_guides / n_guides_per_gene
  true_parameters = dplyr::mutate(true_parameters,
    gene = c(rep(1:n_genes, each = n_guides_per_gene), rep(0, n_control_guides)))
  true_parameters = dplyr::group_by(true_parameters, gene)
  true_parameters = dplyr::mutate(true_parameters, gene_effect_size = mean(effect_left))
  true_parameters
}


classify_effect_sizes = function(effect_sizes, true_effect_size_fname) {
  load(true_effect_size_fname)
  epsilon = 0.01
  equal_mixture_summary = group_by(equal_mixture_summary, type)
  equal_mixture_range = summarize(equal_mixture_summary,
    lower = min(abs(mean)), upper = max(abs(mean)))
  equal_mixture_range = bind_rows(
    data.frame(type = c('none'), lower = 0, upper = min(equal_mixture_range$lower)),
    equal_mixture_range)
  equal_mixture_range$upper[
    which(max(equal_mixture_range$upper) == equal_mixture_range$upper)] = 999
  classification = vector('character', length(effect_sizes))
  # print(equal_mixture_range)
  for (i in seq_along(effect_sizes)) {
    tmp = equal_mixture_range$lower - epsilon <= abs(effect_sizes[i]) &
      abs(effect_sizes[i]) < equal_mixture_range$upper + epsilon
    if (!any(tmp)) {
      print(equal_mixture_range)
      stop(paste0("error! can't classify this effect size: ", effect_sizes[i]))
    }
    classification[i] = equal_mixture_range$type[which(tmp)[1]]
    # classification[i] = equal_mixture_range$type[tmp]

  }
  # TODO: make figures that break it into different groups
  data.frame(gene_effect_sizes = effect_sizes, type = classification)
}

get_effect_size_information = function(info, oracle_mapping) {
  truth_fname = sub('/nimble/', '/truth/sim_', info$filename)
  truth = readRDS(truth_fname)
  tmp = compute_gene_effect_sizes(truth$true_parameters,
    as.integer(info$n_guides_per_target), as.integer(info$n_control_guides))
  effect_size_types = classify_effect_sizes(tmp$gene_effect_size,
    'coverage_effect_sizes.RData')
  effect_size_types = dplyr::mutate(effect_size_types, gene_id = as.character(tmp$gene))
  effect_size_types = dplyr::distinct(dplyr::select(effect_size_types, gene_id, type))

  dplyr::inner_join(oracle_mapping, effect_size_types, by = 'gene_id')
}

# get_effect_size_information(metadata[1, , drop = FALSE], data.frame())
