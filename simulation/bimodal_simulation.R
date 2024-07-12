#' the workhorse for generating bimodal distributions
#'
#' @param n_cells the total number of cells
#' @param n_guides the total number of guides
#' @param f_n_guides_per_cell a function that simply takes one parameter
bimodal_simulation = function(
  experiment_params,
  bimodal_description,
  seed = 1) {
  # TODO: check parameters in `experiment_params`
  # set.seed(seed)

  bins = format_bins(experiment_params$bins, experiment_params$n_replicates)

  true_parameters = data.frame(guide_id = 1:(experiment_params$n_guides))
  n_cells = rpois(experiment_params$n_replicates, experiment_params$n_cells)

  # generate the probability of being in the 'right'  distribution
  mixing_pi = bimodal_description$f_pi(experiment_params, next_seed())
  true_parameters$pi1 = c(mixing_pi$pi_effect, mixing_pi$pi_null)
  true_parameters = dplyr::mutate(true_parameters, pi0 = 1 - pi1)
  # transform into logistic space
  true_parameters = dplyr::mutate(true_parameters, b = safe_logit(pi1))

  initialize_reporter = bimodal_description$f_init_reporter(experiment_params, next_seed())
  true_parameters = data.frame(true_parameters, initialize_reporter$guide_effects)
  null_reporter_distribution = initialize_reporter$null_distribution
  # for now, all of the null reporters share the same pi distribution
  # print(mixing_pi$pi_null)
  null_reporter_distribution$pi1 = unique(mixing_pi$pi_null)
  null_reporter_distribution$pi0 = 1 - null_reporter_distribution$pi1


  p_guide = MCMCpack::rdirichlet(1,
    rep(experiment_params$guide_concentration, experiment_params$n_guides))
  p_guide = as.numeric(p_guide)
  p_guide = p_guide / sum(p_guide)

  # now have enough to simulate
  reporter_burn_in = matrix(0,
    nrow = experiment_params$burn_in,
    ncol = experiment_params$n_replicates)

  # compute the burn in
  # browser()
  for (r in 1:(experiment_params$n_replicates)) {
    # generate the number of guides per each cell
    n_guides_per_cell = bimodal_description$f_n_guides_per_cell(
      experiment_params$burn_in, next_seed())
    # for (i in 1:(n_cells[r])) {

    # browser()
    # reporter_burn_in[i, r] = reporter_sample
    tmp = generate_bimodal_cells(
      n_guides_per_cell,
      true_parameters,
      null_reporter_distribution[r, ],
      p_guide)
    reporter_burn_in[, r] = tmp$obs_reporter

    # for (i in 1:(experiment_params$burn_in)) {
    #   which_guides = rmultinom(1, n_guides_per_cell[i], p_guide)
    #   # for now, simply on or off
    #   which_guides = which_guides > 0
    #   mode_probability = logistic(sum(true_parameters$b[which_guides]))
    #   which_mode = sample(c(0, 1), 1, p = c(1 - mode_probability, mode_probability))
    #   reporter_sample = NULL
    #   if (which_mode) {
    #     # right
    #     total_effect = sum(true_parameters$effect_right[which_guides])
    #     reporter_sample = rnorm(1,
    #       mean = null_reporter_distribution$ss_right_mean[r] + total_effect,
    #       sd = null_reporter_distribution$ss_right_sd[r])
    #   } else {
    #     # left
    #     total_effect = sum(true_parameters$effect_left[which_guides])
    #     reporter_sample = rnorm(1,
    #       mean = null_reporter_distribution$ss_left_mean[r] + total_effect,
          # sd = null_reporter_distribution$ss_left_sd[r])
      # }
      # reporter_burn_in[i, r] = reporter_sample
    # }
  }

  # browser()
  bin_quantiles = lapply(bins,
    function(b) {
      as.numeric(as.matrix(b[, c('lower', 'upper')]))
    })
  # bin_quantiles = as.numeric(as.matrix(bins[, c('lower', 'upper')]))
  cutoffs = lapply(1:(experiment_params$n_replicates),
    function(r) {
      cur = quantile(reporter_burn_in[, r], probs = bin_quantiles[[r]])
      cur[bin_quantiles[[r]] == 0] = -Inf
      cur[bin_quantiles[[r]] == 1] = Inf
      matrix(cur, ncol = 2)
    })


  n_bins = nrow(bins[[1]])
  guide_counts = array(0,
    dim = c(experiment_params$n_replicates, experiment_params$n_guides, n_bins + 1))
  dimnames(guide_counts)[[3]] = c(paste0('bin_', 1:n_bins), 'unobserved')
  for (r in 1:(experiment_params$n_replicates)) {
    n_guides_per_cell = bimodal_description$f_n_guides_per_cell(n_cells[r], next_seed())

    # browser()
    tmp = generate_bimodal_cells(
      n_guides_per_cell,
      true_parameters,
      null_reporter_distribution[r, ],
      p_guide,
      cutoffs[[r]])
    guide_counts[r, , ] = tmp$bin_counts
    # for (i in 1:(n_cells[r])) {
    #   if (i %% 1e6 == 0) {
    #     message(paste0('replicate: ', r, ' cell number ', i))
    #   }
    #   which_guides = as.numeric(rmultinom(1, n_guides_per_cell[i], p_guide))
    #   # for now, simply on or off
    #   which_guides = which_guides > 0
    #   mode_probability = logistic(sum(true_parameters$b[which_guides]))
    #   which_mode = sample(c(0, 1), 1, p = c(1 - mode_probability, mode_probability))
    #   reporter_sample = NULL
    #   if (which_mode) {
    #     # right
    #     total_effect = sum(true_parameters$effect_right[which_guides])
    #     reporter_sample = rnorm(1,
    #       mean = null_reporter_distribution$ss_right_mean[r] + total_effect,
    #       sd = null_reporter_distribution$ss_right_sd[r])
    #   } else {
    #     # left
    #     total_effect = sum(true_parameters$effect_left[which_guides])
    #     reporter_sample = rnorm(1,
    #       mean = null_reporter_distribution$ss_left_mean[r] + total_effect,
    #       sd = null_reporter_distribution$ss_left_sd[r])
    #   }
    #   which_bin = cutoffs[[r]][, 1] <= reporter_sample & reporter_sample  < cutoffs[[r]][, 2]
    #   which_bin = which(which_bin)
    #   if (length(which_bin)) {
    #     guide_counts[r, which_guides, which_bin] = guide_counts[
    #       r, which_guides, which_bin] + 1
    #   } else {
    #     guide_counts[r, which_guides, n_bins + 1] = guide_counts[
    #       r, which_guides, n_bins + 1] + 1
    #   }
    # }
    colnames(guide_counts[r, , ]) = c(paste0('bin_', 1:n_bins), 'unobserved')
  }

  list(
    true_parameters = true_parameters,
    null_reporter_distribution = null_reporter_distribution,
    p_guide = p_guide,
    reporter_burn_in = reporter_burn_in,
    cutoffs = cutoffs,
    guide_counts = guide_counts,
    n_cells = n_cells
    )
}


next_seed = function() {
  sample.int(.Machine$integer.max, 1)
}

format_bins = function(bins, n_replicates) {
  if (length(bins) == 1) {
    bins = lapply(1:n_replicates,
      function(x) {
        bins[[1]]
      })
  }
  bins = lapply(seq_along(bins),
    function(r) {
      ret = lapply(bins[[r]],
        function(b) {
          stopifnot(length(b) == 2)
          data.frame(matrix(b, ncol = 2))
        })
      ret = dplyr::bind_rows(ret)
      colnames(ret) = c('lower', 'upper')
      dplyr::mutate(ret, bin_id = 1:nrow(ret))
    })
  # colnames(bins) = c('lower', 'upper')
  # bins = dplyr::mutate(bins, bin_id = 1:nrow(bins))

  bins
}

safe_logit = function(p) {
  x = qlogis(p)
  if (any(!is.finite(x))) {
    pos_on_bounds = x > 0 & !is.finite(x)
    neg_on_bounds = x < 0 & !is.finite(x)
    p[pos_on_bounds] = p[pos_on_bounds] - .Machine$double.eps
    p[neg_on_bounds] = p[neg_on_bounds] + .Machine$double.eps
    x = qlogis(p)
  }

  x
}


logistic = function(x) {
  plogis(x)
}

generate_1_guide_per_cell = function(n_cells) {
  rep.int(1L, n_cells)
}



#' simulate the number of guides per cell
#'
#' uses the zero-truncated Poisson distribution to generate the number of guides.
#' currently uses the `actuar` package.
#'
#' @param n_cells the total number of samples
#' @param lambda the mean under the _unconditional_ distribution
#' @param seed the seed
#' @return an integer vector of the number of guides per cell
generate_n_guides_per_cell = function(
  n_cells,
  lambda,
  seed = NULL
  ) {
  stopifnot(n_cells == round(n_cells))
  stopifnot(lambda > 0)

  if (is.null(seed)) {
    seed = 932
    warning(paste0('seed not set. defaulting to: ', seed))
  }
  set.seed(seed)

  n_cells = as.integer(n_cells)
  n_guides = actuar::rztpois(n_cells, lambda)

  n_guides
}

#' generate a mixture probability using the beta distribution
#'
#' @param n_total the total number of guides
#' @param n_effect the total number of guides that will have an effect
#' @param a_null shape1 of the beta distribution as parameterized by `rbeta`
#' @param b_null shape2 of the beta distribution as parameterized by `rbeta`
#' @param a_effect shape1 of the beta distribution as parameterized by `rbeta`
#' @param b_effect shape2 of the beta distribution as parameterized by `rbeta`
#' @param seed the seed to set. if not set, will default to something
#' @param config NOT IMPLEMENTED YET. TODO: if set, will override the remaining arguments
#' @return n_total - n_effect null probabilities and n_effect effect probabilities.
#' Also returns the configuration
generate_pi_beta = function(
  n_null,
  n_effect,
  a_null,
  b_null,
  a_effect,
  b_effect,
  seed = NULL,
  config = NULL
  ) {
  # TODO: if config is not NULL, fill in all of the values
  stopifnot(round(n_null) == n_null)
  stopifnot(round(n_effect) == n_effect)
  stopifnot(a_null >= 0 && b_null >= 0)
  stopifnot(a_effect >= 0 && b_effect >= 0)

  n_null = as.integer(n_null)
  n_effect = as.integer(n_effect)

  if (is.null(seed)) {
    seed = 326
    warning(paste0('generate_pi_beta: seed not set. using: ', seed))
  }
  set.seed(seed)

  n_total = n_null + n_effect
  p_effects = n_effect / n_total

  pi_null = rbeta(n_null, a_null, b_null)
  pi_effect = rbeta(n_effect, a_effect, b_effect)

  config = list(
    n_total = n_total,
    n_effect = n_effect,
    a_null = a_null,
    b_null = b_null,
    a_effect = a_effect,
    b_effect = b_effect,
    seed = seed)

  list(
    pi_null = pi_null,
    pi_effect = pi_effect,
    config = config)
}

# the output probability is the probability of being 'on the right'
generate_pi_unimodal = function(
  n_total,
  seed = NULL,
  config = NULL) {
  rep.int(0, n_total)
}

cell_simulation_to_nimble = function(parameters, sim) {
  guide_names = 1:parameters$n_guides
  n_targeting = parameters$n_guides - parameters$n_control
  guide_to_gene = c(
    rep(1:(n_targeting / parameters$n_guides_per_target),
      each = parameters$n_guides_per_target),
    rep('control', parameters$n_control))
  gene_mapping = data.frame(guide = as.character(guide_names),
    i = guide_names, gene = guide_to_gene)

  nt_guide_names = dplyr::filter(gene_mapping, grepl('control', gene))
  test_guide_names = dplyr::filter(gene_mapping, !grepl('control', gene))
  test_guide_names = dplyr::mutate(test_guide_names,
    mapping = as.integer(factor(as.integer(gene))))
  test_guide_names = dplyr::arrange(test_guide_names, i)

  data = list(x = sim$counts)
  const = list(
    N = dim(sim$counts)[1],
    N_guides = nrow(test_guide_names),
    x_total = apply(sim$counts, c(1, 2), sum),
    N_genes = length(unique(test_guide_names$mapping)),
    N_nt = nrow(nt_guide_names),
    N_bins = dim(sim$counts)[3],
    nt_index = 1:nrow(nt_guide_names),
    nt_data_index = nt_guide_names$i,
    guide_index = 1:nrow(test_guide_names),
    guide_data_index = test_guide_names$i,
    guide_to_gene = test_guide_names$mapping,
    dispersion_prior_mean = 210
    )
  const$N_cutoffs = const$N_bins - 1
  dispersion_init = 210

  if (length(parameters$bin_sizes) == 1) {
    bin_alpha = matrix(rep(unlist(parameters$bin_sizes), const$N), nrow = const$N, byrow = TRUE)
    parameters$bin_sizes = lapply(1:const$N, function(x) parameters$bin_sizes[[1]])
  } else {
    bin_alpha = matrix(unlist(parameters$bin_sizes), nrow = const$N, byrow = TRUE)
  }


  init = list(
    bin_alpha = bin_alpha,
    psi = 0.2,
    sigma_gene = 5,
    sigma_guide = 5,
    dispersion = dispersion_init
    )

  init$cutoffs = matrix(
    sapply(1:nrow(init$bin_alpha),
      function(i) {
        dirichlet_to_normal_bins(init$bin_alpha[i, ])
      }),
      nrow = const$N, byrow = TRUE)
  # init$cutoffs = matrix(rep(dirichlet_to_normal_bins(init$bin_alpha[1, ]), const$N),
  # nrow = const$N, byrow = TRUE)
  q_init = array(0, dim = c(const$N, const$N_guides, const$N_bins))
  for (n in 1:dim(q_init)[1]) {
    for (g in 1:dim(q_init)[2]) {
      q_init[n, g, ] = parameters$bin_sizes[[n]]
    }
  }
  init$q = q_init
  init$gene_inclusion = rep(0, const$N_genes)
  init$gene_shift = rep(0, const$N_genes)
  init$guide_shift = rep(0, const$N_guides)
  init$total_shift = rep(0, const$N_guides)
  const$bin_alpha_prior = init$bin_alpha

  list(
    init = init,
    const = const,
    data = data
    )
}

