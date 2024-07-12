
# things to tune
########################################################################
# seed
# number of guides per target
# number of total genes
# number of control guides
# proportion of genes with effects
# number of replicates
# number of MCMC samples
# number of MCMC warm up samples
# bin configuration
# multiplicity of infection
# number of cells
# type of effect sizes (easy, moderate, difficult)
# dispersion (two significant digit)

# Error in job vary_coverage while creating output files /oak/stanford/groups/pritch/users/jake_harold/simulation/results/nimble/vary_coverage_9_1_20000_100_0.10_3_20000_10000_1_500000_equalmixture_200.00_4.rds, /oak/stanford/groups/pritch/users/jake_harold/simulation/results/truth/sim_vary_coverage_9_1_20000_100_0.10_3_20000_10000_1_500000_equalmixture_200.00_4.rds, /oak/stanford/groups/pritch/users/jake_harold/simulation/results/truth/nimble_vary_coverage_9_1_20000_100_0.10_3_20000_10000_1_500000_equalmixture_200.00_4.rds.
# ClusterJobException in line 500 of /oak/stanford/groups/pritch/users/jake_harold/simulation/Snakefile:
# Error executing rule vary_coverage on cluster (jobid: 10563, external: Submitted batch job 48263806, jobscript: /oak/stanford/groups/pritch/users/jake_harold/simulation/.snakemake/tmp.f1e77gks/snakejob.vary_coverage.10563.sh). For detailed error see the cluster log.


args = c(
  '9',
  '1',
  '20000',
  '100',
  '0.10',
  '3',
  '20000',
  '10000',
  '1',
  '500000',
  'equalmixture',
  '200.00',
  '4')

args = c('42',
  '4',
  '250',
  '100',
  '0.10',
  '2',
  '50000',
  '25000',
  '3',
  '410000',
  'easy',
  '200.00',
  '4')

# this one
args = c('42',
  '4',
  '1000',
  '100',
  '0.10',
  '3',
  # '50000',
  # '25000',
  '25000',
  '12500',
  '3',
  '2050000',
  'equalmixture',
  '200.00',
  '4')

# debugging NA
# [1, 2373, ]
args = c('3',
  '4',
  '1000',
  '100',
  '0.10',
  '3',
  '20000',
  '10000',
  '1',
  '410000',
  'difficult',
  '200.00',
  '4'
)

args = c('2',
  '4',
  '1000',
  '100',
  '0.50',
  '3',
  '10000',
  '5000',
  '1',
  '7333333',
  'easy',
  '200.00',
  '1')


n_mcmc_samples = 1000
n_mcmc_warmup = 500
n_mcmc_chains = 2

args = commandArgs(trailingOnly = TRUE)

seed = as.integer(args[1])
n_guides_per_target = as.integer(args[2])
n_total_genes = as.integer(args[3])
n_control_guides = as.integer(args[4])
p_gene_effects = as.numeric(args[5])
n_replicates = as.integer(args[6])
n_mcmc_samples = as.integer(args[7])
n_mcmc_warmup = as.integer(args[8])
bin_configuration = as.integer(args[9])
n_cells = as.integer(args[10])
effect_size_type = as.character(args[11])
dispersion = as.numeric(args[12])
n_mcmc_chains = as.integer(args[13])

fbase = paste(
  'vary_coverage',
  seed,
  n_guides_per_target,
  n_total_genes,
  n_control_guides,
  sprintf('%.2f', round(p_gene_effects, 2)),
  n_replicates,
  n_mcmc_samples,
  n_mcmc_warmup,
  bin_configuration,
  n_cells,
  effect_size_type,
  sprintf('%.2f', round(dispersion, 2)),
  n_mcmc_chains,
  sep = '_')

source('bimodal_simulation.R')

# computed by a_posterior * number of guides (5998)
guide_rate_total_mass = 38782.1

bimodal_parameters = list()
bimodal_parameters$n_guides = (n_guides_per_target * n_total_genes) + n_control_guides
bimodal_parameters$n_control = n_control_guides
bimodal_parameters$n_guides_per_target = n_guides_per_target
bimodal_parameters$n_cells = n_cells
bimodal_parameters$n_genes = n_total_genes
bimodal_parameters$p_gene_effects = p_gene_effects
bimodal_parameters$n_gene_effects = bimodal_parameters$n_genes * bimodal_parameters$p_gene_effects
bimodal_parameters$n_replicates = n_replicates
bimodal_parameters$burn_in = 1e6

bins = switch(bin_configuration,
  list(list(
    c(0, 0.15), c(0.15, 0.5),
    c(0.5, 0.85), c(1 - 0.15, 1))),
  list(list(
    c(0, 0.2), c(0.2, 0.5),
    c(0.5, 0.8), c(0.8, 1))),
  list(list(
    c(0, 0.25), c(0.25, 0.5),
    c(0.5, 0.75), c(0.75, 1))),
  list(list(
    c(0, 0.05), c(0.05, 0.5),
    c(0.5, 0.95), c(0.95, 1))),
  list(list(
    c(0, 0.01), c(0.01, 0.5),
    c(0.5, 0.99), c(0.99, 1)))
)

bimodal_parameters$bins = bins
bimodal_parameters$bin_sizes = lapply(bimodal_parameters$bins,
  function(b) {
    sapply(
      b,
      function(x) x[2] - x[1])
  })
bimodal_parameters$null_reporter_left_mean = 0
bimodal_parameters$null_reporter_left_sd = 1
bimodal_parameters$null_reporter_right_mean = NA
bimodal_parameters$null_reporter_right_sd = NA
bimodal_parameters$guide_concentration = guide_rate_total_mass / bimodal_parameters$n_guides
bimodal_parameters$sigma_guide = 0.08732188

bimodal_init = list()
bimodal_init$f_n_guides_per_cell = function(n_cells, seed) {
  generate_1_guide_per_cell(n_cells)
}
bimodal_init$f_pi = function(bimodal_parameters, seed) {
  # the pi parameter is shared among all the null guides
  n_total = bimodal_parameters$n_guides
  res = list()
  n_effect = with(bimodal_parameters, n_gene_effects * n_guides_per_target)
  n_null = bimodal_parameters$n_guides - n_effect
  res$pi_null = generate_pi_unimodal(n_null)
  res$pi_effect = generate_pi_unimodal(n_effect)
  res
}

load('../simulation/coverage_effect_sizes.RData')
effect_distribution = switch(effect_size_type,
  easy = all_easy_summary$mean,
  moderate = all_moderate_summary$mean,
  difficult = all_difficult_summary$mean,
  equalmixture = equal_mixture_summary$mean)

bimodal_init$f_init_reporter = function(bimodal_parameters, seed) {
  set.seed(seed)
  gene_effect_sizes = sample(effect_distribution, bimodal_parameters$n_gene_effects, replace = TRUE)
  guide_effect_sizes = sapply(1:bimodal_parameters$n_gene_effects,
    function(i) {
      es = rnorm(
        bimodal_parameters$n_guides_per_target,
        mean = gene_effect_sizes[i],
        sd = bimodal_parameters$sigma_guide)
      correct_sign = sign(gene_effect_sizes[i])
      if (!all(sign(es) == correct_sign)) {
        es = abs(es) * correct_sign
      }
      es
    })

  n_null = with(bimodal_parameters, n_guides - n_gene_effects * n_guides_per_target)
  guide_effects = data.frame(
    effect_left = c(guide_effect_sizes, rep(0, n_null)), effect_right = NA)
  null_distribution = data.frame(
    ss_left_mean = rep(0, bimodal_parameters$n_replicates),
    ss_right_mean = rep(NA, bimodal_parameters$n_replicates),
    ss_left_sd = rep(1, bimodal_parameters$n_replicates),
    ss_right_sd = rep(NA, bimodal_parameters$n_replicates)
    )

  list(guide_effects = guide_effects, null_distribution = null_distribution)
}

library('Rcpp')

sourceCpp('bimodal_simulation.cpp')

set.seed(seed)

sim = bimodal_simulation(bimodal_parameters, bimodal_init)

source('dm_simple.R')

counts_to_dirichlet = function(counts, phi, eps = 0.01) {
  print(dim(counts))
  totals_by_replicate = apply(counts, c(1, 2), sum)
  proportion = counts
  sim = counts
  for (i in 1:dim(counts)[1]) {
    totals = matrix(rep(totals_by_replicate[i, ], times = dim(counts)[3]), ncol = dim(counts)[3])
    proportion[i, , ] = counts[i, , ] / totals

    # print(totals)
    for (j in 1:dim(counts)[2]) {
      # browser()
      alpha = proportion[i, j, ] * phi
      if (any(is.na(alpha))) {
        alpha[which(is.na(alpha))] = eps
      }
      if (any(alpha == 0)) {
        # NB: this is to prevent issues from sampling when alpha == 0
        alpha = alpha + eps
      }
      sim[i, j, ] = rdirmnom(1, totals_by_replicate[i, j], alpha)
    }
  }
  list(proportion = proportion, sim_counts = sim, totals = totals_by_replicate)
}

print(paste0('n_bins: ', length(bimodal_parameters$bins[[1]])))
print(bimodal_parameters$bins)
sim$counts = sim$guide_counts[, , 1:length(bimodal_parameters$bins[[1]])]

# debugonce(counts_to_dirichlet)
dm_counts = counts_to_dirichlet(sim$counts, dispersion)

sim$counts = dm_counts$sim_counts
sim_nimble = cell_simulation_to_nimble(bimodal_parameters, sim)


n_model = nimbleModel(gg_code,
  data = sim_nimble$data,
  constants = sim_nimble$const,
  inits = sim_nimble$init)

n_configuration = configureMCMC(n_model)
n_configuration$addMonitors(
  c('gene_inclusion',
    'sample_dispersion',
    'gene_shift',
    'guide_shift',
    'total_shift',
    'bin_alpha'))

n_mcmc_build = buildMCMC(n_configuration)
C_n_model = compileNimble(n_model)
C_n_mcmc = compileNimble(n_mcmc_build, project = n_model)

chain_seed = sample.int(.Machine$integer.max, 1)

total_time = system.time({
samples_n = runMCMC(C_n_mcmc,
  niter = n_mcmc_samples,
  nburnin = n_mcmc_warmup,
  nchains = n_mcmc_chains,
  thin = 1,
  setSeed = chain_seed,
  summary = TRUE)
})

saveRDS(samples_n, file = paste0('results/nimble/', fbase, '.rds'))
saveRDS(samples_n$summary, file = paste0('results/nimble/', fbase, '_summary.rds'))
saveRDS(sim, file = paste0('results/truth/sim_', fbase, '.rds'))
saveRDS(sim_nimble, file = paste0('results/truth/nimble_', fbase, '.rds'))
