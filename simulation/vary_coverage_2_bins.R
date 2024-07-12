
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

args = c('42',
  '4',
  '250',
  '100',
  '0.10',
  '3',
  '10000',
  '5000',
  '1',
  '7333333',
  'easy',
  '200.00',
  '1')

# this is the basic case trying to debug
args = c(
  '2',
  '4',
  '1000',
  '100',
  '0.10',
  # '0.01',
  '2',
  '50000',
  '25000',
  '2',
  '410000',
  'easy',
  '200.00',
  '4'
)

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
  'vary_coverage_2_bins',
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

counts_to_dirichlet = function(counts, phi) {
  print(dim(counts))
  totals_by_replicate = apply(counts, c(1, 2), sum)
  proportion = counts
  sim = counts
  for (i in 1:dim(counts)[1]) {
    totals = matrix(rep(totals_by_replicate[i, ], times = dim(counts)[3]), ncol = dim(counts)[3])
    proportion[i, , ] = counts[i, , ] / totals
    for (j in 1:dim(counts)[2]) {
      alpha = proportion[i, j, ] * phi
      sim[i, j, ] = rdirmnom(1, totals_by_replicate[i, j], alpha)
    }
  }
  list(proportion = proportion, sim_counts = sim, totals = totals_by_replicate)
}

print(paste0('n_bins: ', length(bimodal_parameters$bins[[1]])))
print(bimodal_parameters$bins)
sim$counts = sim$guide_counts[, , 1:length(bimodal_parameters$bins[[1]])]

which_replace = which(sim$counts == 0)
sim$counts[which_replace] = 1

dm_counts = counts_to_dirichlet(sim$counts, dispersion)

sim$counts = dm_counts$sim_counts

sim_nimble = cell_simulation_to_nimble(bimodal_parameters, sim)

gene_mapping = data.frame(
  guide = as.character(c(sim_nimble$const$guide_data_index, sim_nimble$const$nt_data_index)),
  gene = c(as.character(sim_nimble$const$guide_to_gene), rep('Non-targeting', sim_nimble$const$N_nt))
)

guide_counts = apply(sim$counts, c(1, 2), sum)
guide_counts = lapply(1:nrow(guide_counts),
  function(i) {
    as.integer(guide_counts[i, ])
  })

tmp_counts = sim$counts[, , c(1, 4)]
dimnames(tmp_counts)[[2]] = gene_mapping$guide
dimnames(tmp_counts)[[3]] = c('low', 'high')

library('waterbear')

# debugonce(wb_estimate_unobserved)

# debugonce(waterbear:::estimate_unobserved_per_sample)

observed_bin_sizes = c(bins[[1]][[1]][2] - bins[[1]][[1]][1],
  bins[[1]][[4]][2] - bins[[1]][[4]][1])


complete_array = wb_estimate_unobserved(tmp_counts, guide_counts, observed_bin_sizes, c('low', 'unobserved', 'high'))

# bin_size_prior = c(
#   bins[[1]][[1]][2] - bins[[1]][[1]][1],
#   # bins[[1]][[4]][1] - bins[[1]][[1]][2],
#   bins[[1]][[4]][2] - bins[[1]][[4]][1])
# bin_size_prior = bin_size_prior / sum(bin_size_prior)
# complete_array = tmp_counts

bin_size_prior = c(
  bins[[1]][[1]][2] - bins[[1]][[1]][1],
  bins[[1]][[4]][1] - bins[[1]][[1]][2],
  bins[[1]][[4]][2] - bins[[1]][[4]][1])

wo = wb_make_object(complete_array, gene_mapping, bin_size_prior = bin_size_prior)

# debugonce(wb_em_start)
# debugonce(waterbear:::estimate_bin_sizes)
# debugonce(waterbear:::estimate_guide_mu)
wo = wb_em_start(wo)
# wo = wo$wo


saveRDS(wo, file = paste0('results/truth/wo_', fbase, '.rds'))

source('../waterbear/R/nimble.R')
n_model = nimbleModel(sample_specific_dispersion_model,
  data = wo$data, constants = wo$const, inits = wo$init)

nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
nimbleOptions(MCMCsaveHistory = TRUE)

n_configuration = configureMCMC(n_model)
# sampler_control = list(order = wo$const$order)
# n_configuration$removeSamplers('gene_inclusion')
# n_configuration$addSampler(target = 'gene_inclusion', type = 'rank_RW',
#   control = sampler_control)

n_configuration$addMonitors(
  c('gene_inclusion',
    'total_shift',
    'guide_shift',
    'gene_shift',
    'dispersion',
    'sample_dispersion',
    'psi'
    ))

n_mcmc_build = buildMCMC(n_configuration)
C_n_model = compileNimble(n_model, showCompilerOutput = TRUE)
C_n_mcmc = compileNimble(n_mcmc_build, project = n_model, showCompilerOutput = TRUE)

chain_seed = sample.int(.Machine$integer.max, 1)

#     user   system  elapsed
# 5197.118   11.399 5220.962
# n_samples = 1000
# n_burnin = 500
# n_mcmc_chains = 1
# seed = 42
system.time({
samples_n = runMCMC(
  C_n_mcmc,
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

# # TODO: practice
# gi = waterbear:::wb_recode(samples_n, wo)
# complete_array[1, 300, ]
