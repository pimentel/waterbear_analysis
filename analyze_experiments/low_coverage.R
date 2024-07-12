# when debugging
cur_seed = 326
n_samples = 11000
n_burnin = 1000
thin = 1

cur_seed = as.integer(argv[1])
n_samples = as.integer(argv[2])
n_burnin = as.integer(argv[3])
thin = as.integer(argv[4])

library('dplyr')
library('tidyr')
library('nimble')


########################################################################
# munge the low coverage data

raw_counts = read.table('../data/IL2RA_low_coverage_screen_2020-12-09.count.txt',
  stringsAsFactors = FALSE, header = TRUE, row.names = 1,
  sep = '\t')

low_coverage = dplyr::select(raw_counts, Gene, contains('50x'))
tmp_names = grep('50x', colnames(low_coverage), value = TRUE)
data.table::setnames(low_coverage, tmp_names, sub('_50x', '', tmp_names))

ordering = c('low', 'low_mid', 'mid_high', 'high')

low_coverage_by_donor = lapply(paste0('D', c(1, 3)),
  function(d) {
    tmp = dplyr::select(low_coverage, starts_with(d))
    tmp_names = colnames(tmp)
    data.table::setnames(tmp, tmp_names, sub(paste0(d, '_'), '', tmp_names))
    as.matrix(tmp[, ordering])
  })

low_coverage_array = array(NA, dim = c(length(low_coverage_by_donor),
    nrow(low_coverage_by_donor[[1]]),
    ncol(low_coverage_by_donor[[1]])))

for (i in 1:length(low_coverage_by_donor)) {
  low_coverage_array[i, , ] = (low_coverage_by_donor[[i]])
}

dimnames(low_coverage_array) = list(NULL,
  rownames(low_coverage_by_donor[[1]]),
  colnames(low_coverage_by_donor[[1]]))

low_coverage_genes = data.frame(guide = rownames(raw_counts), gene = raw_counts$Gene)


counts_to_nimble = function(
  counts_array,
  gene_mapping,
  bin_size_prior = NULL
  ) {
  guide_names = data.frame(guide = dimnames(counts_array)[[2]])
  guide_names = dplyr::mutate(guide_names, i = 1:nrow(guide_names))
  N_bins = dim(counts_array)[3]

  if (is.null(bin_size_prior)) {
    bin_size_prior = rep(1 / N_bins, length.out = N_bins)
  }

  nt_guide_names = inner_join(guide_names,
    dplyr::filter(gene_mapping, grepl('Non-', gene)), by = 'guide')
  test_guide_names = dplyr::inner_join(guide_names,
    dplyr::filter(gene_mapping, !grepl('Non-', gene)), by = c('guide'))
  test_guide_names = dplyr::mutate(test_guide_names, mapping = as.integer(factor(gene)))
  test_guide_names = dplyr::arrange(test_guide_names, i)

  gg_data = list(x = counts_array)
  dispersion_init = 300
  gg_const = list(
    N = dim(gg_data$x)[1],
    N_guides = nrow(test_guide_names),
    x_total = apply(gg_data$x, c(1, 2), sum),
    N_genes = length(unique(test_guide_names$mapping)),
    N_nt = nrow(nt_guide_names),
    N_bins = N_bins,
    nt_index = 1:nrow(nt_guide_names),
    nt_data_index = nt_guide_names$i,
    guide_index = 1:nrow(test_guide_names),
    guide_data_index = test_guide_names$i,
    guide_to_gene = test_guide_names$mapping,
    dispersion_prior_mean = dispersion_init
    )

  gg_const$bin_alpha_prior = matrix(
    rep(bin_size_prior, gg_const$N),
      nrow = gg_const$N, byrow = TRUE)
  gg_const$N_cutoffs = gg_const$N_bins - 1

  gg_init = list(
    psi = 0.2,
    sigma_gene = 5,
    sigma_guide = 5,
    dispersion = dispersion_init
    # bin_alpha = matrix(rep(c(0.15, 0.7, 0.15), gg_const$N),
    #   nrow = gg_const$N, byrow = TRUE)
    )
  # gg_const$bin_alpha_prior = gg_init$bin_alpha
  gg_init$bin_alpha = gg_const$bin_alpha_prior
  gg_init$cutoffs = matrix(rep(dirichlet_to_normal_bins(gg_init$bin_alpha[1, ]), gg_const$N),
    nrow = gg_const$N, byrow = TRUE)
  q_init = array(0, dim = c(gg_const$N, gg_const$N_guides, N_bins))
  for (n in 1:dim(q_init)[1]) {
    for (g in 1:dim(q_init)[2]) {
      q_init[n, g, ] = bin_size_prior
    }
  }

  gg_init$q = q_init
  gg_init$gene_inclusion = rep(0, gg_const$N_genes)
  gg_init$gene_shift = rep(0, gg_const$N_genes)
  gg_init$guide_shift = rep(0, gg_const$N_guides)
  gg_init$total_shift = rep(0, gg_const$N_guides)

  list(data = gg_data, init = gg_init, const = gg_const, test_guide_names = test_guide_names)
}

# debugonce(counts_to_nimble)
tmp = counts_to_nimble(low_coverage_array, low_coverage_genes)

gg_model = nimbleModel(gg_code, data = tmp$data, constants = tmp$const, inits = tmp$init)

gg_configuration = configureMCMC(gg_model)

# gg_configuration$monitors
gg_configuration$addMonitors(c('gene_inclusion', 'total_shift', 'guide_shift', 'gene_shift'))

gg_mcmc_build = buildMCMC(gg_configuration)
C_gg_model = compileNimble(gg_model)
C_gg_mcmc = compileNimble(gg_mcmc_build, project = gg_model)

n_samples = 50000
n_burnin = 25000
n_chains = 4
thin = 1
cur_seed = 326

system.time({
  samples_gg1 = runMCMC(C_gg_mcmc, niter = n_samples, nburnin = n_burnin,
    nchains = n_chains,
    thin = thin,
    setSeed = cur_seed,
    summary = TRUE)
})

saveRDS(samples_gg1, file = '50x.rds')


########################################################################

counts_array = readRDS('../data/cd25_counts_array.rds')
gene_mapping = readRDS('../data/gene_mapping.rds')


ddirchmulti <- nimbleFunction(
  run = function(x = double(1), alpha = double(1), size = double(0),
    log = integer(0, default = 0)) {
    returnType(double(0))
    logProb <- lgamma(size + 1) - sum(lgamma(x + 1)) + lgamma(sum(alpha)) -
      sum(lgamma(alpha)) + sum(lgamma(alpha + x)) - lgamma(sum(alpha) +
        size)
    if(log) return(logProb)
    else return(exp(logProb))
  })

rdirchmulti <- nimbleFunction(
  run = function(n = integer(0), alpha = double(1), size = double(0)) {
    returnType(double(1))
    if(n != 1) print("rdirchmulti only allows n = 1; using n = 1.")
    p <- rdirch(1, alpha)
    return(rmulti(1, size = size, prob = p))
  })

dirichlet_to_normal_bins = nimbleFunction(
  run = function(alpha = double(1)) {
    returnType(double(1))
    a0 = sum(alpha)
    p = alpha / a0
    # qp = double(1)
    qp = numeric(length(p))
    qp[1] = p[1]
    for (i in 2:length(p)) {
      qp[i] = qp[i - 1] + p[i]
    }
    cutoff = numeric(length(qp) - 1)
    for (i in 1:(length(qp) - 1)) {
      cutoff[i] = probit(qp[i])
    }
    return(cutoff)
  })

cutoff_to_probability = nimbleFunction(
  run = function(cutoff = double(1), offset = double(0)) {
    returnType(double(1))
    # cutoff = c(-Inf, cutoff, Inf)
    # p = double(1)
    p = numeric(length(cutoff) + 1)
    p[1] = phi(cutoff[1] - offset)
    total = p[1]
    for (i in 2:length(cutoff)) {
      p[i] = phi(cutoff[i] - offset) - phi(cutoff[i - 1] - offset)
      total = total + p[i]
    }
    p[length(cutoff) + 1] = 1 - total
    return(p)
  })
gg_code = nimbleCode({
  dispersion ~ dexp(1.0 / dispersion_prior_mean)
  for (n in 1:N) {
    bin_alpha[n, 1:N_bins] ~ ddirch(bin_alpha_prior[n, 1:N_bins])
    cutoffs[n, 1:N_cutoffs] <- dirichlet_to_normal_bins(bin_alpha[n, 1:N_bins])
  }
  psi ~ dbeta(10, 10)
  sigma_gene ~ dgamma(1, 0.10)
  sigma_guide ~ dgamma(1, 0.10)
  for (gene in 1:N_genes) {
    gene_inclusion[gene] ~ dbern(psi)
    gene_shift[gene] ~ dnorm(0, sd = sigma_gene)
  }
  for (g in 1:N_nt) {
    for (n in 1:N) {
      x[n, nt_data_index[g], 1:N_bins] ~ ddirchmulti(dispersion * bin_alpha[n, 1:N_bins],
        x_total[n, nt_data_index[g]])
    }
  }
  for (g in 1:N_guides) {
    guide_shift[g] ~ dnorm(gene_shift[guide_to_gene[g]], sd = sigma_guide)
    total_shift[g] <- gene_inclusion[guide_to_gene[g]] * guide_shift[g]
    for (n in 1:N) {
      q[n, g, 1:N_bins] <- cutoff_to_probability(cutoffs[n, 1:N_cutoffs], total_shift[g])
      x[n, guide_data_index[g], 1:N_bins] ~ ddirchmulti(
        dispersion * q[n, g, 1:N_bins],
        x_total[n, guide_data_index[g]])
    }
  }
})


########################################################################
# normal coverage

normal_coverage = dplyr::select(raw_counts, Gene, contains('500x'))
tmp_names = grep('500x', colnames(normal_coverage), value = TRUE)
data.table::setnames(normal_coverage, tmp_names, sub('_500x', '', tmp_names))

ordering = c('low', 'low_mid', 'mid_high', 'high')

normal_coverage_by_donor = lapply(paste0('D', 1:3),
  function(d) {
    tmp = dplyr::select(normal_coverage, starts_with(d))
    tmp_names = colnames(tmp)
    data.table::setnames(tmp, tmp_names, sub(paste0(d, '_'), '', tmp_names))
    as.matrix(tmp[, ordering])
  })

normal_coverage_array = array(NA, dim = c(length(normal_coverage_by_donor),
    nrow(normal_coverage_by_donor[[1]]),
    ncol(normal_coverage_by_donor[[1]])))

for (i in 1:length(normal_coverage_by_donor)) {
  normal_coverage_array[i, , ] = (normal_coverage_by_donor[[i]])
}

dimnames(normal_coverage_array) = list(NULL,
  rownames(normal_coverage_by_donor[[1]]),
  colnames(normal_coverage_by_donor[[1]]))

normal_coverage_genes = data.frame(guide = rownames(raw_counts), gene = raw_counts$Gene)

debugonce(counts_to_nimble)
tmp = counts_to_nimble(normal_coverage_array, normal_coverage_genes)

gg_model = nimbleModel(gg_code, data = tmp$data, constants = tmp$const, inits = tmp$init)

gg_configuration = configureMCMC(gg_model)

# gg_configuration$monitors
gg_configuration$addMonitors(c('gene_inclusion', 'total_shift', 'guide_shift', 'gene_shift'))

gg_mcmc_build = buildMCMC(gg_configuration)
C_gg_model = compileNimble(gg_model)
C_gg_mcmc = compileNimble(gg_mcmc_build, project = gg_model)

n_samples = 50000
n_burnin = 25000
n_chains = 4
thin = 1
cur_seed = 326

system.time({
  samples_gg1 = runMCMC(C_gg_mcmc, niter = n_samples, nburnin = n_burnin,
    nchains = n_chains,
    thin = thin,
    setSeed = cur_seed,
    summary = TRUE)
})

saveRDS(samples_gg1, file = '500x.rds')

########################################################################
# now compare the 2

library('dplyr')
library('tidyr')
library('nimble')
library('ggplot2')
library('cowplot')
library('ggrepel')

theme_set(theme_cowplot())

samples_low = readRDS('50x.rds')
samples_normal = readRDS('500x.rds')


gene_nimble_mapping = dplyr::select(tmp$test_guide_names, gene, mapping)
gene_nimble_mapping = distinct(gene_nimble_mapping)

rownames_to_index = function(x) {
  x = sub('^[a-zA-Z_]+\\[', '', x)
  x = sub('\\]', '', x)
  as.integer(x)
}

low_summary = samples_low$summary$all.chains
gi_low = low_summary[grepl('gene_inclusion', rownames(low_summary)), ]
gi_low = data.frame(gi_low)
gi_low = mutate(gi_low, mapping = rownames_to_index(rownames(gi_low)))
gi_low = inner_join(gi_low, gene_nimble_mapping, by = 'mapping')

fdr = 0.10
ppi = 1 - fdr
low_fdr_10 = gi_low[gi_low[, 'Mean'] > ppi, ]

normal_summary = samples_normal$summary$all.chains
gi_normal = normal_summary[grepl('gene_inclusion', rownames(normal_summary)), ]
gi_normal = data.frame(gi_normal)
gi_normal = mutate(gi_normal, mapping = rownames_to_index(rownames(gi_normal)))
gi_normal = inner_join(gi_normal, gene_nimble_mapping, by = 'mapping')

gi_500 = gi_normal

fdr = 0.10
ppi = 1 - fdr
normal_fdr_10 = gi_normal[gi_normal[, 'Mean'] > ppi, ]

ppi_to_factor = function(x, y, ppi, label_x, label_y) {
  out = rep('neither', length(x))
  both = x >= ppi & y >= ppi
  x_only = x >= ppi & y < ppi
  y_only = y >= ppi & x < ppi
  out[both] = 'both'
  out[x_only] = label_x
  out[y_only] = label_y
  as.factor(out)
}

ppi_merge = inner_join(
  dplyr::select(gi_low, low_mean = Mean, low_sd = St.Dev., gene),
  dplyr::select(gi_normal, normal_mean = Mean, normal_sd = St.Dev., gene),
  by = 'gene'
  )


ppi_merge = dplyr::mutate(ppi_merge,
  status = ppi_to_factor(low_mean, normal_mean, ppi, '50x', '500x'))
ppi_merge = dplyr::mutate(ppi_merge,
  gene_label = ifelse(status != 'neither', gene, ''))

p = ggplot(ppi_merge, aes(low_mean, normal_mean, color = status))
p = p + geom_point(alpha = 0.3)
p = p + geom_text_repel(aes(label = gene_label), max.overlaps = Inf)
save_plot('low_normal_ppi.pdf', p, base_width = 10, base_height = 10)

# get old experiment

seed = 42:(42 + 4 - 1)
fnames = paste0('../results/mcmc/dmult_flexible_bins', seed, '_100000_50000_1.rds')


j = lapply(fnames,
  function(f) {
    x = readRDS(f)
    x = x$summary
    x = x[grepl('gene_inclusion', rownames(x)), 'Mean']
  })

k = simplify2array(j)
k = apply(k, 1, mean)
old_gene_summary = data.frame(gene_label = names(k), mean = k)
old_gene_summary = mutate(old_gene_summary, mapping = rownames_to_index(gene_label))


# do this to get the gene names
counts_array = readRDS('../data/cd25_counts_array.rds')
gene_mapping = readRDS('../data/gene_mapping.rds')

guide_names = data.frame(guide = dimnames(counts_array)[[2]])
guide_names = dplyr::mutate(guide_names, i = 1:nrow(guide_names))

nt_guide_names = inner_join(guide_names,
  dplyr::filter(gene_mapping, grepl('Non-', gene)), by = 'guide')
test_guide_names = dplyr::inner_join(guide_names,
  dplyr::filter(gene_mapping, !grepl('Non-', gene)), by = c('guide'))

test_guide_names = dplyr::mutate(test_guide_names, mapping = as.integer(factor(gene)))
test_guide_names = dplyr::arrange(test_guide_names, i)

old_gene_names = distinct(select(test_guide_names, gene, mapping))
old_gene_summary = inner_join(old_gene_summary, old_gene_names, by = 'mapping')
old_gene_summary = dplyr::select(old_gene_summary, gene, mean)


ppi_500_old = inner_join(
  select(gi_normal, normal_mean = Mean, gene),
  select(old_gene_summary, old_mean = mean, gene),
  by = 'gene')

ppi_500_old = mutate(ppi_500_old, status = ppi_to_factor(normal_mean, old_mean, 0.90, '500x', 'old'))
ppi_500_old = mutate(ppi_500_old, gene_label = ifelse(status != 'neither', gene, ''))

p = ggplot(ppi_500_old, aes(old_mean, normal_mean, color = status))
p = p + geom_point(alpha = 0.3)
p = p + geom_text_repel(aes(label = gene_label), max.overlaps = Inf)
save_plot('old_normal_ppi.pdf', p, base_width = 10, base_height = 10)

table(ppi_500_old$status)

# now looking at the effect sizes
shift_normal = normal_summary[grepl('gene_shift', rownames(normal_summary)), ]
shift_normal = data.frame(shift_normal)
shift_normal = mutate(shift_normal, mapping = rownames_to_index(rownames(shift_normal)))
shift_normal = inner_join(shift_normal, gene_nimble_mapping, by = 'mapping')

shift_500 = shift_normal

shift_low = low_summary[grepl('gene_shift', rownames(low_summary)), ]
shift_low = data.frame(shift_low)
shift_low = mutate(shift_low, mapping = rownames_to_index(rownames(shift_low)))
shift_low = inner_join(shift_low, gene_nimble_mapping, by = 'mapping')

shift_low_normal = inner_join(
  dplyr::select(shift_low, gene, low_mean = Mean),
  dplyr::select(shift_normal, gene, normal_mean = Mean),
  by = 'gene'
  )

shift_low_normal = inner_join(
  shift_low_normal,
  dplyr::select(ppi_merge, gene, gene_label, status),
  by = 'gene')

p = ggplot(shift_low_normal, aes(low_mean, normal_mean, color = status))
p = p + geom_point(alpha = 0.3)
p = p + geom_text_repel(aes(label = gene_label), max.overlaps = Inf)
save_plot('shift_low_normal.pdf', p, base_width = 10, base_height = 10)

j = lapply(fnames,
  function(f) {
    x = readRDS(f)
    x = x$summary
    x = x[grepl('gene_shift', rownames(x)), 'Mean']
  })
k = simplify2array(j)
k = apply(k, 1, mean)
old_gene_shift = data.frame(gene_label = names(k), mean = k)
old_gene_shift = mutate(old_gene_shift, mapping = rownames_to_index(gene_label))

old_gene_shift = inner_join(old_gene_shift, old_gene_names, by = 'mapping')
old_gene_shift = dplyr::select(old_gene_shift, gene, mean)

shift_old_normal = inner_join(
  dplyr::select(shift_normal, gene, normal_mean = Mean),
  dplyr::select(old_gene_shift, gene, old_mean = mean),
  by = 'gene')

shift_old_normal = dplyr::inner_join(
  shift_old_normal,
  dplyr::select(ppi_500_old, gene, gene_label, status),
  by = 'gene')

p = ggplot(shift_old_normal, aes(old_mean, normal_mean, color = status))
p = p + geom_point(alpha = 0.3)
p = p + geom_text_repel(aes(label = gene_label), max.overlaps = Inf)
save_plot('shift_old_normal.pdf', p, base_width = 10, base_height = 10)

# assumes ran some stuff from low coverage
