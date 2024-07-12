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
library('ggplot2')
library('cowplot')
library('ggrepel')
theme_set(theme_cowplot())


########################################################################
# munge the low coverage data

raw_counts = read.table('../data/IL2RA_high_moi_screen_2021-02-04.count.txt',
  stringsAsFactors = FALSE, header = TRUE, row.names = 1,
  sep = '\t')

select_string = 'MOI5'
moi5 = dplyr::select(raw_counts, Gene, contains(select_string))
tmp_names = grep(select_string, colnames(moi5), value = TRUE)
data.table::setnames(moi5, tmp_names, sub(paste0('_', select_string), '', tmp_names))

ordering = c('Low', 'Low_Mid', 'Mid_High', 'High')

moi5_by_donor = lapply(paste0('D', c(1, 2)),
  function(d) {
    tmp = dplyr::select(moi5, starts_with(d))
    tmp_names = colnames(tmp)
    data.table::setnames(tmp, tmp_names, sub(paste0(d, '_'), '', tmp_names))
    as.matrix(tmp[, ordering])
  })

moi5_array = array(NA, dim = c(length(moi5_by_donor),
    nrow(moi5_by_donor[[1]]),
    ncol(moi5_by_donor[[1]])))

for (i in 1:length(moi5_by_donor)) {
  moi5_array[i, , ] = (moi5_by_donor[[i]])
}

dimnames(moi5_array) = list(NULL,
  rownames(moi5_by_donor[[1]]),
  colnames(moi5_by_donor[[1]]))

moi5_genes = data.frame(guide = rownames(raw_counts), gene = raw_counts$Gene)


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

debugonce(counts_to_nimble)
library('waterbear')
tmp = counts_to_nimble(moi5_array, moi5_genes)

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

saveRDS(samples_gg1, file = 'moi5.rds')

########################################################################

select_string = 'MOI2'
moi2 = dplyr::select(raw_counts, Gene, contains(select_string))
tmp_names = grep(select_string, colnames(moi2), value = TRUE)
data.table::setnames(moi2, tmp_names, sub(paste0('_', select_string), '', tmp_names))

ordering = c('Low', 'Low_Mid', 'Mid_High', 'High')

moi2_by_donor = lapply(paste0('D', c(1, 2)),
  function(d) {
    tmp = dplyr::select(moi2, starts_with(d))
    tmp_names = colnames(tmp)
    data.table::setnames(tmp, tmp_names, sub(paste0(d, '_'), '', tmp_names))
    as.matrix(tmp[, ordering])
  })

moi2_array = array(NA, dim = c(length(moi2_by_donor),
    nrow(moi2_by_donor[[1]]),
    ncol(moi2_by_donor[[1]])))

for (i in 1:length(moi2_by_donor)) {
  moi2_array[i, , ] = (moi2_by_donor[[i]])
}

dimnames(moi2_array) = list(NULL,
  rownames(moi2_by_donor[[1]]),
  colnames(moi2_by_donor[[1]]))

moi2_genes = data.frame(guide = rownames(raw_counts), gene = raw_counts$Gene)

tmp = counts_to_nimble(moi2_array, moi2_genes)

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
  samples_moi2 = runMCMC(C_gg_mcmc, niter = n_samples, nburnin = n_burnin,
    nchains = n_chains,
    thin = thin,
    setSeed = cur_seed,
    summary = TRUE)
})

saveRDS(samples_moi2, file = 'moi2.rds')
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

normal_coverage_by_donor = lapply(paste0('D', c(1, 2, 3)),
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

shift_low = low_summary[grepl('gene_shift', rownames(low_summary)), ]
shift_low = data.frame(shift_low)
shift_low = mutate(shift_low, mapping = rownames_to_index(rownames(shift_low)))
shift_low = inner_join(shift_low, gene_nimble_mapping, by = 'mapping')
shift_low = dplyr::select(shift_low, low_gene_shift = Mean, gene)

normal_summary = samples_normal$summary$all.chains
gi_normal = normal_summary[grepl('gene_inclusion', rownames(normal_summary)), ]
gi_normal = data.frame(gi_normal)
gi_normal = mutate(gi_normal, mapping = rownames_to_index(rownames(gi_normal)))
gi_normal = inner_join(gi_normal, gene_nimble_mapping, by = 'mapping')


fdr = 0.10
ppi = 1 - fdr
normal_fdr_10 = gi_normal[gi_normal[, 'Mean'] > ppi, ]

validated_targets = read.csv('../data/cd25_validated_targets.csv', header = TRUE,
  stringsAsFactors = FALSE)

validated_targets = dplyr::mutate(validated_targets, validated = TRUE)


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
ppi_merge = inner_join(ppi_merge, shift_low, by = 'gene')
ppi_merge = inner_join(ppi_merge, shift_normal, by = 'gene')

ppi_merge = left_join(ppi_merge, validated_targets, by = c('gene' = 'Gene_knocked_out'))
ppi_merge = dplyr::mutate(ppi_merge, validated = ifelse(is.na(validated), FALSE, validated))
ppi_merge = dplyr::mutate(ppi_merge, validated_gene = ifelse(validated, gene, NA))

coverage_validated = dplyr::filter(ppi_merge, validated)
alpha = 0.10
summarize(coverage_validated,
  x500_sensitivity = mean(1 - normal_mean < alpha),
  x500_tp = sum(1 - normal_mean < alpha),
  x50_sensitivity = mean(1 - low_mean < alpha),
  x50_tp = sum(1 - low_mean < alpha)
  )  %>%
pivot_longer(cols = everything())

ppi_merge = dplyr::mutate(ppi_merge,
  status = ppi_to_factor(low_mean, normal_mean, ppi, '50x', '500x'))
ppi_merge = dplyr::mutate(ppi_merge,
  gene_label = ifelse(status != 'neither', gene, ''))

p = ggplot(ppi_merge, aes(low_gene_shift, normal_gene_shift, color = status))
p = p + geom_vline(xintercept = 0, linetype = 2, alpha = 0.2)
p = p + geom_hline(yintercept = 0, linetype = 2, alpha = 0.2)
p = p + geom_point(alpha = 0.3)
p = p + geom_text_repel(aes(label = validated_gene), max.overlaps = Inf)
p = p + xlab('50x coverage')
p = p + ylab('500x coverage')
save_plot('50x_500x_validated.pdf', p, base_width = 10, base_height = 10)


shift_normal = normal_summary[grepl('gene_shift', rownames(normal_summary)), ]
shift_normal = data.frame(shift_normal)
shift_normal = mutate(shift_normal, mapping = rownames_to_index(rownames(shift_normal)))
shift_normal = inner_join(shift_normal, gene_nimble_mapping, by = 'mapping')
shift_normal = dplyr::select(shift_normal, normal_gene_shift = Mean, gene)

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

# compare moi
moi5_samples = readRDS('moi5.rds')


moi5_summary = moi5_samples$summary$all.chains
gi_moi5 = moi5_summary[grepl('gene_inclusion', rownames(moi5_summary)), ]
gi_moi5 = data.frame(gi_moi5)
gi_moi5 = mutate(gi_moi5, mapping = rownames_to_index(rownames(gi_moi5)))
gi_moi5 = inner_join(gi_moi5, gene_nimble_mapping, by = 'mapping')
gi_moi5 = mutate(gi_moi5, mcmc_moi5_fdr = 1 - Mean)

mean(gi_moi5$Mean > 0.90)

moi2_samples = readRDS('moi2.rds')

moi2_summary = moi2_samples$summary$all.chains
gi_moi2 = moi2_summary[grepl('gene_inclusion', rownames(moi2_summary)), ]
gi_moi2 = data.frame(gi_moi2)
gi_moi2 = mutate(gi_moi2, mapping = rownames_to_index(rownames(gi_moi2)))
gi_moi2 = inner_join(gi_moi2, gene_nimble_mapping, by = 'mapping')
gi_moi2 = mutate(gi_moi2, mcmc_moi2_fdr = 1 - Mean)

#
#

all_moi = inner_join(
  dplyr::select(gi_moi2, gene, mcmc_moi2_fdr),
  dplyr::select(gi_moi5, gene, mcmc_moi5_fdr),
  by = 'gene')

shift_moi5 = moi5_summary[grepl('gene_shift', rownames(moi5_summary)), ]
shift_moi5 = data.frame(shift_moi5)
shift_moi5 = mutate(shift_moi5, mapping = rownames_to_index(rownames(shift_moi5)))
shift_moi5 = inner_join(shift_moi5, gene_nimble_mapping, by = 'mapping')
shift_moi5 = dplyr::select(shift_moi5, moi5_gene_shift = Mean, gene)

shift_moi2 = moi2_summary[grepl('gene_shift', rownames(moi2_summary)), ]
shift_moi2 = data.frame(shift_moi2)
shift_moi2 = mutate(shift_moi2, mapping = rownames_to_index(rownames(shift_moi2)))
shift_moi2 = inner_join(shift_moi2, gene_nimble_mapping, by = 'mapping')
shift_moi2 = dplyr::select(shift_moi2, moi2_gene_shift = Mean, gene)

all_moi = inner_join(all_moi, shift_moi2, by = 'gene')
all_moi = inner_join(all_moi, shift_moi5, by = 'gene')

all_moi = mutate(all_moi,
  grouping = ppi_to_factor(1 - mcmc_moi2_fdr, 1 - mcmc_moi5_fdr, 0.9, 'mcmc_moi2', 'mcmc_moi5') )


p = ggplot(all_moi, aes(moi2_gene_shift, moi5_gene_shift, color = grouping))
p = p + geom_point(alpha = 0.3)
save_plot('mcmc_moi2_moi5.pdf', p, base_width = 10, base_height = 10)

# load the other mageck results

mageck_moi2 = read.table(
  '../data/high_moi_screen/moi2_results/IL2RA_moi2_low_high_pos_enrichment_2021-02-04.gene_summary.txt',
  header = TRUE, stringsAsFactors = FALSE)


mageck_moi5 = read.table(
  '../data/high_moi_screen/moi5_results/IL2RA_moi5_low_high_pos_enrichment_2021-02-04.gene_summary.txt',
  header = TRUE, stringsAsFactors = FALSE)
mageck_moi5 = dplyr::mutate(mageck_moi5, min_fdr = pmin(neg.fdr, pos.fdr))
mageck_moi5 = dplyr::select(mageck_moi5, gene = id, lfc = pos.lfc, min_fdr)

mageck_standard = read.table(
  '../data/CD25_D1_D2_D3_low_high_pos_enrichment_2019-04-23.gene_summary.txt',
  header = TRUE, stringsAsFactors = FALSE)
mageck_standard = dplyr::mutate(mageck_standard, min_fdr = pmin(neg.fdr, pos.fdr))
mageck_standard = dplyr::select(mageck_standard, gene = id, standard_lfc = pos.lfc, standard_min_fdr = min_fdr)

all_moi = inner_join(all_moi, mageck_standard, by = 'gene')

all_moi = mutate(all_moi,
  mageck_5_grouping = ppi_to_factor(1 - mcmc_moi5_fdr, 1 - standard_min_fdr, 0.9, 'mcmc_moi5', 'MAGeCK'),
  mageck_5_text = ifelse(mageck_5_grouping != 'neither', gene, ''))

p = ggplot(all_moi, aes(standard_lfc, moi5_gene_shift, color = mageck_5_grouping))
p = p + geom_point(alpha = 0.7)
p = p + geom_text_repel(aes(label = mageck_5_text), alpha = 0.7, max.overlaps = 50)
save_plot('mcmc_moi5_mageck.pdf', p, base_width = 10, base_height = 10)

# compare high 500x to mageck
gi_500 = dplyr::mutate(gi_500, c500_fdr = 1 - Mean)
gi_500 = dplyr::select(gi_500, gene, c500_fdr)

shift_500 = dplyr::select(shift_500, gene, shift_500 = Mean)

all_moi = inner_join(
  all_moi,
  gi_500, by = 'gene'
  )

all_moi = inner_join(
  all_moi,
  shift_500, by = 'gene')

all_moi = mutate(all_moi,
  mageck_500_grouping = ppi_to_factor(1 - c500_fdr, 1 - standard_min_fdr, 0.9, 'mcmc_500x', 'MAGeCK'),
  mageck_500_text = ifelse(mageck_500_grouping != 'neither', gene, '')
  )

p = ggplot(all_moi, aes(standard_lfc, shift_500, color = mageck_500_grouping))
p = p + geom_point(alpha = 0.7)
p = p + geom_text_repel(aes(label = mageck_500_text), alpha = 0.7, max.overlaps = 50)
save_plot('mcmc_500x_mageck.pdf', p, base_width = 10, base_height = 10)

all_moi = mutate(all_moi,
  moi5_500x_grouping = ppi_to_factor(1 - mcmc_moi5_fdr, 1 - c500_fdr, 0.9, 'mcmc_moi5', 'mcmc_500x'),
  moi5_500x_text = ifelse(moi5_500x_grouping != 'neither', gene, '')
  )

p = ggplot(all_moi, aes(shift_500, moi5_gene_shift, color = moi5_500x_grouping))
p = p + geom_point(alpha = 0.7)
p = p + geom_text_repel(aes(label = moi5_500x_text), alpha = 0.7, max.overlaps = 50)
save_plot('mcmc_500x_moi5.pdf', p, base_width = 10, base_height = 10)

table(all_moi$moi5_500x_grouping)

########################################################################
# newer data analysis
########################################################################

raw_counts = read.table('../data/IL2RA_high_moi_screen_2021-04-15.count.txt', header = TRUE,
  stringsAsFactors = FALSE, sep = '\t')

ordering = c('Low', 'Low_Mid', 'Mid_High', 'High')

c_name = grep('^D', colnames(raw_counts), value = TRUE)
sample_mapping = data.frame(c_name = c_name)
sample_mapping = mutate(sample_mapping, sample = sub('_.*$', '', c_name))
sample_mapping = mutate(sample_mapping, bin = sub('D._MOI._', '', c_name))
sample_mapping = mutate(sample_mapping, moi = sub('_.*$', '', sub('D._MOI', '', c_name)))

moi5_sample_mapping = dplyr::filter(sample_mapping, moi == 5)

moi5_array = counts_to_array(raw_counts, moi5_sample_mapping, ordering)

gene_mapping = raw_counts[, c('sgRNA', 'Gene')]
colnames(gene_mapping) = c('guide', 'gene')

wb = counts_to_wb(moi5_array, gene_mapping, c(0.15, 0.35, 0.35, 0.15))

gg_model = nimbleModel(gg_code, data = wb$data, constants = wb$const, inits = wb$init)

gg_configuration = configureMCMC(gg_model)

# gg_configuration$monitors
gg_configuration$addMonitors(c('gene_inclusion', 'total_shift', 'guide_shift', 'gene_shift',
    'sample_dispersion'))

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

samples_gg1$summary

saveRDS(samples_gg1, file = 'high_moi5_05_07_21.rds')

# compare MAGeCK and water bear
gene_nimble_mapping = dplyr::select(wb$test_guide_names, gene, mapping)
gene_nimble_mapping = distinct(gene_nimble_mapping)

moi5_samples = readRDS('high_moi5_05_07_21.rds')

moi5_summary = moi5_samples$summary$all.chains
gi_moi5 = moi5_summary[grepl('gene_inclusion', rownames(moi5_summary)), ]
gi_moi5 = data.frame(gi_moi5)
gi_moi5 = mutate(gi_moi5, mapping = rownames_to_index(rownames(gi_moi5)))
gi_moi5 = inner_join(gi_moi5, gene_nimble_mapping, by = 'mapping')
gi_moi5 = mutate(gi_moi5, mcmc_moi5_fdr = 1 - Mean)

shift_moi5 = moi5_summary[grepl('gene_shift', rownames(moi5_summary)), ]
shift_moi5 = data.frame(shift_moi5)
shift_moi5 = mutate(shift_moi5, mapping = rownames_to_index(rownames(shift_moi5)))
shift_moi5 = inner_join(shift_moi5, gene_nimble_mapping, by = 'mapping')
shift_moi5 = dplyr::select(shift_moi5, moi5_gene_shift = Mean, gene)

wb_moi5 = inner_join(
  dplyr::select(gi_moi5, gene, mcmc_moi5_fdr),
  shift_moi5,
  by = 'gene')

mageck_moi5 = read.table('../data/IL2RA_moi5_low_high_pos_enrichment_2021-04-15.gene_summary.txt',
  header = TRUE, sep = '\t', stringsAsFactors = FALSE)

mageck_moi5 = dplyr::mutate(mageck_moi5, min_fdr = pmin(neg.fdr, pos.fdr))
mageck_moi5 = dplyr::select(mageck_moi5, gene = id, moi5_lfc = pos.lfc, moi5_min_fdr = min_fdr)

all_moi5 = inner_join(wb_moi5, mageck_moi5, by = 'gene')

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

all_moi5 = dplyr::mutate(all_moi5,
  status = ppi_to_factor(1 - moi5_min_fdr, 1 - mcmc_moi5_fdr, 0.90, 'mageck', 'wb_dispersion'))
all_moi5 = dplyr::mutate(all_moi5,
  gene_label = ifelse(status == 'both', gene, NA))

p = ggplot(all_moi5, aes(moi5_gene_shift, -moi5_lfc, color = status))
p = p + geom_vline(xintercept = 0, linetype = 2, alpha = 0.2)
p = p + geom_hline(yintercept = 0, linetype = 2, alpha = 0.2)
p = p + geom_point(alpha = 0.3)
p = p + geom_text_repel(aes(label = gene_label), max.overlaps = Inf)
save_plot('wb_mageck_moi5_ssd.pdf', p, base_width = 10, base_height = 10)


validated_targets = read.csv('../data/cd25_validated_targets.csv', header = TRUE,
  stringsAsFactors = FALSE)

validated_targets = dplyr::mutate(validated_targets, validated = TRUE)

all_moi5 = left_join(all_moi5, validated_targets, by = c('gene' = 'Gene_knocked_out'))
all_moi5 = dplyr::mutate(all_moi5, validated = ifelse(is.na(validated), FALSE, validated))

all_moi5 = dplyr::mutate(all_moi5, validated_gene = ifelse(validated, gene, NA))

p = ggplot(all_moi5, aes(moi5_gene_shift, -moi5_lfc, color = status))
p = p + geom_vline(xintercept = 0, linetype = 2, alpha = 0.2)
p = p + geom_hline(yintercept = 0, linetype = 2, alpha = 0.2)
p = p + geom_point(alpha = 0.3)
p = p + geom_text_repel(aes(label = validated_gene), max.overlaps = Inf)
save_plot('wb_mageck_moi5_ssd_validated.pdf', p, base_width = 10, base_height = 10)



dplyr::filter(all_moi5, gene == 'IRF4' | gene == 'GATA3')

########################################################################
# newer data analysis, shared dispersion
########################################################################

raw_counts = read.table('../data/IL2RA_high_moi_screen_2021-04-15.count.txt', header = TRUE,
  stringsAsFactors = FALSE, sep = '\t')

ordering = c('Low', 'Low_Mid', 'Mid_High', 'High')

c_name = grep('^D', colnames(raw_counts), value = TRUE)
sample_mapping = data.frame(c_name = c_name)
sample_mapping = mutate(sample_mapping, sample = sub('_.*$', '', c_name))
sample_mapping = mutate(sample_mapping, bin = sub('D._MOI._', '', c_name))
sample_mapping = mutate(sample_mapping, moi = sub('_.*$', '', sub('D._MOI', '', c_name)))

moi5_sample_mapping = dplyr::filter(sample_mapping, moi == 5)

moi5_array = counts_to_array(raw_counts, moi5_sample_mapping, ordering)

gene_mapping = raw_counts[, c('sgRNA', 'Gene')]
colnames(gene_mapping) = c('guide', 'gene')

wb = counts_to_wb(moi5_array, gene_mapping, c(0.15, 0.35, 0.35, 0.15))

gg_model = nimbleModel(gg_code_shared_dispersion, data = wb$data, constants = wb$const, inits = wb$init)

gg_configuration = configureMCMC(gg_model)

# gg_configuration$monitors
gg_configuration$addMonitors(c('gene_inclusion', 'total_shift', 'guide_shift', 'gene_shift',
    'dispersion'))

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

samples_gg1$summary

saveRDS(samples_gg1, file = 'high_moi5_shared_dispersion_05_07_21.rds')


moi5_shared_samples = readRDS('high_moi5_shared_dispersion_05_07_21.rds')

moi5_shared_summary = moi5_shared_samples$summary$all.chains
gi_moi5_shared = moi5_shared_summary[grepl('gene_inclusion', rownames(moi5_shared_summary)), ]
gi_moi5_shared = data.frame(gi_moi5_shared)
gi_moi5_shared = mutate(gi_moi5_shared, mapping = rownames_to_index(rownames(gi_moi5_shared)))
gi_moi5_shared = inner_join(gi_moi5_shared, gene_nimble_mapping, by = 'mapping')
gi_moi5_shared = mutate(gi_moi5_shared, mcmc_moi5_shared_fdr = 1 - Mean)

shift_moi5_shared = moi5_shared_summary[grepl('gene_shift', rownames(moi5_shared_summary)), ]
shift_moi5_shared = data.frame(shift_moi5_shared)
shift_moi5_shared = mutate(shift_moi5_shared, mapping = rownames_to_index(rownames(shift_moi5_shared)))
shift_moi5_shared = inner_join(shift_moi5_shared, gene_nimble_mapping, by = 'mapping')
shift_moi5_shared = dplyr::select(shift_moi5_shared, moi5_shared_gene_shift = Mean, gene)

wb_moi5_shared = inner_join(
  dplyr::select(gi_moi5_shared, gene, mcmc_moi5_shared_fdr),
  shift_moi5_shared,
  by = 'gene')

all_moi5 = inner_join(all_moi5, wb_moi5_shared, by = 'gene')

all_moi5 = dplyr::mutate(all_moi5,
  status_shared = ppi_to_factor(1 - mcmc_moi5_fdr, 1 - mcmc_moi5_shared_fdr, 0.90,
    'wb', 'wb_shared'))


p = ggplot(all_moi5, aes(moi5_gene_shift, moi5_shared_gene_shift, color = status_shared))
p = p + geom_vline(xintercept = 0, linetype = 2, alpha = 0.2)
p = p + geom_hline(yintercept = 0, linetype = 2, alpha = 0.2)
p = p + geom_point(alpha = 0.3)
p = p + geom_text_repel(aes(label = validated_gene), max.overlaps = Inf)
save_plot('wb_wbshared_moi5_ssd_validated.pdf', p, base_width = 10, base_height = 10)

# who gets more?

moi5_validated = dplyr::filter(all_moi5, validated)
alpha = 0.10
summarize(moi5_validated,
  wb_sensitivity = mean(mcmc_moi5_fdr < alpha),
  wb_tp = sum(mcmc_moi5_fdr < alpha),
  wb_shared_sensitivity = mean(mcmc_moi5_shared_fdr < alpha),
  wb_shared_tp = sum(mcmc_moi5_shared_fdr < alpha),
  mageck_sensitivity = mean(moi5_min_fdr < 0.10),
  mageck_tp = sum(moi5_min_fdr < 0.10),
  )  %>%
pivot_longer(cols = everything())


summarize(all_moi5,
  wb_total = sum(mcmc_moi5_fdr < alpha),
  wb_shared_total = sum(mcmc_moi5_shared_fdr < alpha),
  mageck_total = sum(moi5_min_fdr < alpha),
  wb_wbshared_spearman = cor(moi5_gene_shift, moi5_shared_gene_shift, method = 'spearman'),
  wb_mageck_spearman = cor(moi5_gene_shift, -moi5_lfc, method = 'spearman'),
  wbshared_mageck_spearman = cor(moi5_shared_gene_shift, -moi5_lfc, method = 'spearman')) %>%
  pivot_longer(cols = everything())


########################################################################
# now compare to the original screen
########################################################################

original_ssd = readRDS('original_screen_mcmc_ssd.rds')


original_ssd_summary = original_ssd$summary$all.chains
gi_original_ssd = original_ssd_summary[grepl('gene_inclusion', rownames(original_ssd_summary)), ]
gi_original_ssd = data.frame(gi_original_ssd)
gi_original_ssd = mutate(gi_original_ssd, mapping = rownames_to_index(rownames(gi_original_ssd)))
gi_original_ssd = inner_join(gi_original_ssd, gene_nimble_mapping, by = 'mapping')
gi_original_ssd = mutate(gi_original_ssd, mcmc_original_ssd_fdr = 1 - Mean)

shift_original_ssd = original_ssd_summary[grepl('gene_shift', rownames(original_ssd_summary)), ]
shift_original_ssd = data.frame(shift_original_ssd)
shift_original_ssd = mutate(shift_original_ssd, mapping = rownames_to_index(rownames(shift_original_ssd)))
shift_original_ssd = inner_join(shift_original_ssd, gene_nimble_mapping, by = 'mapping')
shift_original_ssd = dplyr::select(shift_original_ssd, original_ssd_gene_shift = Mean, gene)

wb_original_ssd =
  inner_join(
  dplyr::select(gi_original_ssd, gene, mcmc_original_ssd_fdr),
  shift_original_ssd,
  by = 'gene')

moi5_original = inner_join(all_moi5, wb_original_ssd, by = 'gene')

moi5_original = dplyr::mutate(moi5_original,
  status_wb_original = ppi_to_factor(1 - mcmc_original_ssd_fdr, 1 - mcmc_moi5_fdr, 0.90,
    'wb_original_data', 'wb_moi5'))

p = ggplot(moi5_original, aes(original_ssd_gene_shift, moi5_gene_shift, color = status_wb_original))
p = p + geom_vline(xintercept = 0, linetype = 2, alpha = 0.2)
p = p + geom_hline(yintercept = 0, linetype = 2, alpha = 0.2)
p = p + geom_point(alpha = 0.3)
p = p + geom_text_repel(aes(label = validated_gene), max.overlaps = Inf)
save_plot('wbmoi5_wb_original_validated.pdf', p, base_width = 10, base_height = 10)

summarize(filter(moi5_original, validated),
  wb_sensitivity = mean(mcmc_moi5_fdr < alpha),
  wb_tp = sum(mcmc_moi5_fdr < alpha),
  wb_original_sensitivity = mean(mcmc_original_ssd_fdr < alpha),
  wb_original_tp = sum(mcmc_original_ssd_fdr < alpha),
  mageck_original_sensitivity = mean(standard_min_fdr < alpha),
  mageck_original_tp = sum(standard_min_fdr < alpha),
  )  %>%
pivot_longer(cols = everything())

mageck_standard = read.table(
  '../data/CD25_D1_D2_D3_low_high_pos_enrichment_2019-04-23.gene_summary.txt',
  header = TRUE, stringsAsFactors = FALSE)
mageck_standard = dplyr::mutate(mageck_standard, min_fdr = pmin(neg.fdr, pos.fdr))
mageck_standard = dplyr::select(mageck_standard, gene = id, standard_lfc = pos.lfc, standard_min_fdr = min_fdr)


moi5_original = inner_join(moi5_original, mageck_standard, by = 'gene')

moi5_original = dplyr::mutate(moi5_original,
  status_wb_original_mageck = ppi_to_factor(1 - mcmc_original_ssd_fdr, 1 - standard_min_fdr, 0.90,
    'wb_original_data', 'mageck_original'))

p = ggplot(moi5_original, aes(original_ssd_gene_shift, -standard_lfc, color = status_wb_original_mageck))
p = p + geom_vline(xintercept = 0, linetype = 2, alpha = 0.2)
p = p + geom_hline(yintercept = 0, linetype = 2, alpha = 0.2)
p = p + geom_point(alpha = 0.3)
p = p + geom_text_repel(aes(label = validated_gene), max.overlaps = Inf)
save_plot('wb_mageck_original.pdf', p, base_width = 10, base_height = 10)

moi5_original = dplyr::mutate(moi5_original,
  status_mageck = ppi_to_factor(1 - moi5_min_fdr, 1 - standard_min_fdr, 0.90,
    'moi5', 'original'))

p = ggplot(moi5_original, aes(moi5_lfc, standard_lfc, color = status_mageck))
p = p + geom_vline(xintercept = 0, linetype = 2, alpha = 0.2)
p = p + geom_hline(yintercept = 0, linetype = 2, alpha = 0.2)
p = p + geom_point(alpha = 0.3)
p = p + geom_text_repel(aes(label = validated_gene), max.overlaps = Inf)
save_plot('mageck_moi5_original.pdf', p, base_width = 10, base_height = 10)

summarize(moi5_original,
  wb_original_total = sum(mcmc_original_ssd_fdr < alpha),
  mageck_original_total = sum(standard_min_fdr < alpha),
  wb_original_moi5_spearman = cor(original_ssd_gene_shift, moi5_gene_shift, method = 'spearman'),
  mageck_moi5_original_spearman = cor(moi5_lfc, standard_lfc, method = 'spearman'),
  wb_original_mageck_original_spearman = cor(original_ssd_gene_shift, -standard_lfc, method = 'spearman'),
  wb_moi5_mageck_original_spearman = cor(moi5_gene_shift, -standard_lfc, method = 'spearman')
  )  %>%
pivot_longer(cols = everything())

# with new implementation

devtools::install('../waterbear')

library('waterbear')

raw_counts = read.table('../data/IL2RA_high_moi_screen_2021-04-15.count.txt', header = TRUE,
  stringsAsFactors = FALSE, sep = '\t')

ordering = c('Low', 'Low_Mid', 'Mid_High', 'High')

c_name = grep('^D', colnames(raw_counts), value = TRUE)
sample_mapping = data.frame(c_name = c_name)
sample_mapping = mutate(sample_mapping, sample = sub('_.*$', '', c_name))
sample_mapping = mutate(sample_mapping, bin = sub('D._MOI._', '', c_name))
sample_mapping = mutate(sample_mapping, moi = sub('_.*$', '', sub('D._MOI', '', c_name)))

moi5_sample_mapping = dplyr::filter(sample_mapping, moi == 5)

moi5_array = wb_counts_to_array(raw_counts, moi5_sample_mapping, ordering)

gene_mapping = raw_counts[, c('sgRNA', 'Gene')]
colnames(gene_mapping) = c('guide', 'gene')

wo = wb_make_object(moi5_array, gene_mapping, bin_size_prior = c(0.15, 0.35, 0.35, 0.15))

(counts, gene_mapping)


# debugonce(wb_em_start)
wo = wb_em_start(wo)
wo = wo$wo

source('../waterbear/R/nimble.R')
n_model = nimbleModel(sample_specific_dispersion_model,
  data = wo$data, constants = wo$const, inits = wo$init)

nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
nimbleOptions(MCMCsaveHistory = TRUE)

n_configuration = configureMCMC(n_model)
sampler_control = list(order = wo$const$order)
n_configuration$removeSamplers('gene_inclusion')
n_configuration$addSampler(target = 'gene_inclusion', type = 'rank_RW',
  control = sampler_control)

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
C_n_model = compileNimble(n_model)
C_n_mcmc = compileNimble(n_mcmc_build, project = n_model, showCompilerOutput = TRUE)


#     user   system  elapsed
# 5197.118   11.399 5220.962
n_samples = 10000
n_burnin = 5000
n_chains = 2
seed = 42
system.time({
samples = runMCMC(
  C_n_mcmc, niter = n_samples, nburnin = n_burnin,
  nchains = n_chains,
  setSeed = seed,
  summary = TRUE)
})

# > C_n_mcmc$samplerFunctions[[which_node]]$getAcceptanceHistory()
# [1] 1.000e+04 3.420e+02 4.539e+03 4.483e+03 9.780e+02 3.420e-02
which_node = n_configuration$findSamplersOnNodes('gene_inclusion')
C_n_mcmc$samplerFunctions[[which_node]]$getAcceptanceHistory()

rc = waterbear:::wb_recode(samples, wo)
rc %>% arrange(desc(Mean)) %>% head(30)

rc = mutate(rc, rank_fdr = 1 - Mean)

all_moi5 = inner_join(wb_moi5,
  dplyr::select(rc, gene, rank_fdr), by = 'gene')


rank_shift_moi5 = samples$summary$all.chain[grepl('gene_shift', rownames(moi5_summary)), ]
rank_shift_moi5 = data.frame(rank_shift_moi5)
rank_shift_moi5 = mutate(rank_shift_moi5, mapping = rownames_to_index(rownames(rank_shift_moi5)))
rank_shift_moi5 = inner_join(rank_shift_moi5, gene_nimble_mapping, by = 'mapping')
rank_shift_moi5 = dplyr::select(rank_shift_moi5, rank_gene_shift = Mean, gene)

rank_moi5 = inner_join(
  dplyr::select(rc, gene, rank_fdr),
  rank_shift_moi5,
  by = 'gene'
)


mcmc_moi5 = inner_join(rank_moi5, wb_moi5, by = 'gene')


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

mcmc_moi5 = dplyr::mutate(mcmc_moi5,
  status = ppi_to_factor(1 - rank_fdr, 1 - mcmc_moi5_fdr, 0.90, 'wb_rank', 'wb_dispersion'))
mcmc_moi5 = dplyr::mutate(mcmc_moi5,
  gene_label = ifelse(status == 'both', gene, NA))


p = ggplot(mcmc_moi5, aes(moi5_gene_shift, rank_gene_shift, color = status))
p = p + geom_vline(xintercept = 0, linetype = 2, alpha = 0.2)
p = p + geom_hline(yintercept = 0, linetype = 2, alpha = 0.2)
p = p + geom_point(alpha = 0.3)
p = p + geom_text_repel(aes(label = gene_label), max.overlaps = Inf)
save_plot('wb_rank_moi5_ssd.pdf', p, base_width = 10, base_height = 10)


mcmc_moi5 = left_join(mcmc_moi5, validated_targets, by = c('gene' = 'Gene_knocked_out'))
mcmc_moi5 = dplyr::mutate(mcmc_moi5, validated = ifelse(is.na(validated), FALSE, validated))

mcmc_moi5 = dplyr::mutate(mcmc_moi5, validated_gene = ifelse(validated, gene, NA))

# p = ggplot(mcmc_moi5, aes(moi5_gene_shift, -moi5_lfc, color = status))
p = ggplot(mcmc_moi5, aes(moi5_gene_shift, rank_gene_shift, color = status))
p = p + geom_vline(xintercept = 0, linetype = 2, alpha = 0.2)
p = p + geom_hline(yintercept = 0, linetype = 2, alpha = 0.2)
p = p + geom_point(alpha = 0.3)
p = p + geom_text_repel(aes(label = validated_gene), max.overlaps = Inf)
save_plot('wb_rank_moi5_ssd_validated.pdf', p, base_width = 10, base_height = 10)

mcmc_validated = dplyr::filter(mcmc_moi5, validated)
alpha = 0.10
summarize(mcmc_validated,
  wb_sensitivity = mean(mcmc_moi5_fdr < alpha),
  wb_tp = sum(mcmc_moi5_fdr < alpha),
  wb_rank_sensitivity = mean(rank_fdr < alpha),
  wb_rank_tp = sum(rank_fdr < alpha),
  )  %>%
pivot_longer(cols = everything())
