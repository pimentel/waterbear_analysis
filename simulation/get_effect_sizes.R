library('dplyr')
library('ggplot2')
library('cowplot')
library('tidyr')

theme_set(theme_cowplot())

fnames = Sys.glob('../data/coverage_comparison/mageck_results/*.gene*')
all_results = lapply(fnames, read.table, stringsAsFactors = FALSE,
  header = TRUE)

which_coverage = basename(fnames)
which_coverage = sub('coverage_test_', '', which_coverage)
which_coverage = sub('_.+$', '', which_coverage)

alpha = 0.05
significant_results = lapply(seq_along(all_results),
  function(i) {
    x = all_results[[i]]
    x = dplyr::filter(x, neg.fdr < alpha | pos.fdr < alpha)
    x = dplyr::mutate(x, coverage = which_coverage[i])
  })
names(significant_results) = which_coverage

coverage_comparison = lapply(2:length(significant_results),
  function(i) {
    x = significant_results[[i]]
    list(
      both = intersect(significant_results[[1]]$id, x$id),
      upper_only = setdiff(significant_results[[1]]$id, x$id),
      lower_only = setdiff(x$id, significant_results[[1]]$id),
      comparison = which_coverage[[i]]
      )
  })

# these are the ones that you can get at every coverage level including the highest coverage
easy = Reduce(
  function(x, y) {
    z = list()
    z$both = intersect(x$both, y$both)
    z
  }, coverage_comparison)
easy = easy$both

# the ones that you can get at high coverage and 200x coverage but not 50x
moderate = intersect(significant_results[['1000x']]$id, significant_results[['200x']]$id)
moderate = setdiff(moderate, significant_results[['50x']]$id)

# the ones that you get add high coverage but not at 200x coverage
difficult = setdiff(significant_results[['1000x']]$id, significant_results[['200x']]$id)
difficult

counts_array = readRDS('../data/cd25_counts_array.rds')
gene_mapping = readRDS('../data/gene_mapping.rds')

# now look at the CD25 data
fnames = Sys.glob('../results/mcmc/dmult_flexible_bins*_100000_50000_*.rds')

guide_names = data.frame(guide = dimnames(counts_array)[[2]])
guide_names = dplyr::mutate(guide_names, i = 1:nrow(guide_names))

nt_guide_names = inner_join(guide_names,
  dplyr::filter(gene_mapping, grepl('Non-', gene)), by = 'guide')
test_guide_names = dplyr::inner_join(guide_names,
  dplyr::filter(gene_mapping, !grepl('Non-', gene)), by = c('guide'))

test_guide_names = dplyr::mutate(test_guide_names, mapping = as.integer(factor(gene)))
test_guide_names = dplyr::arrange(test_guide_names, i)

# after this, select the genes that are in the different categories and then pull the effect sizes from
# the summaries
gene_index_mapping = dplyr::select(test_guide_names, gene, mapping)
gene_index_mapping = distinct(gene_index_mapping)
gene_index_mapping = dplyr::mutate(gene_index_mapping, mapping = as.character(mapping))

x = readRDS(fnames[1])

all_effects = lapply(seq_along(fnames),
  function(i) {
    x = readRDS(fnames[i])
    x = x$summary
    x = data.frame(variable = rownames(x), x)
    x = dplyr::filter(x, grepl('gene_shift', variable))
    x = dplyr::mutate(x, mapping = sub('gene_shift\\[', '', variable))
    x = dplyr::mutate(x, mapping = sub('\\]', '', mapping))
    x = dplyr::left_join(x, gene_index_mapping, by = 'mapping')
    x = dplyr::mutate(x, chain = i)
    easy_effects = inner_join(x, data.frame(gene = easy), by = 'gene')
    moderate_effects = inner_join(x, data.frame(gene = moderate), by = 'gene')
    difficult_effects = inner_join(x, data.frame(gene = difficult), by = 'gene')
    list(easy = easy_effects,
      moderate = moderate_effects,
      difficult = difficult_effects)
  })

all_easy_effects = lapply(all_effects,
  function(x) {
    x$easy
  })
all_easy_effects = bind_rows(all_easy_effects)
all_easy_summary = group_by(all_easy_effects, gene)
all_easy_summary = summarize(all_easy_summary, mean = mean(Mean))

all_moderate_effects = lapply(all_effects,
  function(x) {
    x$moderate
  })
all_moderate_effects = bind_rows(all_moderate_effects)
all_moderate_summary = group_by(all_moderate_effects, gene)
all_moderate_summary = summarize(all_moderate_summary, mean = mean(Mean))

all_difficult_effects = lapply(all_effects,
  function(x) {
    x$difficult
  })
all_difficult_effects = bind_rows(all_difficult_effects)
all_difficult_summary = group_by(all_difficult_effects, gene)
all_difficult_summary = summarize(all_difficult_summary, mean = mean(Mean))

par(mfrow = c(1, 3))
hist(abs(all_easy_summary$mean), breaks = 10, xlim = c(0, 1))
hist(abs(all_moderate_summary$mean), breaks = 10, xlim = c(0, 1))
hist(abs(all_difficult_summary$mean), breaks = 10, xlim = c(0, 1))
dev.off()


save(all_easy_summary, all_moderate_summary, all_difficult_summary,
  file = 'coverage_effect_sizes.RData')

# adapted from `effect_size_figure.R`
effect_sizes = bind_rows(
  dplyr::mutate(all_easy_summary, type = 'easy'),
  dplyr::mutate(all_moderate_summary, type = 'moderate'),
  dplyr::mutate(all_difficult_summary, type = 'difficult')
)
effect_sizes = mutate(effect_sizes,
  type = factor(type, levels = c('easy', 'moderate', 'difficult')))

cbb = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

set.seed(42)
p = ggplot(effect_sizes, aes(type, abs(mean), color = type))
p = p + geom_boxplot(outlier.shape = NA, size = 1.5)
p = p + geom_jitter(height = 0, size = 2, alpha = 0.5, fill = 'black', shape = 21)
p = p + scale_color_manual(values = cbb)
p = p + ylab('absolute effect size')
p = p + xlab('effect size class')
p = p + theme(legend.position = 'none')
save_plot('img/effect_size_distribution.pdf', p, base_width = 7, base_height = 5)

table(sign(effect_sizes$mean)) / length(effect_sizes$mean)


es_summary = tapply(abs(effect_sizes$mean), effect_sizes$type, summary)
tapply(abs(effect_sizes$mean), effect_sizes$type, summary)

gamma_mle = function(data) {
  ab = optim(c(1, 1),
    function(x) {
      y = abs(data)
      ll = dgamma(y, shape = x[1], rate = x[2], log = TRUE)
      -sum(ll)
    }, method = 'Nelder-Mead')
  ab_hat = ab$par
  ab_hat
}

gamma_fits = tapply(abs(effect_sizes$mean), effect_sizes$type, gamma_mle)

x = seq(es_summary[['difficult']]['1st Qu.'], es_summary[['difficult']]['Max.'], length.out = 1000)

plot(x,
  dgamma(x, shape = gamma_fits[['difficult']][1], rate = gamma_fits[['difficult']][2]),
  col = 'white', xlim = c(min(effect_sizes$mean), max(effect_sizes$mean)))

for (i in 1:length(gamma_fits)) {
  g = gamma_fits[[i]]
  x = seq(es_summary[[i]]['1st Qu.'], es_summary[[i]]['Max.'], length.out = 1000)
  lines(x,  dgamma(x, shape = g[1], rate = g[2]), col = i)
}

gamma_all = gamma_mle(effect_sizes$mean)
lines(x, dgamma(x, shape = gamma_all[1], rate = gamma_all[2]), col = 'blue')

dev.off()

# look at the quantiles of these distributions and then make a choice about how to cut

summary(abs(effect_sizes$mean))
lower_cutoff = quantile(abs(effect_sizes$mean), probs = 0.05)

threshold_mean = effect_sizes$mean[abs(effect_sizes$mean) >= lower_cutoff]
threshold_mean = abs(threshold_mean)

dens = density(abs(threshold_mean))

n = 1000
tmp = sample(threshold_mean, n, replace = TRUE) + rnorm(n, 0, dens$bw)

cutoffs = seq(0, 1, length.out = 4)

cutoff_list = list()
for (i in 1:3) {
  cutoff_list[[i]] = quantile(threshold_mean, probs = c(cutoffs[i], cutoffs[i + 1]))
}

sample_from_density = function(n_samples, population, bandwidth, lower_cutoff, upper_cutoff) {
  samples = vector('numeric', n_samples)
  for (i in 1:n_samples) {
    current_sample = NULL
    repeat {
      current_sample = sample(population, 1, replace = TRUE) + rnorm(1, 0, bandwidth)
      if (current_sample >= lower_cutoff && current_sample <= upper_cutoff) {
        break
      }
    }
    samples[i] = current_sample
  }
  samples
}

n_samples_per_group = 10000

set.seed(102)
effect_size_samples = lapply(cutoff_list,
  function(cutoffs) {
    effect_sizes = Filter(function(x) x <= cutoffs[2] && cutoffs[1] <= x, threshold_mean)
    s = sample_from_density(n_samples_per_group, effect_sizes, dens$bw, cutoffs[1], cutoffs[2])
    effect_direction = sample(rep(c(-1, 1), length.out = length(s)))
    s * effect_direction
  })

# pdf('img/effect_size_ecdf.pdf', width = 5, height = 5, units = 'in', pointsize = 12, res = 1200)
pdf('img/effect_size_ecdf.pdf', width = 5, height = 5, pointsize = 12)
plot(ecdf(abs(effect_sizes$mean)), main = 'absolute value of effect sizes')
lines(dens$x, cumsum(dens$y) / sum(dens$y), col = cbb[1], lwd = 2)
lines(ecdf(abs(unlist(effect_size_samples))), col = cbb[2], lwd = 2)
legend(x = 'bottomright',
  legend = c('empirical', 'original smoothed', 'thresholded for simulations'),
  col = c('black', cbb[1:2]), lwd  = c(1, 2, 2))
dev.off()

sample_classes = c('difficult', 'moderate', 'easy')

smooth_effect_sizes = lapply(seq_along(sample_classes),
  function(i) {
    data.frame(mean = effect_size_samples[[i]],
      type = rep(sample_classes[i], length(effect_size_samples[[i]])))
  })
smooth_effect_sizes = bind_rows(smooth_effect_sizes)
smooth_effect_sizes = dplyr::mutate(smooth_effect_sizes, type = factor(type, levels = c('easy', 'moderate', 'difficult')))

p = ggplot(smooth_effect_sizes, aes(type, abs(mean), color = type))
p = p + geom_boxplot(outlier.shape = NA, size = 1.5)
p = p + scale_color_manual(values = cbb)
p = p + ylab('absolute effect size')
p = p + xlab('effect size class')
p = p + theme(legend.position = 'none')
save_plot('img/smooth_effect_size_distribution.pdf', p, base_width = 7, base_height = 5)
