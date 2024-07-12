# when debugging
cur_seed = 326
n_samples = 20000
n_burnin = 10000
thin = 1

# cur_seed = as.integer(argv[1])
# n_samples = as.integer(argv[2])
# n_burnin = as.integer(argv[3])
# thin = as.integer(argv[4])

library('stringr')
library('dplyr')
library('tidyr')
library('nimble')
library('ggplot2')
library('cowplot')
library('ggrepel')
theme_set(theme_cowplot())
library('waterbear')
library('survcomp')
library('xtable')


# install('~/Dropbox/postdoc/waterbear/')


DATA_HOME = '/oak/stanford/groups/pritch/users/jake/waterbear/high_moi/count_files'
raw_counts = read.table(paste0(DATA_HOME, '/high_low_MOI_screens_D1_D2_D3_2022-04-23.count.txt'), header = TRUE,
                        stringsAsFactors = FALSE, sep = '\t')


ordering = c('Q1', 'Q2', 'Q3', 'Q4')

c_name = grep('Low_Coverage_High_MOI', colnames(raw_counts), value = TRUE)
sample_mapping = data.frame(c_name = c_name)
sample_mapping = mutate(sample_mapping, sample = str_extract(c_name, 'D1|D2|D3'))
sample_mapping = mutate(sample_mapping, bin = str_extract(c_name, 'Q[1-4]'))

wb_array = wb_counts_to_array(raw_counts, sample_mapping, ordering)

gene_mapping = raw_counts[, c('sgRNA', 'Gene')]
colnames(gene_mapping) = c('guide', 'gene')

wo = wb_make_object(wb_array, gene_mapping, bin_size_prior = c(0.2, 0.3, 0.3, 0.2))

# debugonce(wb_make_object)
# wo_tmp = wb_make_object(array, gene_mapping, bin_size_prior = c(0.2, 0.3, 0.3, 0.2))


#debugonce(wb_em_start)
wo = wb_em_start(wo)
#wo = wo$wo

source('../waterbear/R/nimble.R')
n_model = nimbleModel(sample_specific_dispersion_model,
                      data = wo$data, constants = wo$const, inits = wo$init)

# nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
# nimbleOptions(MCMCsaveHistory = TRUE)

n_configuration = configureMCMC(n_model)
# sampler_control = list(order = wo$const$order)
# n_configuration$removeSamplers('gene_inclusion')
# n_configuration$addSampler(target = 'gene_inclusion', type = 'rank_RW',
#                            control = sampler_control)

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
n_chains = 4
seed = 42
system.time({
  samples = runMCMC(
    C_n_mcmc, niter = n_samples, nburnin = n_burnin,
    nchains = n_chains,
    setSeed = seed,
    summary = TRUE)
})

save(wo, samples, file = 'low_coverage_high_moi.RData')

saveRDS(samples, 'low_coverage_high_moi.rds')

saveRDS(samples, 'new_moi.rds')

samples = readRDS('new_moi.rds')

# XXX: these help with getting guide level posterior summaries
get_guide_summaries = function(wo, samples) {
  sample_summary = samples$summary$all.chains
  guide_to_guide_mapping = data.frame(guide_data_index = wo$const$guide_data_index,
    wb_guide_index = wo$const$guide_index)
  sample_summary = data.frame(variable = rownames(sample_summary), sample_summary)
  guides = dplyr::filter(sample_summary, grepl('guide_shift', variable))
  guides = dplyr::mutate(guides, wb_guide_index = as.integer(str_extract(variable, '(\\d)+')))
  guides = inner_join(guides, guide_to_guide_mapping, by = 'wb_guide_index')

  guides = inner_join(guides, wo$test_guide_names, by = c('guide_data_index' = 'i'))
  gene_inclusion = dplyr::filter(sample_summary, grepl('gene_inclusion', variable))
  gene_inclusion = dplyr::mutate(gene_inclusion,
    gene_mapping = as.integer(str_extract(variable, '(\\d)+')))
  gene_inclusion = inner_join(gene_inclusion,
    dplyr::distinct(dplyr::select(wo$test_guide_names, -c(i, guide))),
    by = c('gene_mapping' = 'mapping'))
  gene_shift = dplyr::filter(sample_summary, grepl('gene_shift', variable))
  gene_shift = dplyr::mutate(gene_shift,
    gene_mapping = as.integer(str_extract(variable, '(\\d)+')))
  gene_shift = inner_join(gene_shift,
    dplyr::distinct(dplyr::select(wo$test_guide_names, -c(i, guide))),
    by = c('gene_mapping' = 'mapping'))

  list(guides = guides, gene_inclusion = gene_inclusion, gene_shift = gene_shift)
}

# XXX: this is necessary for getting the local false sign rate
gene_posterior_mass = function(wo, samples, alpha) {
  stopifnot(is.list(samples$samples))
  # assumes chains are organizing the same exact manner
  # this should happen assuming it is the exact same model
  gene_columns = grep('gene_shift', colnames(samples$samples$chain1))
  n_chains = length(samples$samples)
  n_samples = nrow(samples$samples$chain1)
  tmp_values = vector('numeric', n_chains * n_samples)
  one_chain = vector('numeric', n_samples)
  n_genes = sum(gene_columns)

  gene_summary = lapply(gene_columns,
    function(gene) {
      tmp = sapply(samples$samples,
        function(chain) {
          chain[, gene]
          # }, one_chain)
        }, simplify = FALSE)
      gene_name = colnames(samples$samples$chain1)[gene]
      tmp = unlist(tmp)
      qs = quantile(tmp, probs = c(0 + alpha/2, 1 - alpha/2))
      # gt0 = mean(tmp > .Machine$double.eps)
      # lt0 = mean(tmp < .Machine$double.eps)
      mu = mean(tmp)
      gt0 = mean(tmp >= 0)
      lt0 = mean(tmp <= 0)
      data.frame(nimble_id = gene_name, lower = qs[1], upper = qs[2],
        gt0 = gt0, lt0 = lt0, mu = mu, row.names = NULL)
    })

  gene_summary = dplyr::bind_rows(gene_summary)
  gene_summary = dplyr::mutate(gene_summary, mapping = as.integer(stringr::str_extract(nimble_id, '(\\d)+')))
  gene_summary = inner_join(gene_summary,
    dplyr::distinct(dplyr::select(wo$test_guide_names, -c(i, guide))),
    by = c('mapping' = 'mapping'))
  gene_summary = dplyr::mutate(gene_summary, lfsr = pmin(gt0, lt0))
  dplyr::select(gene_summary, -nimble_id, -mapping)
}

samples$summary$all.chains[grepl('dispersion', rownames(samples$summary$all.chains)), ]
samples$summary$all.chains[grepl('psi', rownames(samples$summary$all.chains)), ]
samples$summary$all.chains[grepl('sigma', rownames(samples$summary$all.chains)), ]


# run this
level = 0.10

guide_summaries = get_guide_summaries(wo, samples)
gene_posterior_summary = gene_posterior_mass(wo, samples, level)

gene_columns = grep('gene_shift', colnames(samples$samples$chain1))

chain = samples$samples$chain1

tmp = sapply(samples$samples,
  function(chain) {
    chain[, gene_columns[1]]
  }, simplify = FALSE)

interesting_genes = data.frame(gene =
  c('IL2RA',
    'STAT5B',
    'MYB',
    'JAK3',
    'MED3',
    'CBFB'))

hchm_d1d2$gene_inclusion = hchm_d1d2$guide_summaries$gene_inclusion

gene_shift = guide_summaries$gene_shift
gene_shift = inner_join(gene_shift, gene_posterior_summary, by = 'gene')
gene_shift = inner_join(gene_shift,
  dplyr::select(gene_inclusion, gene, pip = Mean), by = 'gene')

# gene_shift = inner_join(guide_summaries$gene_shift, interesting_genes, by = 'gene')
guide_shift = inner_join(guide_summaries$guides, interesting_genes, by = 'gene')

gene_shift = mutate(gene_shift, high_pip = pip > 1 - level)

# get 'significant' genes
sig_genes_wb = dplyr::select(filter(gene_shift, high_pip & lfsr > 1 - level / 2), gene)
sig_genes_wb =  dplyr::mutate(sig_genes_wb, wb_sig = TRUE)

p = ggplot(gene_shift, aes(lfsr, Mean, color = pip))
p = p + geom_point(alpha = 0.4)
p = p + geom_vline(xintercept = c(level / 2, 1 - level / 2), linetype = 2)
p = p + xlab('proportion of gene_shift > 0')
p = p + ylab('gene shift')
save_plot(p, file = 'lfsr_mean.pdf')

# this is done for sanity check, you can ignore this
gene_inclusion = tmp$gene_inclusion
high_pip = filter(gene_inclusion, Mean > 0.90)

high_pip_genes = inner_join(
  select(high_pip, gene_mapping),
  tmp$gene_shift,
  by = 'gene_mapping')
high_pip_genes = arrange(high_pip_genes, abs(Mean))

high_pip_and_not_zero = filter(high_pip_genes, !(X95.CI_low < 0 & 0 < X95.CI_upp))

high_pip_include_zero = filter(high_pip_genes, X95.CI_low < 0 & 0 < X95.CI_upp)

# this is to plot guides -- you can ignore this for now
guides_to_plot = inner_join(
  select(high_pip_include_zero[1:10, ], gene),
  tmp$guides)

p = ggplot(guides_to_plot, aes(gene, Mean, color = gene))
p = p + geom_pointrange(
  aes(ymin = X95.CI_low, ymax = X95.CI_upp),
  position=position_jitter(width=0.25))
p = p + geom_hline(yintercept = 0, linetype = 'dotted')
p = p + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
  legend.position = 'none')
save_plot(p, file = 'small_effect_genes.pdf')

p = ggplot(
  inner_join(
    select(high_pip_and_not_zero[1:10,], gene),
    tmp$guides),
  aes(gene, Mean, color = gene))
p = p + geom_pointrange(
  aes(ymin = X95.CI_low, ymax = X95.CI_upp),
  position=position_jitter(width=0.25))
p = p + geom_hline(yintercept = 0, linetype = 'dotted')
p = p + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
  legend.position = 'none')
save_plot(p, file = 'large_effect_genes.pdf')

genes = waterbear:::wb_recode(samples, wo)
gene_shift = waterbear:::wb_recode(samples, wo, 'gene_shift')

guide_level = inner_join(tmp$guides, top_genes, by = 'gene')


# END: this is to plot guides -- you can ignore this for now

# load mageck results
library('readr')
high_moi_results <- read_tsv('../data/High_MOI_CD25_Q1_Q4_pos_enrichment_2022-03-02.gene_summary.txt')
high_moi_results = dplyr::select(high_moi_results,
  id, neg_fdr = `neg|fdr`, pos_fdr = `pos|fdr`, lfc = `pos|lfc`)

gene_shift = left_join(gene_shift, sig_genes_wb, by = 'gene')
gene_shift = mutate(gene_shift, wb_sig = ifelse(is.na(wb_sig), FALSE, TRUE))

compare = inner_join(
  high_moi_results,
  gene_shift,
  by = c('id' = 'gene'))
compare = mutate(compare, sig = case_when(
  wb_sig & (pos_fdr < 0.1 | neg_fdr < 0.1) ~ 'Both',
  wb_sig ~ 'waterbear',
  (pos_fdr < 0.1 | neg_fdr < 0.1 ) ~ 'mageck',
  TRUE ~ 'neither'
  ))
compare = mutate(compare,
  gene2 = case_when(sig == 'Both' ~ id,
                           TRUE ~ ''))
compare <- inner_join(compare, compare %>% group_by(sig) %>% summarise(n = n())) %>%
  mutate(label = paste0(sig, ' (', n, ')'))

library('cowplot')
library('colorspace')


p = ggplot(compare, aes(x = -lfc, y = Mean, color = label, fill = label)) +
  geom_hline(yintercept = 0, linetype = 2, color = '#bdbdbd') +
  geom_vline(xintercept = 0, linetype = 2, color = '#bdbdbd') +
  xlab('MAGeCK guide RNA enrichment\nhigh/low FACS bins (log2 fold change)') +
  ylab('Water bear gene shift') +
  geom_point(size = 1, pch = 21, alpha = 0.4) +
  # geom_text_repel(aes(label = label),
  #                 seed = 123,
  #                 show.legend = F,
  #                 segment.color = '#bdbdbd',
  #                 color = 'black',
  #                 segment.size	= 0.25,
  #                 force = 5,
  #                 size = 1.75,
  #                 min.segment.length = 0.1) +
  # geom_point(data = filter(compare, sig == 'Both'),
  #            fill = '#2166ac', color = darken('#2166ac', 0.3), size = 1, pch = 21) +
  # geom_point(data = filter(compare, sig == 'waterbear'),
  #            fill = '#D55E00', color = darken('#D55E00', 0.3), size = 1, pch = 21) +
  # scale_fill_manual(name = "Significant:", values = c('#2166ac', '#999999', '#D55E00')) +
  # scale_color_manual(name = "Significant:", values = darken(c('#2166ac', '#999999', '#D55E00'), .3)) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.05, 1),
        legend.justification = c("left", "top"),
        legend.spacing.x = unit(0, 'cm'),
        legend.background = element_blank())
save_plot(p, file = 'lfsr_wb_mageck_moi5.pdf')


interesting_genes = dplyr::select(dplyr::filter(compare, sig == 'waterbear'), gene = id)

guide_shift = inner_join(guide_summaries$guides, interesting_genes, by = 'gene')

# TODO: look at mageck effect size
# TODO: pool

# end mageck

p = ggplot(guide_shift, aes(guide, Mean, color = gene))
p = p + geom_pointrange(
  aes(ymin = X95.CI_low, ymax = X95.CI_upp),
  position=position_jitter(width=0.25))
p = p + geom_hline(yintercept = 0, linetype = 'dotted')
p = p + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
  legend.position = 'none')
save_plot(p, file = 'small_effect_genes.pdf')


test <- samples$summary$all.chains
test2 <- tibble(gene = rownames(test), mean = test[, 'Mean'])
gene_inclusion <- filter(test2, grepl('gene_inclusion', gene)) %>%
  mutate(id = str_extract(gene, '(\\d)+')) %>%
  dplyr::select(id, gene_inclusion = mean)
gene_inclusion
gene_shift <- filter(test2, grepl('gene_shift', gene)) %>%
  mutate(id = str_extract(gene, '(\\d)+')) %>%
  dplyr::select(id, gene_shift = mean)
results <- inner_join(gene_inclusion, gene_shift, by = 'id')

a <- wo$test_guide_names %>%
  dplyr::select(gene, mapping) %>%
  distinct() %>%
  mutate(mapping = as.character(mapping))

results <- inner_join(results, a, by = c('id' = 'mapping'))
write_csv(x = results, file = '~/Desktop/waterbear_high_moi_results.csv')



# visualize raw data
# add integers to gene names
# paragraphs describing waterbear


c_name = grep('High_Coverage_Low_MOI', colnames(raw_counts), value = TRUE)
sample_mapping = data.frame(c_name = c_name)
sample_mapping = mutate(sample_mapping, sample = str_extract(c_name, 'D1|D2|D3'))
sample_mapping = mutate(sample_mapping, bin = str_extract(c_name, 'Q[1-4]'))
sample_mapping = sample_mapping[complete.cases(sample_mapping), ]

wb_array = wb_counts_to_array(raw_counts, sample_mapping, ordering)

gene_mapping = raw_counts[, c('sgRNA', 'Gene')]
colnames(gene_mapping) = c('guide', 'gene')

wo = wb_make_object(wb_array, gene_mapping, bin_size_prior = c(0.2, 0.3, 0.3, 0.2))

# debugonce(wb_make_object)
# wo_tmp = wb_make_object(array, gene_mapping, bin_size_prior = c(0.2, 0.3, 0.3, 0.2))


#debugonce(wb_em_start)
wo = wb_em_start(wo)

source('../waterbear/R/nimble.R')
n_model = nimbleModel(sample_specific_dispersion_model,
                      data = wo$data, constants = wo$const, inits = wo$init)

n_configuration = configureMCMC(n_model)

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
n_chains = 4
seed = 42
system.time({
  samples = runMCMC(
    C_n_mcmc, niter = n_samples, nburnin = n_burnin,
    nchains = n_chains,
    setSeed = seed,
    summary = TRUE)
})
save(wo, samples, file = 'high_coverage_low_moi.RData')

c_name = grep('Low_Coverage_Low_MOI', colnames(raw_counts), value = TRUE)
sample_mapping = data.frame(c_name = c_name)
sample_mapping = mutate(sample_mapping, sample = str_extract(c_name, 'D1|D2|D3'))
sample_mapping = mutate(sample_mapping, bin = str_extract(c_name, 'Q[1-4]'))
sample_mapping = sample_mapping[complete.cases(sample_mapping), ]
wb_array = wb_counts_to_array(raw_counts, sample_mapping, ordering)

gene_mapping = raw_counts[, c('sgRNA', 'Gene')]
colnames(gene_mapping) = c('guide', 'gene')

wo = wb_make_object(wb_array, gene_mapping, bin_size_prior = c(0.2, 0.3, 0.3, 0.2))
wo = wb_em_start(wo)

source('../waterbear/R/nimble.R')
n_model = nimbleModel(sample_specific_dispersion_model,
                      data = wo$data, constants = wo$const, inits = wo$init)

n_configuration = configureMCMC(n_model)

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
n_chains = 4
seed = 42
system.time({
  samples = runMCMC(
    C_n_mcmc, niter = n_samples, nburnin = n_burnin,
    nchains = n_chains,
    setSeed = seed,
    summary = TRUE)
})
save(wo, samples, file = 'low_coverage_low_moi.RData')

c_name = grep('High_Coverage_High_MOI', colnames(raw_counts), value = TRUE)
sample_mapping = data.frame(c_name = c_name)
sample_mapping = mutate(sample_mapping, sample = str_extract(c_name, 'D1|D2|D3'))
sample_mapping = mutate(sample_mapping, bin = str_extract(c_name, 'Q[1-4]'))
sample_mapping = sample_mapping[complete.cases(sample_mapping), ]
wb_array = wb_counts_to_array(raw_counts, sample_mapping, ordering)

gene_mapping = raw_counts[, c('sgRNA', 'Gene')]
colnames(gene_mapping) = c('guide', 'gene')

wo = wb_make_object(wb_array, gene_mapping, bin_size_prior = c(0.2, 0.3, 0.3, 0.2))
wo = wb_em_start(wo)

source('../waterbear/R/nimble.R')
n_model = nimbleModel(sample_specific_dispersion_model,
                      data = wo$data, constants = wo$const, inits = wo$init)

n_configuration = configureMCMC(n_model)

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
n_chains = 4
seed = 42
system.time({
  samples = runMCMC(
    C_n_mcmc, niter = n_samples, nburnin = n_burnin,
    nchains = n_chains,
    setSeed = seed,
    summary = TRUE)
})
save(wo, samples, file = 'high_coverage_high_moi.RData')

# analysis of sampling

low_high = new.env()
load(file = 'low_coverage_high_moi.RData', envir = low_high)

high_low = new.env()
load(file = 'high_coverage_low_moi.RData', envir = high_low)

low_low = new.env()
load(file = 'low_coverage_low_moi.RData', envir = low_low)

high_high = new.env()
load(file = 'high_coverage_high_moi.RData', envir = high_high)

# run this
level = 0.10

# paper_labels = data.frame(gene = c('MYC', 'JAK3', 'RELA', 'ATXN7L3', 'PTEN', 'TP53'))
paper_labels = data.frame(gene = c('MYC', 'JAK3', 'RELA', 'ATXN7L3', 'PTEN'))
paper_labels = mutate(paper_labels, gene_paper = gene)

lh_guide_summaries = get_guide_summaries(low_high$wo, low_high$samples)
lh_gene_posterior_summary = gene_posterior_mass(low_high$wo, low_high$samples, level)

hl_guide_summaries = get_guide_summaries(high_low$wo, high_low$samples)
hl_gene_posterior_summary = gene_posterior_mass(high_low$wo, high_low$samples, level)


ll_gene_posterior_summary = gene_posterior_mass(low_low$wo, low_low$samples, level)

hh_gene_posterior_summary = gene_posterior_mass(high_high$wo, high_high$samples, level)

saveRDS(hh_gene_posterior_summary, 'high_coverage_high_moi.rds')

# gene_columns = grep('gene_shift', colnames(samples$samples$chain1))


lh_hl_summary = inner_join(
  lh_gene_posterior_summary,
  hl_gene_posterior_summary,
  # dplyr::select(lh_gene_posterior_summary, mu, gene, lfsr),
  # dplyr::select(hl_gene_posterior_summary, mu, gene, lfsr),
  suffix = c('_lh', '_hl'),
  by = c('gene')
)

lh_hl_summary = mutate(
  lh_hl_summary,
  sig = case_when(
    lfsr_lh < level & lfsr_hl < level ~ 'Both',
    lfsr_lh < level & lfsr_hl >= level ~ 'Low coverage, high MOI',
    lfsr_lh >= level & lfsr_hl < level ~ 'High coverage, low MOI',
    TRUE ~ 'Neither'
  )
)

interesting_genes = data.frame(gene =
  c('IL2RA',
    'STAT5B',
    'MYB',
    'JAK3',
    'MED30',
    'CBFB',
    'ZNF626',
    'DUX4',
    'MED11'))
interesting_genes = mutate(interesting_genes, gene_text = gene)
lh_hl_summary = left_join(lh_hl_summary, interesting_genes, by = 'gene')
lh_hl_summary = left_join(lh_hl_summary, paper_labels, by = 'gene')

lh_hl_counts = lh_hl_summary %>% group_by(sig) %>% summarize(n = length(gene))
lh_hl_counts = mutate(lh_hl_counts, label_text = sort(paste0(sig, ' (', n, ')')),
  color = c('#009E73', '#0072B2', '#E69F00', '#999999'))
names(lh_hl_counts$color) = lh_hl_counts$sig
names(lh_hl_counts$label_text) = lh_hl_counts$sig




all_gene_summary = bind_rows(
  mutate(lh_gene_posterior_summary, method = 'low_high'),
  mutate(ll_gene_posterior_summary, method = 'low_low'),
  mutate(hl_gene_posterior_summary, method = 'high_low'),
  mutate(hh_gene_posterior_summary, method = 'high_high')
)

group_by(all_gene_summary, method) %>%
  summarize(sum(lfsr < level))

effect_size = pivot_wider(all_gene_summary, gene, names_from = method, values_from = mu)

cor(as.matrix(effect_size[, -1]), method = 'spearman')


data = as.matrix(effect_size[, -1])

library(PerformanceAnalytics)

chart.Correlation(data, histogram = TRUE, method = "spearman")

p = ggplot(all_gene_summary, aes(lfsr))
p = p + geom_histogram()
p = p + facet_wrap(~method)
p

gene_filter = 'DUX4'
p = ggplot(dplyr::filter(all_gene_summary, gene == gene_filter), aes(factor(method), mu))
p = p + geom_point()
p = p + geom_pointrange(aes(ymin = lower, ymax = upper))
p


all_lfsr = pivot_wider(all_gene_summary, gene, names_from = method, values_from = lfsr)
all_significant = all_lfsr[, -1] < 0.1
mode(all_significant) = 'integer'
all_significant = data.frame(all_lfsr$gene, all_significant)



library('UpSetR')
dir.create('img')

# upset(all_significant, sets = c('low_low', 'low_high', 'high_low', 'high_high'))
pdf('img/upset_moi.pdf', width = 7, height = 4)
upset(all_significant, order.by = 'freq')
dev.off()

lh_coverage = t(apply(low_high$wo$data$x, c(1, 2), sum))
ll_coverage = t(apply(low_low$wo$data$x, c(1, 2), sum))
hl_coverage = t(apply(high_low$wo$data$x, c(1, 2), sum))
hh_coverage = t(apply(high_high$wo$data$x, c(1, 2), sum))

all_coverage = lapply(
  list('low_high', 'low_low', 'high_low', 'high_high'),
  function(x) {
    y = get(x)
    coverage = t(apply(y$wo$data$x, c(1, 2), sum))
    current_coverage = pivot_longer(
      data.frame(guide = rownames(coverage), coverage),
      cols = -c(guide))
    mutate(current_coverage, experiment = x)
  })
all_coverage = bind_rows(all_coverage)

group_by(all_coverage, name, experiment) %>%
  summarize(med = median(value), q05 = )

p = ggplot(all_coverage, aes(name, value))
p = p + geom_boxplot(aes(color = experiment), outlier.shape = NA, coef = 0)
p = p + ylim(6000, 10000)
save_plot(p, file = 'img/guide_coverage.pdf')

pivot_longer(data.frame(guide = rownames(lh_coverage), lh_coverage), cols = -c(guide))

(lapply(list(lh_coverage, ll_coverage, hl_coverage, hh_coverage),))
  function(x) {


  })

summary(hl_coverage)
summary(hh_coverage)

# Supplementary figures of validation
# MAGeCK
mageck_lh = read.table('~/Dropbox/Water_bear_paper/low_high_moi_seq_data/mageck_results/Low_Coverage_High_MOI_CD25_Q1_Q4_pos_enrichment_2022-04-23.gene_summary.txt',
  header = TRUE)
mageck_lh = rename(mageck_lh, gene = id, effect_size = pos.lfc)
mageck_lh = mutate(mageck_lh, fdr = pmin(neg.fdr, pos.fdr))
mageck_lh = select(mageck_lh, gene, mageck_effect_size = effect_size, mageck_fdr = fdr)

lh_waterbear_mageck = inner_join(lh_gene_posterior_summary, mageck_lh, by = 'gene')
lh_waterbear_mageck = mutate(lh_waterbear_mageck,
  sig = case_when(
    lfsr < level & mageck_fdr < level ~ 'Both',
    lfsr < level ~ 'waterbear',
    mageck_fdr < level ~ 'MAGeCK',
    TRUE ~ 'neither'
    ))
# lh_waterbear_mageck = left_join(lh_waterbear_mageck, interesting_genes, by = 'gene')
lh_waterbear_mageck = left_join(lh_waterbear_mageck, mutate(validated_targets, gene_text = gene), by = 'gene')

wm_color_mapping = c('Both' = '#E69F00',
  'waterbear' = '#56B4E9',
  'MAGeCK' = '#009E73',
  'neither' = '#999999')

p = ggplot(lh_waterbear_mageck, aes(-mageck_effect_size, mu))
p = p + geom_hline(yintercept = 0, linetype = 2, color = '#bdbdbd')
p = p + geom_vline(xintercept = 0, linetype = 2, color = '#bdbdbd')
p = p + geom_point(aes(color = sig), alpha = 0.3)
p = p + geom_text_repel(aes(label = gene_text, color = sig), min.segment.length = 0,
  box.padding = 1, max.overlaps = Inf)
p = p + theme(legend.position = c(0.8, 0.1))
p = p + labs(color = 'Significance')
p = p + xlab('MAGeCK effect size')
p = p + ylab('waterbear effect size')
p = p + scale_color_manual(values = wm_color_mapping)
p
save_plot(p, file = 'img/scatter_lh_waterbear_mageck.pdf', base_height = 7, base_asp = 1)

# validated_targets = read.csv('../data/cd25_validated_targets.csv', header = TRUE,
validated_targets = read.csv('../data/cd25_validated_targets\ revised.csv', header = TRUE,
  stringsAsFactors = FALSE)

arrayed_flow = read.csv('../data/waterbear_IL2RA_screen_validation_values.csv', header = TRUE,
  stringsAsFactors = FALSE)
arrayed_flow = mutate(arrayed_flow, gene = sub(' KO', '', Gene))


validated_targets = dplyr::rename(validated_targets, gene = Gene_knocked_out)
validated_targets = filter(validated_targets, gene != 'Non-Targeting')
validated_targets = select(validated_targets, -Screen_target)
validated_targets = mutate(validated_targets,
  validation_direction = case_when(
    grepl('Decrease', Validation_status) ~ 'Decrease',
    grepl('Increase', Validation_status) ~ 'Increase',
    grepl('not-validated', Validation_status) ~ 'Not validated',
    TRUE ~ 'null'
    ))

validated_targets = mutate(validated_targets, validation_direction = factor(validation_direction,
    levels = c('Decrease', 'null', 'Increase', 'Not validated')))
validated_targets = mutate(validated_targets,
  validation_sign = ifelse('Increase' == validation_direction, 1, -1)
    )
validated_targets = mutate(validated_targets,
  validation_sign = ifelse('null' == validation_direction, 0, validation_sign))
validated_targets = mutate(validated_targets,
  validation_sign = ifelse('Not validated' == validation_direction, NA, validation_sign))
validated_targets = dplyr::filter(validated_targets,
  !grepl('Not validated', validation_direction))
validated_targets = left_join(validated_targets, select(arrayed_flow, gene, mfi_avg))
validated_targets = mutate(validated_targets, mfi_rank = rank(mfi_avg))
validated_targets = mutate(validated_targets, mfi_rank_centered = mfi_rank - 14.5)

write.table(validated_targets, file = 'df_validated_targets_revised.tsv',
  row.names = FALSE, sep = '\t')

print(xtable(select(validated_targets, -Validation_status, -validation_sign, mfi_avg), type = 'latex'),
  file = 'validated_targets.tex')

lh_validated = inner_join(lh_gene_posterior_summary, validated_targets, by = 'gene')
lh_validated = arrange(lh_validated, mu)
lh_validated = mutate(lh_validated, posterior_rank = 1:length(mu))

color_mapping = c(
  'true positive' = "#56B4E9",
  'false positive' = "#D55E00",
  'true negative' = "#0072B2",
  'false negative' = "#CC79A7"
)
shape_mapping = c(
  'true positive' = 21,
  'false positive' = 22,
  'true negative' = 24,
  'false negative' = 23
)

lh_validated = mutate(lh_validated,
  status = case_when(
    sign(mu) == validation_sign & lfsr < 0.10 ~ 'true positive',
    lfsr > 0.10 & validation_sign == 0 ~ 'true negative',
    lfsr < 0.10 & validation_sign == 0 ~ 'false positive',
    lfsr > 0.10 & validation_sign != 0 ~ 'false negative',
  ))

# p = ggplot(lh_validated, aes(mfi_avg, mu))
p = ggplot(lh_validated, aes((mfi_rank_centered), mu))
p = p + geom_hline(yintercept = 0, linetype = 2, color = '#bdbdbd')
p = p + geom_vline(xintercept = 0, linetype = 2, color = '#bdbdbd')
p = p + geom_pointrange(aes(ymin = lower, ymax = upper,
    color = status, shape = status)
  # position = position_dodge2(width = 0.75)
  )
p = p + scale_x_discrete()
p = p + xlab('Arrayed knockout direction')
p = p + ylab('waterbear gene posterior')
p = p + scale_color_manual(values = color_mapping)
p = p + scale_shape_manual(values = shape_mapping)
p = p + geom_text_repel(aes(label = gene, color = status), min.segment.length = 0,
  max.overlaps = Inf, box.padding = 0.75)
p = p + theme(legend.position = c(0.65, 0.1), legend.title = element_blank())
p
save_plot(p, file = 'img/lh_waterbear.pdf', base_height = 7, base_asp = 1)
# save_plot(p, file = 'img/lh_waterbear2.pdf', base_height = 7, base_asp = 1)

lh_validated = left_join(lh_validated, paper_labels)
lh_validated = mutate(lh_validated, gene_paper = ifelse(is.na(gene_paper),
    '', gene_paper))

theme_set(theme_cowplot(20))
# p = ggplot(lh_validated, aes(validation_direction, mu))

p = ggplot(lh_validated, aes((mfi_rank_centered), mu))
# p = ggplot(lh_validated, aes((mfi_avg), mu))
p = p + geom_hline(yintercept = 0, linetype = 2, color = '#bdbdbd')
p = p + geom_vline(xintercept = 0, linetype = 2, color = '#bdbdbd')
p = p + geom_pointrange(aes(ymin = lower, ymax = upper,
    color = status, shape = status)
  # position = position_dodge2(width = 0.75)
  )
p = p + scale_x_discrete()
p = p + xlab('Arrayed knockout rank')
p = p + ylab('Waterbear gene posterior')
p = p + scale_color_manual(values = color_mapping)
p = p + scale_shape_manual(values = shape_mapping)
p = p + geom_text_repel(aes(label = gene_paper, color = status), min.segment.length = 0,
# p = p + geom_text_repel(aes(label = gene_paper), min.segment.length = 0,
  max.overlaps = Inf, box.padding = 3
  # position = position_dodge2(width = 0.75)
)
p = p + theme(legend.position = c(0.70, 0.1), legend.title = element_blank())
p
save_plot(p, file = 'img/lh_waterbear_paper.pdf', base_height = 7, base_asp = 1)

lh_mageck_validated = inner_join(validated_targets, mageck_lh, by = 'gene')
lh_mageck_validated = mutate(lh_mageck_validated, mageck_effect_size = -mageck_effect_size)
lh_mageck_validated = mutate(lh_mageck_validated,
  status = case_when(
    sign(mageck_effect_size) == validation_sign & mageck_fdr < 0.10 ~ 'true positive',
    # sign(mageck_effect_size) != validation_sign & mageck_fdr < 0.10 & validation_sign != 0~ 'true positive (discordant)',
    mageck_fdr > 0.10 & validation_sign == 0 ~ 'true negative',
    mageck_fdr < 0.10 & validation_sign == 0 ~ 'false positive',
    mageck_fdr > 0.10 & validation_sign != 0 ~ 'false negative',
  ))
lh_mageck_validated = arrange(lh_mageck_validated, mageck_effect_size)

p = ggplot(lh_mageck_validated, aes((mfi_rank_centered), mageck_effect_size))
# p = ggplot(lh_mageck_validated, aes(validation_direction, mageck_effect_size))
p = p + geom_hline(yintercept = 0, linetype = 2, color = '#bdbdbd')
p = p + geom_vline(xintercept = 0, linetype = 2, color = '#bdbdbd')
p = p + geom_point(aes(color = status, shape = status),
  position = position_dodge2(width = 1),
  size = 4,
  stroke = 1
  )
p = p + scale_x_discrete()
p = p + xlab('Arrayed knockout direction')
p = p + ylab('MAGeCK gene log fold change')
p = p + scale_color_manual(values = color_mapping)
p = p + scale_shape_manual(values = shape_mapping)
p = p + geom_text_repel(aes(label = gene, color = status), min.segment.length = 0,
  max.overlaps = Inf, box.padding = 0.75,
  position = position_dodge2(width = 1)
)
p = p + theme(legend.position = c(0.70, 0.1), legend.title = element_blank())
p

# save_plot(p, file = 'img/lh_mageck.pdf', base_height = 7, base_asp = 1)
save_plot(p, file = 'img/lh_mageck_revised.pdf', base_height = 7, base_asp = 1)

# # now comparing solely the coverage
# lh_hl_colors = c(
#   'Both'  = '#009E73',
#   'High coverage, low MOI' = '#0072B2',
#   'Low coverage, high MOI' = '#E69F00',
#   'Neither' = '#999999'
# )

p = ggplot(lh_hl_summary, aes(mu_lh, mu_hl))
p = p + geom_hline(yintercept = 0, linetype = 2, color = '#bdbdbd', alpha = 0.75)
p = p + geom_vline(xintercept = 0, linetype = 2, color = '#bdbdbd', alpha = 0.75)
p = p + geom_abline(alpha = 0.05)
p = p + geom_point(aes(fill = sig, alpha = sig, size = sig), color = 'black', shape = 21)
p = p + geom_text_repel(aes(label = gene_text, color = sig), min.segment.length = 0)
p = p + theme(legend.position = c(0.6, 0.2))
p = p + xlab('effect size (low coverage, high MOI)')
p = p + ylab('effect size (high coverage, low MOI)')
p = p + scale_fill_manual(values = lh_hl_counts$color, labels = lh_hl_counts$label_text)
p = p + scale_color_manual(values = lh_hl_counts$color, labels = lh_hl_counts$label_text)
p = p + scale_size_manual(values = c(rep(2, 3), 0.4), labels = lh_hl_counts$label_text)
p = p + scale_alpha_manual(values = c(rep(0.5, 3), 0.1), labels = lh_hl_counts$label_text)
p
save_plot(p, file = 'img/scatter_lh_hl.pdf', base_height = 7, base_asp = 1)

theme_set(theme_cowplot(20))
p = ggplot(lh_hl_summary, aes(mu_lh, mu_hl))
p = p + geom_hline(yintercept = 0, linetype = 2, color = '#bdbdbd', alpha = 0.75)
p = p + geom_vline(xintercept = 0, linetype = 2, color = '#bdbdbd', alpha = 0.75)
p = p + geom_abline(alpha = 0.05)
p = p + geom_point(aes(fill = sig, alpha = sig, size = sig), color = 'black', shape = 21)
p = p + geom_text_repel(aes(label = gene_paper, color = sig), min.segment.length = 0,
  box.padding = 2.5, max.overlaps = Inf, size = 5)
# p = p + theme(legend.position = c(0.5, 0.2))
p = p + theme(legend.position = c(0.05, 0.92))

p = p + xlab('Effect size (low coverage, high MOI)')
p = p + ylab('Effect size (high coverage, low MOI)')
p = p + scale_fill_manual(values = lh_hl_counts$color, labels = lh_hl_counts$label_text, name = '')
p = p + scale_color_manual(values = lh_hl_counts$color, labels = lh_hl_counts$label_text, name = '')
p = p + scale_size_manual(values = c(rep(2, 3), 0.4), labels = lh_hl_counts$label_text, name = '')
p = p + scale_alpha_manual(values = c(rep(0.5, 3), 0.1), labels = lh_hl_counts$label_text, name = '')
p
save_plot(p, file = 'img/scatter_lh_hl_paper.pdf', base_height = 7, base_asp = 1)


## grant fig
##
theme_set(theme_cowplot(40))
p = ggplot(lh_hl_summary, aes(mu_lh, mu_hl))
p = p + geom_hline(yintercept = 0, linetype = 2, color = '#bdbdbd', alpha = 0.75)
p = p + geom_vline(xintercept = 0, linetype = 2, color = '#bdbdbd', alpha = 0.75)
p = p + geom_abline(alpha = 0.05)
p = p + geom_point(aes(fill = sig, alpha = sig, size = sig),
  color = 'black', shape = 21, size = 5)
# p = p + geom_text_repel(aes(label = gene_paper, color = sig), min.segment.length = 0,
  # box.padding = 2.5, max.overlaps = Inf, size = 5)
# p = p + theme(legend.position = c(0.5, 0.2))
p = p + theme(legend.position = 'none')
p = p + xlab('Effect size (High MOI)')
p = p + ylab('Effect size (Low MOI)')
p = p + scale_fill_manual(values = lh_hl_counts$color, labels = lh_hl_counts$label_text, name = '')
p = p + scale_color_manual(values = lh_hl_counts$color, labels = lh_hl_counts$label_text, name = '')
p = p + scale_size_manual(values = c(rep(2, 3), 0.4), labels = lh_hl_counts$label_text, name = '')
p = p + scale_alpha_manual(values = c(rep(0.5, 3), 0.1), labels = lh_hl_counts$label_text, name = '')
p
save_plot(p, file = 'img/scatter_lh_hl_grant.pdf', base_height = 7, base_asp = 1.2)
## end grant fig

lh_hl_disagreements = filter(lh_hl_summary, sig != 'Both' & sig != 'Neither')
lh_hl_disagreements = filter(lh_hl_disagreements, sign(mu_lh) != sign(mu_hl))

lh_hl_summary %>% group_by(sig) %>% summarize(total = length(gene))
lh_hl_disagreements %>% group_by(sig) %>% summarize(total = length(gene))

p = ggplot(lh_hl_disagreements, aes(mu_lh, mu_hl))
p = p + geom_hline(yintercept = 0, linetype = 2, color = '#bdbdbd', alpha = 0.75)
p = p + geom_vline(xintercept = 0, linetype = 2, color = '#bdbdbd', alpha = 0.75)
p = p + geom_abline(alpha = 0.05)
p = p + geom_pointrange(aes(fill = sig, xmin = lower_lh, xmax = upper_lh),
  color = 'black', shape = 21, alpha = 0.2)
p = p + geom_pointrange(aes(fill = sig, ymin = lower_hl, ymax = upper_hl), color = 'black', shape = 21, alpha = 0.2)
p = p + theme(legend.position = c(0.6, 0.2))
p = p + xlab('effect size (low coverage, high MOI)')
p = p + ylab('effect size (high coverage, low MOI)')
p = p + coord_cartesian(ylim = c(-0.5, 0.5))
p

# p = p + geom_text_repel(aes(label = gene_text, color = sig), min.segment.length = 0)
# p = p + scale_fill_manual(values = lh_hl_counts$color, labels = lh_hl_counts$label_text)
# p = p + scale_color_manual(values = lh_hl_counts$color, labels = lh_hl_counts$label_text)
# p = p + scale_size_manual(values = c(rep(2, 3), 0.4), labels = lh_hl_counts$label_text)
# p = p + scale_alpha_manual(values = c(rep(0.5, 3), 0.1), labels = lh_hl_counts$label_text)


mageck_hl = read.table('~/Dropbox/Water_bear_paper/low_high_moi_seq_data/mageck_results/High_Coverage_Low_MOI_CD25_Q1_Q4_pos_enrichment_2022-04-23.gene_summary.txt',
  header = TRUE)
mageck_hl = rename(mageck_hl, gene = id, effect_size = pos.lfc)
mageck_hl = mutate(mageck_hl, fdr = pmin(neg.fdr, pos.fdr))
mageck_hl = select(mageck_hl, gene, mageck_effect_size = effect_size, mageck_fdr = fdr)

mageck_hl_lh = inner_join(
  mageck_hl,
  mageck_lh,
  by = 'gene',
  suffix = c('_hl', '_lh')
)

mageck_hl_lh = mutate(
  mageck_hl_lh,
  sig = case_when(
    mageck_fdr_lh < level & mageck_fdr_hl < level ~ 'Both',
    mageck_fdr_lh < level & mageck_fdr_hl >= level ~ 'Low coverage, high MOI',
    mageck_fdr_lh >= level & mageck_fdr_hl < level ~ 'High coverage, low MOI',
    TRUE ~ 'neither'
  )
)
mageck_hl_lh = left_join(mageck_hl_lh, interesting_genes, by = 'gene')

p = ggplot(mageck_hl_lh, aes(-mageck_effect_size_lh, -mageck_effect_size_hl))
p = p + geom_hline(yintercept = 0, linetype = 2, color = '#bdbdbd')
p = p + geom_vline(xintercept = 0, linetype = 2, color = '#bdbdbd')
p = p + geom_abline(alpha = 0.1)
p = p + geom_point(aes(color = sig), alpha = 0.3)
p = p + geom_text_repel(aes(label = gene_text, color = sig), min.segment.length = 0)
p = p + theme(legend.position = c(0.6, 0.2))
p = p + xlab('effect size (low coverage, high MOI)')
p = p + ylab('effect size (high coverage, low MOI)')
p = p + scale_color_manual(values = lh_hl_colors)
p
save_plot(p, file = 'img/mageck_scatter_lh_hl.pdf', base_height = 7, base_asp = 1)

hl_hh_summary = inner_join(
  dplyr::select(hl_gene_posterior_summary, mu, gene, lfsr),
  dplyr::select(hh_gene_posterior_summary, mu, gene, lfsr),
  suffix = c('_hl', '_hh'),
  by = c('gene')
)
hl_hh_summary = mutate(
  hl_hh_summary,
  sig = case_when(
    lfsr_hl < level & lfsr_hh < level ~ 'Both',
    lfsr_hl < level & lfsr_hh >= level ~ 'High coverage, low MOI',
    lfsr_hl >= level & lfsr_hh < level ~ 'High coverage, high MOI',
    TRUE ~ 'neither'
  )
)

hl_hh_summary = left_join(hl_hh_summary, interesting_genes, by = 'gene')

hl_hh_colors = c(
  'Both'  = '#56B4E9',
  'High coverage, low MOI' = '#009E73',
  'High coverage, high MOI' = '#D55E00',
  'neither' = '#999999'
)

p = ggplot(hl_hh_summary, aes(mu_hl, mu_hh))
p = p + geom_hline(yintercept = 0, linetype = 2, color = '#bdbdbd')
p = p + geom_vline(xintercept = 0, linetype = 2, color = '#bdbdbd')
p = p + geom_abline(alpha = 0.1)
p = p + geom_point(aes(color = sig), alpha = 0.3)
p = p + geom_text_repel(aes(label = gene_text, color = sig), min.segment.length = 0)
p = p + theme(legend.position = c(0.6, 0.2))
p = p + xlab('effect size (high coverage, low MOI)')
p = p + ylab('effect size (high coverage, high MOI)')
p = p + scale_color_manual(values = hl_hh_colors)
p
save_plot(p, file = 'img/scatter_hl_hh.pdf', base_height = 7, base_asp = 1)

# now compare to the original screen
original_ssd = readRDS('original_screen_mcmc_ssd.rds')
original_wo = readRDS('wo_original_2019.rds')


original_gene_posterior_summary = gene_posterior_mass(original_wo, original_ssd, level)


hl_original_summary = inner_join(
  hl_gene_posterior_summary,
  original_gene_posterior_summary,
  suffix = c('_hl', '_19'),
  by = c('gene')
)

hl_original_summary = mutate(
  hl_original_summary,
  sig = case_when(
    lfsr_hl < level & lfsr_19 < level ~ 'Both',
    lfsr_hl < level & lfsr_19 >= level ~ 'New',
    lfsr_hl >= level & lfsr_19 < level ~ 'old',
    TRUE ~ 'Neither'
  )
)

hl_original_totals = hl_original_summary %>%
  group_by(sig) %>%
  summarize(total = length(gene))

hl_original_disagreements = filter(hl_original_summary, sig != 'Both' & sig != 'Neither')
hl_original_disagreements = filter(hl_original_disagreements, sign(mu_hl) != sign(mu_19))

hl_original_disagreements %>%
  group_by(sig) %>%
  summarize(sign_disagreement = length(gene))


p = ggplot(hl_original_disagreements, aes(mu_hl, mu_19))
p = p + geom_hline(yintercept = 0, linetype = 2, color = '#bdbdbd', alpha = 0.75)
p = p + geom_vline(xintercept = 0, linetype = 2, color = '#bdbdbd', alpha = 0.75)
p = p + geom_abline(alpha = 0.05)
p = p + geom_pointrange(aes(fill = sig, xmin = lower_hl, xmax = upper_hl),
  color = 'black', shape = 21, alpha = 0.2)
p = p + geom_pointrange(aes(fill = sig, ymin = lower_19, ymax = upper_19), color = 'black', shape = 21, alpha = 0.2)
p = p + theme(legend.position = c(0.6, 0.2))
p = p + xlab('effect size (high coverage, low MOI)')
p = p + ylab('effect size (2019)')
p = p + coord_cartesian(ylim = c(-0.5, 0.5), xlim = c(-0.5, 0.5))
p

tmp = lh_gene_posterior_summary
index = !grepl('gene', colnames(tmp))
colnames(tmp)[index] = paste0(colnames(tmp)[index], '_lh')

compare_with_original = inner_join(
  hl_original_summary,
  tmp,
  suffix = c('_hl2', '_lh'),
  by = c('gene')
)
compare_with_original = select(compare_with_original, -sig)

is_sig = function(x, level = 0.10) {
  as.integer(x < level)
}

ah = mutate_if(compare_with_original,
  grepl('lfsr', colnames(compare_with_original)),
  list(~identity(.), ~is_sig(.)))
upset(dplyr::select(ah, contains('is_sig')), order.by = 'freq')

compare_suffix = function(df, s1, s2, level = 0.10) {
  c1 = paste0('lfsr_', s1)
  c2 = paste0('lfsr_', s2)
  mu1 = paste0('mu_', s1)
  mu2 = paste0('mu_', s2)

  c1_sig = filter(df, !!sym(c1) < level & !!sym(c2) >= level)
  c2_sig = filter(df, !!sym(c2) < level & !!sym(c1) >= level)
  c1_disagreement = filter(df, !!sym(c1) < level & !!sym(c2) >= level &
    sign(!!sym(mu1)) != sign(!!sym(mu2)))
  c2_disagreement = filter(df, !!sym(c2) < level & !!sym(c1) >= level &
    sign(!!sym(mu1)) != sign(!!sym(mu2)))

  disagreement_summary = data.frame(
    type = c('total', 'total', 'disagreement', 'disagreement'),
    value = c(nrow(c1_sig), nrow(c2_sig), nrow(c1_disagreement), nrow(c2_disagreement)),
    group = c(s1, s2, s1, s2)
  )

  list(c1_sig, c2_sig, c1_disagreement, c2_disagreement, disagreement_summary)
}


compare_suffix(compare_with_original, 'lh', '19')[[5]]
compare_suffix(compare_with_original, 'lh', 'hl')[[5]]
compare_suffix(compare_with_original, 'hl', '19')[[5]]



library('MAUDE')

ordering = c('Q1', 'Q2', 'Q3', 'Q4')
bin_labels <- c('A', 'B', 'C', 'D')
maude_bin_mapping = data.frame(ordering = ordering, bin = bin_labels)

bin_names <- c('Low', 'Low_Mid', 'Mid_High', 'High')

high_moi_gfp = dplyr::select(raw_counts, matches('High_MOI_D[1-3]_GFP'))
high_moi_gfp

low_moi_gfp = dplyr::select(raw_counts, matches('low_MOI_D[1-3]_GFP'))
low_moi_gfp



g_counts = gather(raw_counts, key = 'sample', value = 'count', -sgRNA, -Gene)
g_counts = dplyr::rename(g_counts, gid = sgRNA, element = Gene)
g_counts = mutate(g_counts, sample = sub('Teff_', '', sample))
g_counts = filter(g_counts, !grepl('TF_lib_plasmid', sample))
gfp_counts = filter(g_counts, grepl('GFP', sample))
g_counts = filter(g_counts, !grepl('GFP', sample))

g_counts = separate(g_counts,
  sample, c('coverage', 'discard1', 'moi', 'discard2', 'donor', 'ordering'), sep = '_')
g_counts = select(g_counts, !matches('discard'))

g_counts = left_join(g_counts, maude_bin_mapping, by = 'ordering')
g_counts = select(g_counts, -ordering)

gfp_counts = mutate(gfp_counts, sample = sub('_GFP', '', sample))
gfp_counts = separate(gfp_counts,
  sample, c('coverage', 'discard1', 'moi', 'discard2', 'donor'), sep = '_')
gfp_counts = select(gfp_counts, !matches('discard'))
gfp_counts = dplyr::rename(gfp_counts, NS = count)

gfp_counts = bind_rows(gfp_counts,
  dplyr::filter(gfp_counts, moi == 'High') %>% mutate(coverage = 'Low'),
  dplyr::filter(gfp_counts, moi == 'Low') %>% mutate(coverage = 'Low'))


g_counts = mutate(g_counts, NT = grepl('non-targeting', element, ignore.case = TRUE))

s_counts = spread(g_counts, bin, count)

s_counts = left_join(s_counts, select(gfp_counts, -element), by = c('gid', 'coverage', 'moi', 'donor'))

# this is a practice run only one donor and one sample
bin_bounds = makeBinModel(
  data.frame(Bin = bin_labels, fraction = c(0.2, 0.3, 0.3, 0.2)))
bin_bounds = bin_bounds[1:4, ]
bin_bounds = bind_rows(
  mutate(bin_bounds, donor = 'D1'),
  mutate(bin_bounds, donor = 'D2'),
  mutate(bin_bounds, donor = 'D3'))



s_counts_by_experiment = split(s_counts, f = list(s_counts$coverage, s_counts$moi))

# get guide-level stats

maude_results = lapply(s_counts_by_experiment,
  function(df) {
    print(df[1, ])
    guideLevelStats <- findGuideHitsAllScreens(experiments = unique(df["donor"]),
      countDataFrame = df,
      binStats = bin_bounds,
      # unsortedBin = "notSorted",
      sortBins = bin_labels)
    print('guide done')
    elementLevelStats <- getElementwiseStats(unique(guideLevelStats["donor"]),
      guideLevelStats,
      elementIDs = "element")
    maude_results = summarize(group_by(elementLevelStats, element),
      p_value = combine.test(p.value), mean_effect = mean(meanZ))
    maude_results = mutate(maude_results, fdr = p.adjust(p_value, method = 'BH'))
    list(guide_level_stats = guideLevelStats, element_level_stats = elementLevelStats,
      maude_results = maude_results)
  })

saveRDS(maude_results, file = 'all_moi_maude.rds')

maude_results = readRDS('all_moi_maude.rds')

maude_validated = inner_join(maude_results[['Low.High']]$maude_results, validated_targets,
  by = c('element' = 'gene'))

maude_validated = mutate(maude_validated,
  status = case_when(
    sign(mean_effect) == validation_sign & fdr < 0.10 ~ 'true positive',
    fdr > 0.10 & validation_sign == 0 ~ 'true negative',
    fdr < 0.10 & validation_sign == 0 ~ 'false positive',
    fdr > 0.10 & validation_sign != 0 ~ 'false negative',
  ))


p = ggplot(maude_validated, aes(mfi_rank_centered, mean_effect))
# p = ggplot(maude_validated, aes(validation_direction, mean_effect))
p = p + geom_hline(yintercept = 0, linetype = 2, color = '#bdbdbd')
p = p + geom_vline(xintercept = 0, linetype = 2, color = '#bdbdbd')
p = p + geom_point(aes(color = status, shape = status),
  position = position_dodge2(width = 1),
  size = 4,
  stroke = 1
  )
# p = p + geom_pointrange(aes(ymin = lower, ymax = upper,
#     color = status, shape = status),
#   position = position_dodge2(width = 0.75)
#   )
p = p + scale_x_discrete()
p = p + xlab('Arrayed knockout direction')
p = p + ylab('MAUDE gene meanZ')
p = p + scale_color_manual(values = color_mapping)
p = p + scale_shape_manual(values = shape_mapping)
p = p + geom_text_repel(aes(label = element, color = status), min.segment.length = 0,
  max.overlaps = Inf, box.padding = 0.75,
  position = position_dodge2(width = 1)
)
p = p + theme(legend.position = c(0.70, 0.1), legend.title = element_blank())
p

save_plot(p, file = 'img/lh_maude_revised.pdf', base_height = 7, base_asp = 1)

nrow(filter(maude_results[['Low.High']]$maude_results, fdr < 0.10))

# validate waterbear with fewer replicates
#

DATA_HOME = '/oak/stanford/groups/pritch/users/jake/waterbear/high_moi/count_files'
raw_counts = read.table(paste0(DATA_HOME, '/high_low_MOI_screens_D1_D2_D3_2022-04-23.count.txt'), header = TRUE,
                        stringsAsFactors = FALSE, sep = '\t')


ordering = c('Q1', 'Q2', 'Q3', 'Q4')

c_name = grep('Low_Coverage_High_MOI', colnames(raw_counts), value = TRUE)
sample_mapping = data.frame(c_name = c_name)
sample_mapping = mutate(sample_mapping, sample = str_extract(c_name, 'D1|D2|D3'))
sample_mapping = mutate(sample_mapping, bin = str_extract(c_name, 'Q[1-4]'))

sample_mapping = filter(sample_mapping, sample != 'D3')

wb_array = wb_counts_to_array(raw_counts, sample_mapping, ordering)

gene_mapping = raw_counts[, c('sgRNA', 'Gene')]
colnames(gene_mapping) = c('guide', 'gene')

wo = wb_make_object(wb_array, gene_mapping, bin_size_prior = c(0.2, 0.3, 0.3, 0.2))

# debugonce(wb_make_object)
# wo_tmp = wb_make_object(array, gene_mapping, bin_size_prior = c(0.2, 0.3, 0.3, 0.2))


#debugonce(wb_em_start)
wo = wb_em_start(wo)
#wo = wo$wo

source('../waterbear/R/nimble.R')
n_model = nimbleModel(sample_specific_dispersion_model,
                      data = wo$data, constants = wo$const, inits = wo$init)

# nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
# nimbleOptions(MCMCsaveHistory = TRUE)

n_configuration = configureMCMC(n_model)
# sampler_control = list(order = wo$const$order)
# n_configuration$removeSamplers('gene_inclusion')
# n_configuration$addSampler(target = 'gene_inclusion', type = 'rank_RW',
#                            control = sampler_control)

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
n_chains = 4
seed = 42
system.time({
  samples = runMCMC(
    C_n_mcmc, niter = n_samples, nburnin = n_burnin,
    nchains = n_chains,
    setSeed = seed,
    summary = TRUE)
})

save(wo, samples, file = 'low_coverage_high_moi_d1_d2.RData')
saveRDS(samples, 'low_coverage_high_moi_d1_d2.rds')
saveRDS(samples, 'new_moi_d1_d2.rds')


# 2 donors, high coverge, low MOI

c_name = grep('High_Coverage_Low_MOI', colnames(raw_counts), value = TRUE)
sample_mapping = data.frame(c_name = c_name)
sample_mapping = mutate(sample_mapping, sample = str_extract(c_name, 'D1|D2|D3'))
sample_mapping = mutate(sample_mapping, bin = str_extract(c_name, 'Q[1-4]'))
sample_mapping = sample_mapping[complete.cases(sample_mapping), ]
sample_mapping = filter(sample_mapping, sample == 'D1' | sample == 'D2')

wb_array = wb_counts_to_array(raw_counts, sample_mapping, ordering)

gene_mapping = raw_counts[, c('sgRNA', 'Gene')]
colnames(gene_mapping) = c('guide', 'gene')

wo = wb_make_object(wb_array, gene_mapping, bin_size_prior = c(0.2, 0.3, 0.3, 0.2))

# debugonce(wb_make_object)
# wo_tmp = wb_make_object(array, gene_mapping, bin_size_prior = c(0.2, 0.3, 0.3, 0.2))


#debugonce(wb_em_start)
wo = wb_em_start(wo)

source('../waterbear/R/nimble.R')
n_model = nimbleModel(sample_specific_dispersion_model,
                      data = wo$data, constants = wo$const, inits = wo$init)

n_configuration = configureMCMC(n_model)

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
n_chains = 4
seed = 42
system.time({
  samples = runMCMC(
    C_n_mcmc, niter = n_samples, nburnin = n_burnin,
    nchains = n_chains,
    setSeed = seed,
    summary = TRUE)
})
save(wo, samples, file = 'high_coverage_low_moi_d1_d2.RData')


c_name = grep('Low_Coverage_Low_MOI', colnames(raw_counts), value = TRUE)
sample_mapping = data.frame(c_name = c_name)
sample_mapping = mutate(sample_mapping, sample = str_extract(c_name, 'D1|D2|D3'))
sample_mapping = mutate(sample_mapping, bin = str_extract(c_name, 'Q[1-4]'))
sample_mapping = sample_mapping[complete.cases(sample_mapping), ]
sample_mapping = filter(sample_mapping, sample == 'D1' | sample == 'D2')
wb_array = wb_counts_to_array(raw_counts, sample_mapping, ordering)

gene_mapping = raw_counts[, c('sgRNA', 'Gene')]
colnames(gene_mapping) = c('guide', 'gene')

wo = wb_make_object(wb_array, gene_mapping, bin_size_prior = c(0.2, 0.3, 0.3, 0.2))
wo = wb_em_start(wo)

source('../waterbear/R/nimble.R')
n_model = nimbleModel(sample_specific_dispersion_model,
                      data = wo$data, constants = wo$const, inits = wo$init)

n_configuration = configureMCMC(n_model)

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
n_chains = 4
seed = 42
system.time({
  samples = runMCMC(
    C_n_mcmc, niter = n_samples, nburnin = n_burnin,
    nchains = n_chains,
    setSeed = seed,
    summary = TRUE)
})
save(wo, samples, file = 'low_coverage_low_moi_d1_d2.RData')


c_name = grep('High_Coverage_High_MOI', colnames(raw_counts), value = TRUE)
sample_mapping = data.frame(c_name = c_name)
sample_mapping = mutate(sample_mapping, sample = str_extract(c_name, 'D1|D2|D3'))
sample_mapping = mutate(sample_mapping, bin = str_extract(c_name, 'Q[1-4]'))
sample_mapping = sample_mapping[complete.cases(sample_mapping), ]
sample_mapping = filter(sample_mapping, sample == 'D1' | sample == 'D2')
wb_array = wb_counts_to_array(raw_counts, sample_mapping, ordering)

gene_mapping = raw_counts[, c('sgRNA', 'Gene')]
colnames(gene_mapping) = c('guide', 'gene')

wo = wb_make_object(wb_array, gene_mapping, bin_size_prior = c(0.2, 0.3, 0.3, 0.2))
wo = wb_em_start(wo)

source('../waterbear/R/nimble.R')
n_model = nimbleModel(sample_specific_dispersion_model,
                      data = wo$data, constants = wo$const, inits = wo$init)

n_configuration = configureMCMC(n_model)

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
n_chains = 4
seed = 42
system.time({
  samples = runMCMC(
    C_n_mcmc, niter = n_samples, nburnin = n_burnin,
    nchains = n_chains,
    setSeed = seed,
    summary = TRUE)
})
save(wo, samples, file = 'high_coverage_high_moi_d1_d2.RData')

# d2 d3

DATA_HOME = '/oak/stanford/groups/pritch/users/jake/waterbear/high_moi/count_files'
raw_counts = read.table(
  paste0(DATA_HOME, '/high_low_MOI_screens_D1_D2_D3_2022-04-23.count.txt'), header = TRUE,
                        stringsAsFactors = FALSE, sep = '\t')
ordering = c('Q1', 'Q2', 'Q3', 'Q4')

c_name = grep('Low_Coverage_High_MOI', colnames(raw_counts), value = TRUE)
sample_mapping = data.frame(c_name = c_name)
sample_mapping = mutate(sample_mapping, sample = str_extract(c_name, 'D1|D2|D3'))
sample_mapping = mutate(sample_mapping, bin = str_extract(c_name, 'Q[1-4]'))
sample_mapping = filter(sample_mapping, sample != 'D1')

wb_array = wb_counts_to_array(raw_counts, sample_mapping, ordering)

gene_mapping = raw_counts[, c('sgRNA', 'Gene')]
colnames(gene_mapping) = c('guide', 'gene')

wo = wb_make_object(wb_array, gene_mapping, bin_size_prior = c(0.2, 0.3, 0.3, 0.2))

wo = wb_em_start(wo)

source('../waterbear/R/nimble.R')
n_model = nimbleModel(sample_specific_dispersion_model,
                      data = wo$data, constants = wo$const, inits = wo$init)

n_configuration = configureMCMC(n_model)

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


n_samples = 10000
n_burnin = 5000
n_chains = 4
seed = 42
system.time({
  samples = runMCMC(
    C_n_mcmc, niter = n_samples, nburnin = n_burnin,
    nchains = n_chains,
    setSeed = seed,
    summary = TRUE)
})

save(wo, samples, file = 'low_coverage_high_moi_d2_d3.RData')
# saveRDS(samples, 'low_coverage_high_moi_d2_d3.rds')
# saveRDS(samples, 'new_moi_d1_d2.rds')


# 2 donors, high coverge, low MOI

c_name = grep('High_Coverage_Low_MOI', colnames(raw_counts), value = TRUE)
sample_mapping = data.frame(c_name = c_name)
sample_mapping = mutate(sample_mapping, sample = str_extract(c_name, 'D1|D2|D3'))
sample_mapping = mutate(sample_mapping, bin = str_extract(c_name, 'Q[1-4]'))
sample_mapping = sample_mapping[complete.cases(sample_mapping), ]
sample_mapping = filter(sample_mapping, sample == 'D2' | sample == 'D3')

wb_array = wb_counts_to_array(raw_counts, sample_mapping, ordering)

gene_mapping = raw_counts[, c('sgRNA', 'Gene')]
colnames(gene_mapping) = c('guide', 'gene')

wo = wb_make_object(wb_array, gene_mapping, bin_size_prior = c(0.2, 0.3, 0.3, 0.2))

# debugonce(wb_make_object)
# wo_tmp = wb_make_object(array, gene_mapping, bin_size_prior = c(0.2, 0.3, 0.3, 0.2))


#debugonce(wb_em_start)
wo = wb_em_start(wo)

source('../waterbear/R/nimble.R')
n_model = nimbleModel(sample_specific_dispersion_model,
                      data = wo$data, constants = wo$const, inits = wo$init)

n_configuration = configureMCMC(n_model)

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
n_chains = 4
seed = 42
system.time({
  samples = runMCMC(
    C_n_mcmc, niter = n_samples, nburnin = n_burnin,
    nchains = n_chains,
    setSeed = seed,
    summary = TRUE)
})
save(wo, samples, file = 'high_coverage_low_moi_d2_d3.RData')


c_name = grep('Low_Coverage_Low_MOI', colnames(raw_counts), value = TRUE)
sample_mapping = data.frame(c_name = c_name)
sample_mapping = mutate(sample_mapping, sample = str_extract(c_name, 'D1|D2|D3'))
sample_mapping = mutate(sample_mapping, bin = str_extract(c_name, 'Q[1-4]'))
sample_mapping = sample_mapping[complete.cases(sample_mapping), ]
sample_mapping = filter(sample_mapping, sample == 'D2' | sample == 'D3')
wb_array = wb_counts_to_array(raw_counts, sample_mapping, ordering)

gene_mapping = raw_counts[, c('sgRNA', 'Gene')]
colnames(gene_mapping) = c('guide', 'gene')

wo = wb_make_object(wb_array, gene_mapping, bin_size_prior = c(0.2, 0.3, 0.3, 0.2))
wo = wb_em_start(wo)

source('../waterbear/R/nimble.R')
n_model = nimbleModel(sample_specific_dispersion_model,
                      data = wo$data, constants = wo$const, inits = wo$init)

n_configuration = configureMCMC(n_model)

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
n_chains = 4
seed = 42
system.time({
  samples = runMCMC(
    C_n_mcmc, niter = n_samples, nburnin = n_burnin,
    nchains = n_chains,
    setSeed = seed,
    summary = TRUE)
})
save(wo, samples, file = 'low_coverage_low_moi_d2_d3.RData')


c_name = grep('High_Coverage_High_MOI', colnames(raw_counts), value = TRUE)
sample_mapping = data.frame(c_name = c_name)
sample_mapping = mutate(sample_mapping, sample = str_extract(c_name, 'D1|D2|D3'))
sample_mapping = mutate(sample_mapping, bin = str_extract(c_name, 'Q[1-4]'))
sample_mapping = sample_mapping[complete.cases(sample_mapping), ]
sample_mapping = filter(sample_mapping, sample == 'D2' | sample == 'D3')
wb_array = wb_counts_to_array(raw_counts, sample_mapping, ordering)

gene_mapping = raw_counts[, c('sgRNA', 'Gene')]
colnames(gene_mapping) = c('guide', 'gene')

wo = wb_make_object(wb_array, gene_mapping, bin_size_prior = c(0.2, 0.3, 0.3, 0.2))
wo = wb_em_start(wo)

source('../waterbear/R/nimble.R')
n_model = nimbleModel(sample_specific_dispersion_model,
                      data = wo$data, constants = wo$const, inits = wo$init)

n_configuration = configureMCMC(n_model)

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
n_chains = 4
seed = 42
system.time({
  samples = runMCMC(
    C_n_mcmc, niter = n_samples, nburnin = n_burnin,
    nchains = n_chains,
    setSeed = seed,
    summary = TRUE)
})
save(wo, samples, file = 'high_coverage_high_moi_d2_d3.RData')

# analyze the down sampling results
hchm_d1d2 = new.env()
load('high_coverage_high_moi_d1_d2.RData', hchm_d1d2)

hchm_d1d2$samples$summary$all.chains[grepl('dispersion', rownames(hchm_d1d2$samples$summary$all.chains)), ]
hchm_d1d2$samples$summary$all.chains[grepl('psi', rownames(hchm_d1d2$samples$summary$all.chains)), ]
hchm_d1d2$samples$summary$all.chains[grepl('sigma', rownames(hchm_d1d2$samples$summary$all.chains)), ]

level = 0.10

hchm_d1d2$guide_summaries = get_guide_summaries(hchm_d1d2$wo, hchm_d1d2$samples)
hchm_d1d2$gene_posterior_summary = gene_posterior_mass(hchm_d1d2$wo, hchm_d1d2$samples, level)

gene_columns = grep('gene_shift', colnames(hchm_d1d2$samples$samples$chain1))

interesting_genes = data.frame(gene =
  c('IL2RA',
    'STAT5B',
    'MYB',
    'JAK3',
    'MED3',
    'CBFB'))

hchm_d1d2$gene_inclusion = hchm_d1d2$guide_summaries$gene_inclusion

hchm_d1d2$gene_shift = hchm_d1d2$guide_summaries$gene_shift
hchm_d1d2$gene_shift = inner_join(hchm_d1d2$gene_shift, hchm_d1d2$gene_posterior_summary, by = 'gene')
hchm_d1d2$gene_shift = inner_join(hchm_d1d2$gene_shift,
  dplyr::select(hchm_d1d2$gene_inclusion, gene, pip = Mean), by = 'gene')

# TODO: run the same summary statistics for the other configurations
# TODO: figure out how to compare these effect sizes and perhaps the uncertainty around them

tmp = inner_join(hchm_d1d2$gene_posterior_summary, hh_gene_posterior_summary, suffix = c('_d1d2', '_all'), by = 'gene')

p = ggplot(tmp, aes(mu_d1d2, mu_all))
p = p + geom_point(alpha = 0.1)
save_plot(p, file = 'tmp.pdf')

hchm_d1d2$validated = inner_join(hchm_d1d2$gene_posterior_summary, validated_targets, by = 'gene')
hchm_d1d2$validated = arrange(hchm_d1d2$validated, mu)
hchm_d1d2$validated = mutate(hchm_d1d2$validated, posterior_rank = 1:length(mu))

hchm_d1d2$validated = mutate(hchm_d1d2$validated,
  status = case_when(
    sign(mu) == validation_sign & lfsr < 0.10 ~ 'true positive',
    lfsr > 0.10 & validation_sign == 0 ~ 'true negative',
    lfsr < 0.10 & validation_sign == 0 ~ 'false positive',
    lfsr > 0.10 & validation_sign != 0 ~ 'false negative',
  ))



analyze_configuration = function(env, validated_targets) {
  env$samples$summary$all.chains[grepl('dispersion', rownames(env$samples$summary$all.chains)), ]
  env$samples$summary$all.chains[grepl('psi', rownames(env$samples$summary$all.chains)), ]
  env$samples$summary$all.chains[grepl('sigma', rownames(env$samples$summary$all.chains)), ]

  level = 0.10

  env$guide_summaries = get_guide_summaries(env$wo, env$samples)
  env$gene_posterior_summary = gene_posterior_mass(env$wo, env$samples, level)

  gene_columns = grep('gene_shift', colnames(env$samples$samples$chain1))

  interesting_genes = data.frame(gene =
    c('IL2RA',
      'STAT5B',
      'MYB',
      'JAK3',
      'MED3',
      'CBFB'))

  env$gene_inclusion = env$guide_summaries$gene_inclusion

  env$gene_shift = env$guide_summaries$gene_shift
  env$gene_shift = inner_join(env$gene_shift, env$gene_posterior_summary, by = 'gene')
  env$gene_shift = inner_join(env$gene_shift,
    dplyr::select(env$gene_inclusion, gene, pip = Mean), by = 'gene')
  env$gene_significant = env$gene_shift <= level

  # TODO: run the same summary statistics for the other configurations
  # TODO: figure out how to compare these effect sizes and perhaps the uncertainty around them

  # tmp = inner_join(env$gene_posterior_summary, hh_gene_posterior_summary, suffix = c('_d1d2', '_all'), by = 'gene')

  # p = ggplot(tmp, aes(mu_d1d2, mu_all))
  # p = p + geom_point(alpha = 0.1)
  # save_plot(p, file = 'tmp.pdf')

  env$validated = inner_join(env$gene_posterior_summary, validated_targets, by = 'gene')
  env$validated = arrange(env$validated, mu)
  env$validated = mutate(env$validated, posterior_rank = 1:length(mu))

  env$validated = mutate(env$validated,
    status = case_when(
      sign(mu) == validation_sign & lfsr < 0.10 ~ 'true positive',
      lfsr > 0.10 & validation_sign == 0 ~ 'true negative',
      lfsr < 0.10 & validation_sign == 0 ~ 'false positive',
      lfsr > 0.10 & validation_sign != 0 ~ 'false negative',
      ))

  env$sensitivity = sum(env$validated$status == 'true positive') /
    sum(env$validated$validation_sign != 0)
}

hclm_d1d2 = new.env()
load('high_coverage_low_moi_d1_d2.RData', hclm_d1d2)

lclm_d1d2 = new.env()
load('low_coverage_low_moi_d1_d2.RData', lclm_d1d2)

lchm_d1d2 = new.env()
load('low_coverage_high_moi_d1_d2.RData', lchm_d1d2)

hchm_d1d2 = new.env()
load('high_coverage_high_moi_d1_d2.RData', hchm_d1d2)


analyze_configuration(hclm_d1d2, validated_targets)
analyze_configuration(lclm_d1d2, validated_targets)
analyze_configuration(lchm_d1d2, validated_targets)
analyze_configuration(hchm_d1d2, validated_targets)

analyze_configuration(low_high, validated_targets)
analyze_configuration(high_low, validated_targets)
analyze_configuration(low_low, validated_targets)
analyze_configuration(high_high, validated_targets)

nrow(filter(low_high$gene_posterior_summary, lfsr < 0.10))

data.frame(
  coverage = c('high', 'low', 'low', 'high'),
  moi = c('low', 'low', 'high', 'high')
)

lapply(list(
    hclm_d1d2,
    lclm_d1d2,
    lchm_d1d2,
    hchm_d1d2
    ),
  function(env) env$sensitivity)


lapply(list(
    low_high,
    high_low,
    low_low,
    high_high
    ),
  function(env) env$sensitivity)


significant2 = lapply(
  list(
    hclm_d1d2,
    lclm_d1d2,
    lchm_d1d2,
    hchm_d1d2,
    low_high,
    high_low,
    low_low,
    high_high
  ),
  function(env, labels) {
    df = env$gene_shift
    df = arrange(df, gene)
    res = df$lfsr <= 0.10
    res
  }
)

significant2 = do.call('cbind', significant2)
colnames(significant2) = c(
    'hclm_d1d2',
    'lclm_d1d2',
    'lchm_d1d2',
    'hchm_d1d2',
    'hclm',
    'lclm',
    'lchm',
    'hchm'
    )

tmp = cbind(significant2[[1]], significant2[[2]])


library('UpSetR')
library('ComplexUpset')
dir.create('img')

mode(significant2) = 'integer'
class(significant2) = 'array'

# upset(all_significant, sets = c('low_low', 'low_high', 'high_low', 'high_high'))
pdf('img/upset_downsample.pdf', width = 7, height = 4)
significant2 = data.frame(significant2)
UpSetR::upset(significant2,
  NA,
  intersections =
    list(
      c('hclm'),
      c('hclm_d1d2'),
      c('hclm_d1d2', 'hclm'),
      c('lclm'),
      c('lclm_d1d2'),
      c('lclm_d1d2', 'lclm'),
      c('lchm'),
      c('lchm_d1d2'),
      c('lchm_d1d2', 'lchm'),
      c('hchm_d1d2'),
      c('hchm'),
      c('hchm_d1d2', 'hchm'),
      c(
        'hclm_d1d2',
        'lclm_d1d2',
        'lchm_d1d2',
        'hchm_d1d2',
        'hclm',
        'lclm',
        'lchm',
        'hchm')
    ),
  keep.order = TRUE
)
dev.off()

tmp = inner_join(
  dplyr::select(high_low$gene_shift, Mean, lfsr, gene),
  dplyr::select(hclm_d1d2$gene_shift, Mean, lfsr, gene),
  by = 'gene',
  suffix = c('_full', '_two')
  )

level = 0.10

tmp_summary = mutate(
  tmp,
  sig = case_when(
    lfsr_full < level & lfsr_two < level ~ 'Both',
    lfsr_full < level & lfsr_two >= level ~ 'Three replicates',
    lfsr_full >= level & lfsr_two < level ~ 'Two replicates',
    TRUE ~ 'Neither'
  )
)

tmp_summary_counts = tmp_summary %>% group_by(sig) %>% summarize(n = length(gene))
tmp_summary_counts = mutate(tmp_summary_counts, label_text = sort(paste0(sig, ' (', n, ')')),
  color = c('#009E73', '#999999', '#0072B2', '#E69F00'))
names(tmp_summary_counts$color) = tmp_summary_counts$sig
names(tmp_summary_counts$label_text) = tmp_summary_counts$sig

p = ggplot(tmp_summary, aes(Mean_full, Mean_two))
p = p + geom_hline(yintercept = 0, linetype = 2, color = '#bdbdbd', alpha = 0.75)
p = p + geom_vline(xintercept = 0, linetype = 2, color = '#bdbdbd', alpha = 0.75)
p = p + geom_abline(alpha = 0.05)
p = p + geom_point(aes(fill = sig, alpha = sig, size = sig), color = 'black', shape = 21)
# p = p + geom_text_repel(aes(label = gene_text, color = sig), min.segment.length = 0)
p = p + theme(legend.position = c(0.6, 0.2))
p = p + xlab('effect size (Three replicates)')
p = p + ylab('effect size (Two replicates)')
p = p + scale_fill_manual(values = tmp_summary_counts$color, labels = tmp_summary_counts$label_text)
p = p + scale_color_manual(values = tmp_summary_counts$color, labels = tmp_summary_counts$label_text)
p = p + scale_size_manual(values = c(2, 0.4, 2, 2), labels = tmp_summary_counts$label_text)
p = p + scale_alpha_manual(values = c(0.5, 0.1, 0.5, 0.5), labels = tmp_summary_counts$label_text)
p
save_plot(p, file = 'img/scatter_hclm_2_full.pdf', base_height = 7, base_asp = 1)

plot_two_three_replicates = function(two_reps, three_reps) {
  tmp = inner_join(
    dplyr::select(three_reps$gene_shift, Mean, lfsr, gene),
    dplyr::select(two_reps$gene_shift, Mean, lfsr, gene),
    by = 'gene',
    suffix = c('_full', '_two')
  )

  level = 0.10

  tmp_summary = mutate(
    tmp,
    sig = case_when(
      lfsr_full < level & lfsr_two < level ~ 'Both',
      lfsr_full < level & lfsr_two >= level ~ 'Three replicates',
      lfsr_full >= level & lfsr_two < level ~ 'Two replicates',
      TRUE ~ 'Neither'
    )
  )

  tmp_summary_counts = tmp_summary %>% group_by(sig) %>% summarize(n = length(gene))
  tmp_summary_counts = mutate(tmp_summary_counts, label_text = sort(paste0(sig, ' (', n, ')')),
    color = c('#009E73', '#999999', '#0072B2', '#E69F00'))
  names(tmp_summary_counts$color) = tmp_summary_counts$sig
  names(tmp_summary_counts$label_text) = tmp_summary_counts$sig

  p = ggplot(tmp_summary, aes(Mean_full, Mean_two))
  p = p + geom_hline(yintercept = 0, linetype = 2, color = '#bdbdbd', alpha = 0.75)
  p = p + geom_vline(xintercept = 0, linetype = 2, color = '#bdbdbd', alpha = 0.75)
  p = p + geom_abline(alpha = 0.05)
  p = p + geom_point(aes(fill = sig, alpha = sig, size = sig), color = 'black', shape = 21)
  # p = p + geom_text_repel(aes(label = gene_text, color = sig), min.segment.length = 0)
  p = p + theme(legend.position = c(0.6, 0.2))
  p = p + xlab('effect size (Three replicates)')
  p = p + ylab('effect size (Two replicates)')
  p = p + scale_fill_manual(values = tmp_summary_counts$color, labels = tmp_summary_counts$label_text)
  p = p + scale_color_manual(values = tmp_summary_counts$color, labels = tmp_summary_counts$label_text)
  p = p + scale_size_manual(values = c(2, 0.4, 2, 2), labels = tmp_summary_counts$label_text)
  p = p + scale_alpha_manual(values = c(0.5, 0.1, 0.5, 0.5), labels = tmp_summary_counts$label_text)
  p
}

p = plot_two_three_replicates(high_low, hclm_d1d2)
save_plot(p, file = 'img/scatter_hclm_d1d2.pdf', base_height = 7, base_asp = 1)

p = plot_two_three_replicates(low_high, lchm_d1d2)
save_plot(p, file = 'img/scatter_lchm_d1d2.pdf', base_height = 7, base_asp = 1)

p = plot_two_three_replicates(low_low, lclm_d1d2)
save_plot(p, file = 'img/scatter_lclm_d1d2.pdf', base_height = 7, base_asp = 1)

p = plot_two_three_replicates(high_high, hchm_d1d2)
save_plot(p, file = 'img/scatter_hchm_d1d2.pdf', base_height = 7, base_asp = 1)
