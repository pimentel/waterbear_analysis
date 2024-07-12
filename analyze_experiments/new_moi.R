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



# setwd('~/Desktop/waterbear')

raw_counts = read.table('../data/high_MOI_D1_D2_D3_2022-03-02.count.txt', header = TRUE,
                        stringsAsFactors = FALSE, sep = '\t')


ordering = c('Q1', 'Q2', 'Q3', 'Q4')

c_name = grep('HighMOI_Q', colnames(raw_counts), value = TRUE)
sample_mapping = data.frame(c_name = c_name)
sample_mapping = mutate(sample_mapping, sample = str_extract(c_name, 'D1|D2|D3'))
sample_mapping = mutate(sample_mapping, bin = sub('D._HighMOI_', '', c_name))

array = wb_counts_to_array(raw_counts, sample_mapping, ordering)

gene_mapping = raw_counts[, c('sgRNA', 'Gene')]
colnames(gene_mapping) = c('guide', 'gene')

wo = wb_make_object(array, gene_mapping, bin_size_prior = c(0.2, 0.3, 0.3, 0.2))

debugonce(wb_make_object)
wo_tmp = wb_make_object(array, gene_mapping, bin_size_prior = c(0.2, 0.3, 0.3, 0.2))


#debugonce(wb_em_start)
wo = wb_em_start(wo)
#wo = wo$wo

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
n_chains = 4
seed = 42
system.time({
  samples = runMCMC(
    C_n_mcmc, niter = n_samples, nburnin = n_burnin,
    nchains = n_chains,
    setSeed = seed,
    summary = TRUE)
})

saveRDS(samples, 'new_moi.rds')

samples = readRDS('new_moi.rds')



get_guide_summaries = function(wo, samples) {
  sample_summary = samples$summary$all.chains
  guide_to_guide_mapping = data.frame(guide_data_index = wo$const$guide_data_index,
    wb_guide_index = wo$const$guide_index)
  sample_summary = data.frame(variable = rownames(sample_summary), sample_summary)
  guides = dplyr::filter(sample_summary, grepl('guide_shift', variable))
  guides = dplyr::mutate(guides, wb_guide_index = as.integer(str_extract(variable, '(\\d)+')))
  guides = inner_join(guides, guide_to_guide_mapping, by = 'wb_guide_index')
  # XXX: this mapping is not correct
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
      gt0 = mean(tmp > .Machine$double.eps)
      lt0 = mean(tmp < .Machine$double.eps)
      data.frame(nimble_id = gene_name, lower = qs[1], upper = qs[2],
        gt0 = gt0, lt0 = lt0, row.names = NULL)
    })
  gene_summary = dplyr::bind_rows(gene_summary)
  gene_summary = dplyr::mutate(gene_summary, mapping = as.integer(stringr::str_extract(nimble_id, '(\\d)+')))
  gene_summary = inner_join(gene_summary,
    dplyr::distinct(dplyr::select(wo$test_guide_names, -c(i, guide))),
    by = c('mapping' = 'mapping'))
  gene_summary = dplyr::mutate(gene_summary, lfsr = pmax(gt0, lt0))
  dplyr::select(gene_summary, -nimble_id, -mapping)
}

gene_posterior_summary = gene_posterior_mass(wo, samples, 0.05)

gene_columns = grep('gene_shift', colnames(samples$samples$chain1))

chain = samples$samples$chain1

tmp = sapply(samples$samples,
  function(chain) {
    chain[, gene_columns[1]]
  }, simplify = FALSE)

# debugonce(get_guide_summaries)
guide_summaries = get_guide_summaries(wo, samples)

interesting_genes = data.frame(gene =
  c('IL2RA',
    'STAT5B',
    'MYB',
    'JAK3',
    'MED3',
    'CBFB'))



# top_genes = filter(dplyr::arrange(guide_summaries$gene_inclusion, desc(Mean)), gene == 'IL2RA')
gene_inclusion = guide_summaries$gene_inclusion



gene_shift = guide_summaries$gene_shift
gene_shift = inner_join(gene_shift, gene_posterior_summary, by = 'gene')
gene_shift = inner_join(gene_shift,
  dplyr::select(gene_inclusion, gene, pip = Mean), by = 'gene')

# gene_shift = inner_join(guide_summaries$gene_shift, interesting_genes, by = 'gene')
guide_shift = inner_join(guide_summaries$guides, interesting_genes, by = 'gene')

level = 0.05

gene_shift = mutate(gene_shift, high_pip = pip > 1 - level)

sig_genes_wb = dplyr::select(filter(gene_shift, high_pip & lfsr > 1 - level / 2), gene)
sig_genes_wb =  dplyr::mutate(sig_genes_wb, wb_sig = TRUE)

p = ggplot(gene_shift, aes(lfsr, Mean, color = pip))
p = p + geom_point(alpha = 0.4)
p = p + geom_vline(xintercept = c(level / 2, 1 - level / 2), linetype = 2)
p = p + xlab('proportion of gene_shift > 0')
p = p + ylab('gene shift')
save_plot(p, file = 'lfsr_mean.pdf')

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
