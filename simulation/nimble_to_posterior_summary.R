library('dplyr')
library('tidyr')
library('stringr')

# used for debugging
fnames = Sys.glob('results/nimble/unequal_*equalmixture*.rds')
fnames = fnames[!grepl('summary', fnames)]
current_base = sub('.rds', '', basename(fnames[1]))
sample_fname = fnames[1]

# and not summary!
current_base = basename(sub('.rds', '', fnames[grepl('_4_205000_', fnames)][1]))

sim_nimble = readRDS(paste0('results/truth/nimble_', current_base, '.rds'))
samples_n = readRDS(paste0('results/nimble/', current_base, '.rds'))

# END: debugging

args = commandArgs(trailingOnly = TRUE)

sample_fname = args[1]

source('benchmark_common.R')

if (!file.exists(sample_fname)) {
  stop(paste0("File doesn't exist: ", sample_fname))
}
samples_n = readRDS(sample_fname)

get_guide_summaries = function(samples) {
  sample_summary = samples$summary$all.chains
  # guide_to_guide_mapping = data.frame(guide_data_index = wo$const$guide_data_index,
  # guide_to_guide_mapping = data.frame(guide_data_index = wo$const$guide_data_index,
    # wb_guide_index = wo$const$guide_index)
  sample_summary = data.frame(variable = rownames(sample_summary), sample_summary)
  guides = dplyr::filter(sample_summary, grepl('guide_shift', variable))
  guides = dplyr::mutate(guides, wb_guide_index = as.integer(str_extract(variable, '(\\d)+')))
  # guides = inner_join(guides, guide_to_guide_mapping, by = 'wb_guide_index')
  # XXX: this mapping is not correct
  # guides = inner_join(guides, wo$test_guide_names, by = c('guide_data_index' = 'i'))
  gene_inclusion = dplyr::filter(sample_summary, grepl('gene_inclusion', variable))
  gene_inclusion = dplyr::mutate(gene_inclusion,
    gene_mapping = as.integer(str_extract(variable, '(\\d)+')))
  # gene_inclusion = inner_join(gene_inclusion,
  #   dplyr::distinct(dplyr::select(wo$test_guide_names, -c(i, guide))),
  #   by = c('gene_mapping' = 'mapping'))
  gene_shift = dplyr::filter(sample_summary, grepl('gene_shift', variable))
  gene_shift = dplyr::mutate(gene_shift,
    gene_mapping = as.integer(str_extract(variable, '(\\d)+')))
  # gene_shift = inner_join(gene_shift,
  #   dplyr::distinct(dplyr::select(wo$test_guide_names, -c(i, guide))),
  #   by = c('gene_mapping' = 'mapping'))

  list(guides = guides, gene_inclusion = gene_inclusion, gene_shift = gene_shift)
}


# debugonce(get_guide_summaries)
# tmp = get_guide_summaries(samples_n)

gene_posterior_mass = function(samples, alpha) {
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
      gt0 = mean(tmp >= 0)
      lt0 = mean(tmp <= 0)
      # gt0 = mean(tmp > .Machine$double.eps)
      # lt0 = mean(tmp < .Machine$double.eps)
      data.frame(nimble_id = gene_name, lower = qs[1], upper = qs[2],
        gt0 = gt0, lt0 = lt0, row.names = NULL)
    })
  gene_summary = dplyr::bind_rows(gene_summary)
  gene_summary = dplyr::mutate(gene_summary, mapping = as.integer(stringr::str_extract(nimble_id, '(\\d)+')))
  # gene_summary = inner_join(gene_summary,
  #   dplyr::distinct(dplyr::select(wo$test_guide_names, -c(i, guide))),
  #   by = c('mapping' = 'mapping'))
  gene_summary = dplyr::mutate(gene_summary, lfsr = pmin(gt0, lt0))
  # dplyr::select(gene_summary, -nimble_id, -mapping)
}

gene_posterior_summary = function(samples, alpha) {
  guide_summaries = get_guide_summaries(samples)
  posterior_mass = gene_posterior_mass(samples, alpha)
  merged = inner_join(
    select(mutate(guide_summaries$gene_inclusion, value = 1 - Mean), gene_id = gene_mapping, value),
    select(posterior_mass, gene_id = mapping, lfsr),
    by = 'gene_id'
  )
  merged = mutate(merged, updated_value = ifelse(value < alpha & lfsr < 1 - alpha, 1, value))
  merged = left_join(merged,
    dplyr::select(guide_summaries$gene_shift, gene_shift = Mean, gene_id = gene_mapping))
  merged = rename(merged, old_value = value, value = updated_value)
  # merged = arrange(merged, value, desc(abs(gene_shift)))
  # merged = arrange(merged, value, desc(abs(lfsr)))
  merged = arrange(merged, value, abs(lfsr))
  merged
}

posterior_summary = gene_posterior_summary(samples_n, 0.10)

fbase = sub('.rds', '', basename(sample_fname))
saveRDS(posterior_summary, file = paste0('results/nimble/', fbase, '_posterior_summary.rds'))
