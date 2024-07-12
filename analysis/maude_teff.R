library('dplyr')
library('MAUDE')
library('tidyr')

source('../config.R', verbose = TRUE)

# first, need to aggregate counts
# for the Teff this has already been done by mageck

counts = data.table::fread(
  '../data/screen-2_FILTERED_Teff_D1_D2_D3_2019-04-23.count.txt',
  data.table = FALSE)

# get this into a key-like table so that we can manipulate it easily
g_counts = gather(counts, key = 'sample', value = 'count', -sgRNA, -Gene)
g_counts = dplyr::rename(g_counts, gid = sgRNA, element = Gene)
g_counts = mutate(g_counts, sample = sub('Teff_', '', sample))
gfp_counts = filter(g_counts, grepl('GFP', sample))
g_counts = filter(g_counts, !grepl('GFP', sample))

g_counts = separate(g_counts,
  sample, c('donor', 'reporter', 'bin'), sep = '_')

# use the GFP counts as the number of 'unsorted cells'
gfp_counts = mutate(gfp_counts, donor = sub('_GFP', '', sample))
gfp_counts = dplyr::select(gfp_counts, gid, donor, NS = count)

# maude requires a non-targeting logical
g_counts = mutate(g_counts, NT = grepl('non-targeting', element, ignore.case = TRUE))

# put the counts into columns consistent with the bin
s_counts = spread(g_counts, bin, count)

# put the GFP counts into the gated counts
s_counts = dplyr::select(s_counts, gid:NT, A = low, F = high)
s_counts = left_join(s_counts, gfp_counts, by = c('gid', 'donor'))

# this is a practice run only one donor and one sample
bin_bounds = makeBinModel(
  data.frame(Bin = c('A', 'F'), fraction = c(0.15, 0.15)))

tmp = dplyr::filter(s_counts, donor == 'D1' & reporter == 'CD25')
tmp = dplyr::select(tmp, -reporter)
tmp = dplyr::rename(tmp, screen = donor)
bin_bounds$screen = 'D1'
bin_bounds = bin_bounds[1:2, ]

guide_level = findGuideHitsAllScreens(unique(tmp['screen']), tmp, bin_bounds,
  sortBins = c('A', 'F'))
# end practice run

# now trying to analyze all of them together
all_counts = dplyr::mutate(s_counts, screen = paste0(reporter, '_', donor))
all_counts = dplyr::select(all_counts, -c(donor, reporter))

bin_bounds = makeBinModel(
  data.frame(Bin = c('A', 'F'), fraction = c(0.15, 0.15)))
bin_bounds = bin_bounds[1:2, ]

all_bounds = bind_rows(lapply(unique(all_counts$screen),
  function(s) {
    dplyr::mutate(bin_bounds, screen = s)
  }))

all_maude = findGuideHitsAllScreens(
  unique(all_counts['screen']),
  all_counts, all_bounds, sortBins = c('A', 'F'))

saveRDS(all_maude, file = '../results/Teff/guide_level_maude.rds')




guide_level = dplyr::mutate(guide_level, pvalue = 2 * pnorm(-abs(Z)))
hist(guide_level$pvalue, breaks = 100, xlim = c(0, 1))
hist(guide_level$Z)
plot(ecdf(abs(guide_level$Z)))

element_wise = getElementwiseStats(unique(tmp['screen']), tmp, guide_level)

install('~/forks/MAUDE')

debugonce(findGuideHits)

debugonce(getNBGaussianLikelihood)
guide_level = findGuideHits(tmp, bin_bounds, sortBins = c('A', 'F'), limits = c(-8, 8))


tmp$screen = 'test1'
tmp$gid = 1:nrow(tmp)

mapping = data.frame(element = unique(tmp$element))
mapping = mutate(mapping, id = 1:nrow(mapping))

tmp = inner_join(tmp, mapping, by = 'element')
tmp = select(tmp, -element)
tmp = dplyr::rename(tmp, element = id)

saveRDS(tmp, file = 'test1.rds')



bin_bounds = makeBinModel(
  data.frame(Bin = c('A', 'F'), fraction = c(0.15, 0.15)))
guide_level = findGuideHits(tmp, bin_bounds, sortBins = c('A', 'F'), limits = c(-8, 8))

# TODO: need to break up all of the replicates into their own screen
