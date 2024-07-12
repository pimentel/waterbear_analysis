library('PreciseSums')
library('ggplot2')
library('cowplot')
library('dplyr')

theme_set(theme_cowplot(25))

lpxgt2 = function(n, q, lambda) {
  lpx2 = log1p(-((1 - q)^n + n * q * (1 - q)^(n - 1)))
  lpn = dpois(n, lambda, log = TRUE)
  lpx2 + lpn
}

tmp = sapply(2:1e6, lpxgt2, q = 0.10, lambda = 0.3)
kahanSum(exp(tmp))

moi = seq(0.25, 7, by = 0.25)

p_effect = 0.05
p_05 = lapply(moi,
  function(m) {
    print(m)
    lp = sapply(2:1e6, lpxgt2, q = p_effect, lambda = m)
    data.frame(moi = m, q = p_effect, p_gt_2 = kahanSum(exp(lp)), moi = m)
  })
p_05 = bind_rows(p_05)

p_effect = 0.10
p_10 = lapply(moi,
  function(m) {
    print(m)
    lp = sapply(2:1e6, lpxgt2, q = p_effect, lambda = m)
    data.frame(moi = m, q = p_effect, p_gt_2 = kahanSum(exp(lp)), moi = m)
  })
p_10 = bind_rows(p_10)

p_effect = 0.25
p_25 = lapply(moi,
  function(m) {
    print(m)
    lp = sapply(2:1e6, lpxgt2, q = p_effect, lambda = m)
    data.frame(moi = m, q = p_effect, p_gt_2 = kahanSum(exp(lp)), moi = m)
  })
p_25 = bind_rows(p_25)

p_effect = 0.50
p_50 = lapply(moi,
  function(m) {
    print(m)
    lp = sapply(2:1e6, lpxgt2, q = p_effect, lambda = m)
    data.frame(moi = m, q = p_effect, p_gt_2 = kahanSum(exp(lp)), moi = m)
  })
p_50 = bind_rows(p_50)


pgt = bind_rows(p_05, p_10, p_25, p_50)
pgt = mutate(pgt, pr_effect = factor(q))

cbb <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

theme_set(theme_cowplot(25))

p = ggplot(filter(pgt, pr_effect != 0.5), aes(moi, p_gt_2, group = pr_effect, color = pr_effect))
# p = p + geom_point()
p = p + geom_line(size = 4)
p = p + coord_cartesian(ylim = c(0, 1))
p = p + xlab('MOI')
# p = p + ylab('Pr(# guides with effects in a cell >= 2)')
p = p + ylab('Proportion of cells with two or more effects')
p = p + scale_color_manual(values = cbb, name = 'Proportion of\ngenes with\nan effect')
p = p + theme(legend.position = c(0.2, 0.8))
save_plot('pr_collisions.pdf', p, base_width = 12, base_height = 8)

# talk

lpxgt2 = function(n, q, lambda) {
  lpx2 = log1p(-((1 - q)^n + n * q * (1 - q)^(n - 1)))
  lpn = dpois(n, lambda, log = TRUE)
  lpx2 + lpn
}

tmp = sapply(2:1e6, lpxgt2, q = 0.10, lambda = 0.3)
kahanSum(exp(tmp))

moi = seq(0.25, 7, by = 0.25)

p_effect = 0.05
p_05 = lapply(moi,
  function(m) {
    print(m)
    lp = sapply(2:1e6, lpxgt2, q = p_effect, lambda = m)
    data.frame(moi = m, q = p_effect, p_gt_2 = kahanSum(exp(lp)), moi = m)
  })
p_05 = bind_rows(p_05)

p_effect = 0.10
p_10 = lapply(moi,
  function(m) {
    print(m)
    lp = sapply(2:1e6, lpxgt2, q = p_effect, lambda = m)
    data.frame(moi = m, q = p_effect, p_gt_2 = kahanSum(exp(lp)), moi = m)
  })
p_10 = bind_rows(p_10)

p_effect = 0.25
p_25 = lapply(moi,
  function(m) {
    print(m)
    lp = sapply(2:1e6, lpxgt2, q = p_effect, lambda = m)
    data.frame(moi = m, q = p_effect, p_gt_2 = kahanSum(exp(lp)), moi = m)
  })
p_25 = bind_rows(p_25)

p_effect = 0.50
p_50 = lapply(moi,
  function(m) {
    print(m)
    lp = sapply(2:1e6, lpxgt2, q = p_effect, lambda = m)
    data.frame(moi = m, q = p_effect, p_gt_2 = kahanSum(exp(lp)), moi = m)
  })
p_50 = bind_rows(p_50)


pgt = bind_rows(p_05, p_10, p_25, p_50)
# pgt = bind_rows(p_05)
pgt = mutate(pgt, pr_effect = factor(q))

cbb <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

theme_set(theme_cowplot(25))

p = ggplot(pgt, aes(moi, p_gt_2, group = pr_effect, color = pr_effect))
# p = p + geom_point()
p = p + geom_line(size = 4)
p = p + coord_cartesian(ylim = c(0, 0.5))
p = p + xlab('MOI')
p = p + ylab('Pr(# guides with effects in a cell >= 2)')
# p = p + ylab('')
p = p + scale_color_manual(values = cbb, name = 'Proportion of\nguides\nwith effects')
p = p + theme(legend.position = c(0.05, 0.75))
save_plot('pr_collisions_big_text.pdf', p, base_width = 12, base_height = 8)



rnorm(0.7,
