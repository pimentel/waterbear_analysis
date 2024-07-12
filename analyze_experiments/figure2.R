library('dplyr')
library('ggplot2')
library('tidyr')
library('cowplot')

theme_set(theme_cowplot())

# standard normal discretized
effect_size = 0.8
x = seq(-3, 3, length.out = 1000)
y_null = dnorm(x, 0, 1)
y_effect = dnorm(x, effect_size, 1)
df = data.frame(x, y_null, y_effect)
bin_cutoffs = qnorm(c(0.25, 0.5, 0.75))

bin_labels = data.frame(
  bin_text = c('bin: ', 1:4),
  x = c(-2.5,
    -2,
    bin_cutoffs[1] + (bin_cutoffs[2] - bin_cutoffs[1]) / 2,
    qnorm(0.75) / 2,
    2),
  y = 0.39
)


p = ggplot(df, aes(x, y_null))
p = p + geom_area(fill = 'lightgray')
p = p + geom_vline(xintercept = bin_cutoffs, linetype = 2, color = '#bdbdbd')
p = p + geom_text(aes(x, y, label = bin_text), data = bin_labels)
p = p + xlab('expected bin density')
p = p + ylab(element_blank())
p = p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
  axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.75))
p = p + scale_x_continuous(breaks = bin_labels$x[-1], labels = c(0.25, 0.25, 0.25, 0.25))
save_plot(p, file = 'img/figure2_null.pdf', base_height = 2, base_width = 3)


effect_density = round(
  pnorm(c(bin_cutoffs, Inf), effect_size) - pnorm(c(-Inf, bin_cutoffs), effect_size),
  2)
p = ggplot(df, aes(x, y_effect))
p = p + geom_area(fill = '#1A85FF')
p = p + geom_vline(xintercept = bin_cutoffs, linetype = 2, color = 'gray')
p = p + geom_text(aes(x, y, label = bin_text), data = bin_labels)
p = p + xlab('expected bin density')
p = p + ylab(element_blank())
p = p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
  axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0.75))
p = p + scale_x_continuous(breaks = bin_labels$x[-1], labels = effect_density)
save_plot(p, file = 'img/figure2_effect.pdf', base_height = 2, base_width = 3)
