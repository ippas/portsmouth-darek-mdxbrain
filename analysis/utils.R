library(gplots)
library(devtools)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")


gg_qqplot <- function(ps, ci = 0.95) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), shape = 16, size = 1) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
    geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po)
}

inflation <- function(ps) {
  chisq <- qchisq(1 - ps, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  lambda
}

filter_fold <- function(df, sample_info, real_fold = 2) {
  sample_info_renamed <-
    sample_info %>%
    rename(
      sample_id = id,
      sample_genotype = genotype,
      sample_age_days = age_days,
      sample_tissue = tissue
    )
  df %>%
    pivot_longer(
      matches('sample_[MB][DW][1-5]C?'),
      names_to = 'sample_id',
      values_to = 'value'
    ) %>%
    left_join(sample_info_renamed, by = 'sample_id') %>%
    group_by(gene_stable_id, sample_genotype, sample_age_days, sample_tissue) %>%
    mutate(m = mean(value)) %>%
    group_by(gene_stable_id) %>%
    mutate(
      min_m = min(m),
      max_m = max(m)
    ) %>%
    ungroup %>%
    filter(max_m - min_m > log2(real_fold)) %>%
    select(
      -mouse_id, -sample_genotype, -sample_age_days, -sample_tissue,
      -m, -min_m, -max_m
    ) %>%
    pivot_wider(names_from = 'sample_id', values_from = 'value')
}

pass <- function(x) {
  x
}
