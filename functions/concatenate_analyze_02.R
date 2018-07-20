ml3_dat <- read.csv('../data/study_02/ml3_stroop.csv')
    # Compute the results for the ML3 data
df_ml3 <- plyr::ddply(ml3_dat, .(study_name), function (x) {
  res_t <- t.test(x$StroopEffect, mu = 0, var.equal = TRUE)
  es_r <- sqrt((res_t$statistic / res_t$parameter) / ((res_t$statistic / res_t$parameter) + 1))
  res <- data.frame(type = 'ml3',
    benford_congr_m_p = digit_analysis(x$MC, type = 'benford')$pval,
    benford_congr_sd_p = digit_analysis(x$SDC, type = 'benford')$pval,
    benford_incongr_m_p = digit_analysis(x$MI, type = 'benford')$pval,
    benford_incongr_sd_p = digit_analysis(x$SDI, type = 'benford')$pval,
    terminal_congr_m_p = digit_analysis(x$MC, type = 'terminal')$pval,
    terminal_congr_sd_p = digit_analysis(x$SDC, type = 'terminal')$pval,
    terminal_incongr_m_p = digit_analysis(x$MI, type = 'terminal')$pval,
    terminal_incongr_sd_p = digit_analysis(x$SDI, type = 'terminal')$pval,
    std_congr_p = std_var(sds = x$SDC, n = x$NC, iter = iter, method = 'maxmin', subgroups = rep(0, length(x$SDC))),
    std_incongr_p = std_var(sds = x$SDI, n = x$NI, iter = iter, method = 'maxmin', subgroups = rep(0, length(x$SDI))),
    mult_m_sd_congr = cor(x$MC, x$SDC),
    mult_m_sd_incongr = cor(x$MI, x$SDI),
    mult_m_m_across = cor(x$MC, x$MI),
    mult_sd_sd_across = cor(x$SDC, x$SDI),
    es_r,
    n = length(x$MC))
  return(res)
})

names(df_ml3)[1] <- 'id'
responses <- list.files('../data/study_02/responses')
df_fab <- NULL

for (response in responses) {
  fab_dat <- read.table(sprintf('../data/study_02/responses/%s', response),
   sep = '\t', header = TRUE)
        # For clarity sake
  names(fab_dat) <- c('id',
    'mean_congruent',
    'sd_congruent',
    'congruent_trials',
    'mean_incongruent',
    'sd_incongruent',
    'incongruent_trials')
        # Eliminate overprecision in reporting
        # Instructed participants to fabricate milliseconds and 123.45 is result of copy-paste
        # Spreadsheet automatically removed this, but by copy-pasting this was undone by PP
  fab_dat <- as.data.frame(apply(fab_dat, 2, function(x) round(x, 0)))
        # Stroop effect size calculation
  sdpooled <- sqrt(((fab_dat$congruent_trials - 1) * fab_dat$sd_congruent^2 + (fab_dat$incongruent_trials - 1) * fab_dat$sd_incongruent^2) / ((fab_dat$congruent_trials + fab_dat$incongruent_trials - 2)))
  stroopeffect <- (fab_dat$mean_incongruent - fab_dat$mean_congruent) / sdpooled
  x <- t.test(stroopeffect, mu = 0, var.equal = TRUE)
  es_r <- sqrt((x$statistic / x$parameter) / ((x$statistic / x$parameter) + 1))
  
  binder <- data.frame(id = response,
   type = 'fabricated',
   benford_congr_m_p = digit_analysis(fab_dat$mean_congruent, type = 'benford')$pval,
   benford_congr_sd_p = digit_analysis(fab_dat$sd_congruent, type = 'benford')$pval,
   benford_incongr_m_p = digit_analysis(fab_dat$mean_incongruent, type = 'benford')$pval,
   benford_incongr_sd_p = digit_analysis(fab_dat$sd_incongruent, type = 'benford')$pval,
   terminal_congr_m_p = digit_analysis(fab_dat$mean_congruent, type = 'terminal')$pval,
   terminal_congr_sd_p = digit_analysis(fab_dat$sd_congruent, type = 'terminal')$pval,
   terminal_incongr_m_p = digit_analysis(fab_dat$mean_incongruent, type = 'terminal')$pval,
   terminal_incongr_sd_p = digit_analysis(fab_dat$sd_incongruent, type = 'terminal')$pval,
   std_congr_p = std_var(sds = fab_dat$sd_congruent, n = fab_dat$congruent_trials, iter = iter, method = 'maxmin', subgroups = rep(0, length(fab_dat$congruent_trials))),
   std_incongr_p = std_var(sds = fab_dat$sd_incongruent, n = fab_dat$incongruent_trials, iter = iter, method = 'maxmin', subgroups = rep(0, length(fab_dat$incongruent_trials))),
   mult_m_sd_congr = cor(fab_dat$mean_congruent, fab_dat$sd_congruent),
   mult_m_sd_incongr = cor(fab_dat$mean_incongruent, fab_dat$sd_incongruent),
   mult_m_m_across = cor(fab_dat$mean_congruent, fab_dat$mean_incongruent),
   mult_sd_sd_across = cor(fab_dat$sd_congruent, fab_dat$sd_incongruent),
   es_r = es_r,
   n = length(fab_dat$congruent_trials))

  df_fab <- rbind(df_fab, binder)
}

# Add handcoded Random number generator use
df_ml3$rng <- NA
df_fab$rng <- as.factor(c(1,
  1,
  1,
  0,
  1,
  0,
  0,
  0,
  1,
  0,
  1,
  1,
  1,
  1,
  1,
  0,
  1,
  0,
  1,
  1,
  1,
  1,
  0,
  1,
  1,
  1,
  1,
  0))

df <- rbind(df_ml3, df_fab)

    # Meta-analyze the multivariate associations to acquire the parametric estimate
df$fisher_mult_m_sd_congr = atanh(df$mult_m_sd_congr)
df$fisher_mult_m_sd_incongr = atanh(df$mult_m_sd_incongr)
df$fisher_mult_m_m_across = atanh(df$mult_m_m_across)
df$fisher_mult_sd_sd_across = atanh(df$mult_sd_sd_across)

x <- metafor::rma(df$fisher_mult_m_sd_congr, sei = 1 / sqrt(df$n - 3), method = 'REML', data = df)
fisher_mult_m_sd_congr_mean <- x$b[1]
fisher_mult_m_sd_congr_std <- sqrt(x$tau2[1])

x <- metafor::rma(df$fisher_mult_m_sd_incongr, sei = 1 / sqrt(df$n - 3), method = 'REML', data = df)
fisher_mult_m_sd_incongr_mean <- x$b[1]
fisher_mult_m_sd_incongr_std <- sqrt(x$tau2[1])

x <- metafor::rma(df$fisher_mult_m_m_across, sei = 1 / sqrt(df$n - 3), method = 'REML', data = df)
fisher_mult_m_m_across_mean <- x$b[1]
fisher_mult_m_m_across_std <- sqrt(x$tau2[1])

x <- metafor::rma(df$fisher_mult_sd_sd_across, sei = 1 / sqrt(df$n - 3), method = 'REML', data = df)
fisher_mult_sd_sd_across_mean <- x$b[1]
fisher_mult_sd_sd_across_std <- sqrt(x$tau2[1])

    # Compute the multivariate p-values for the ML data
z_val <- (df$fisher_mult_m_sd_congr - fisher_mult_m_sd_congr_mean) / fisher_mult_m_sd_congr_std
df$p_fisher_mult_m_sd_congr <- 2 * pnorm(q = abs(z_val), mean = 0, sd = 1, lower = FALSE)

z_val <- (df$fisher_mult_m_sd_incongr - fisher_mult_m_sd_incongr_mean) / fisher_mult_m_sd_incongr_std
df$p_fisher_mult_m_sd_incongr <- 2 * pnorm(q = abs(z_val), mean = 0, sd = 1, lower = FALSE)

z_val <- (df$fisher_mult_m_m_across - fisher_mult_m_m_across_mean) / fisher_mult_m_m_across_std
df$p_fisher_mult_m_m_across <- 2 * pnorm(q = abs(z_val), mean = 0, sd = 1, lower = FALSE)

z_val <- (df$fisher_mult_sd_sd_across - fisher_mult_sd_sd_across_mean) / fisher_mult_sd_sd_across_std
df$p_fisher_mult_sd_sd_across <- 2 * pnorm(q = abs(z_val), mean = 0, sd = 1, lower = FALSE)

df$fisher_combination <- apply(cbind(df$terminal_congr_m_p,
  df$terminal_congr_sd_p,
  df$terminal_incongr_m_p,
  df$terminal_incongr_sd_p,
  df$std_congr_p,
  df$std_incongr_p,
  df$p_fisher_mult_m_sd_congr,
  df$p_fisher_mult_m_sd_incongr,
  df$p_fisher_mult_m_m_across,
  df$p_fisher_mult_sd_sd_across), 1, function(x) fisher_method(x)[,3])

write.csv(df, file = '../data/study_02/ml3_fabricated_processed_collated.csv', row.names = FALSE)