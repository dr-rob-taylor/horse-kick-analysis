library(tidyverse)

n_corps <- length(deaths_per_corp)
n_years <- length(deaths_per_year)

# ------------------------------------------------------------------------------
# Step 1: Estimate hyperparameters via marginal maximum likelihood (Empirical Bayes)
#
# Marginalising out lambda_i, the marginal likelihood for each corps is:
#   k_i ~ NegBinom(alpha, beta / (beta + n_years))
#
# We maximise the marginal log-likelihood over (alpha, beta) numerically.
# Note: beta here is the rate parameter of the Gamma hyperprior.
# ------------------------------------------------------------------------------

marginal_log_lik <- function(params, counts, n) {
  alpha <- exp(params[1])  # exponentiate to enforce positivity
  beta  <- exp(params[2])
  
  sum(dnbinom(counts,
              size = alpha,
              prob = beta / (beta + n),
              log  = TRUE))
}

# Optimise over log(alpha) and log(beta)
fit <- optim(
  par     = c(log(2), log(10)),  # sensible starting values
  fn      = marginal_log_lik,
  counts  = deaths_per_corp,
  n       = n_years,
  control = list(fnscale = -1),  # maximise
  method  = "BFGS"
)

alpha_hat <- exp(fit$par[1])
beta_hat  <- exp(fit$par[2])

cat(sprintf("Empirical Bayes estimates:\n  alpha = %.3f\n  beta  = %.3f\n",
            alpha_hat, beta_hat))
cat(sprintf("  Implied prior mean: %.3f\n", alpha_hat / beta_hat))
cat(sprintf("  Implied prior variance: %.3f\n", alpha_hat / beta_hat^2))

# ------------------------------------------------------------------------------
# Step 2: Compute closed-form Gamma posterior for each corps
#
# Given alpha_hat and beta_hat, the posterior for lambda_i is:
#   lambda_i | k_i ~ Gamma(alpha_hat + k_i, beta_hat + n_years)
# ------------------------------------------------------------------------------

posteriors <- tibble(
  corp        = factor(seq_len(n_corps)),
  corp_label  = names(deaths_per_corp),
  deaths      = deaths_per_corp,
  mle         = deaths_per_corp / n_years,
  alpha_post  = alpha_hat + deaths_per_corp,
  beta_post   = beta_hat  + n_years,
  post_mean   = alpha_post / beta_post,
  post_lower  = qgamma(0.025, alpha_post, beta_post),
  post_upper  = qgamma(0.975, alpha_post, beta_post)
)

print(posteriors %>% select(corp, corp_name, deaths, mle, post_mean, post_lower, post_upper),
      n = n_corps)

hypers <- tibble(
  alpha_post = alpha_hat,
  beta_post = beta_hat,
  post_mean = alpha_hat / beta_hat,
  post_lower  = qgamma(0.025, alpha_hat, beta_hat),
  post_upper  = qgamma(0.975, alpha_hat, beta_hat)
)

# ------------------------------------------------------------------------------
# Step 3: Shrinkage plot — MLE vs posterior mean per corps
# ------------------------------------------------------------------------------

grand_mean <- sum(deaths_per_corp) / (n_corps * n_years)

p_shrinkage <- posteriors %>%
  ggplot(aes(x = mle, y = post_mean)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey60", linewidth = 0.6) +
  geom_hline(yintercept = alpha_hat / beta_hat,
             linetype = "dotted", colour = "steelblue", linewidth = 0.6) +
  geom_segment(aes(xend = mle, yend = mle),
               colour = "grey70", linewidth = 0.4) +
  geom_point(size = 3, colour = "steelblue") +
  annotate("text", x = 0.42, y = alpha_hat / beta_hat + 0.015,
           label = "Prior mean", colour = "steelblue", size = 3.5, hjust = 0) +
  labs(
    x     = "MLE (deaths / 20 years)",
    y     = "Posterior mean",
    subtitle = "Shrinkage plot"
  ) +
  theme_minimal(base_size = 13)

# ------------------------------------------------------------------------------
# Step 4: Posterior densities per corps
# ------------------------------------------------------------------------------

lambda_grid <- seq(0.1, 1.8, length.out = 500)

density_df <- posteriors %>%
  rowwise() %>%
  reframe(
    corp    = corp,
    deaths  = deaths,
    lambda  = lambda_grid,
    density = dgamma(lambda_grid, shape = alpha_post, rate = beta_post)
  )

p_densities <- density_df %>%
  ggplot(aes(x = lambda, y = density,
             colour = reorder(corp, deaths),
             group = corp)) +
  geom_line(linewidth = 0.7, alpha = 0.85) +
  scale_colour_viridis_d(option = "plasma", name = "Corps\n(by deaths)") +
  labs(
    x     = expression(lambda[i]),
    y     = "Density",
    subtitle = "Corps-level posterior distributions"
  ) +
  theme_minimal(base_size = 13)

# ------------------------------------------------------------------------------
# Step 5: Combine and save
# ------------------------------------------------------------------------------

library(patchwork)

p_shrinkage + p_densities +
  plot_annotation(
    title = "Hierarchical Gamma-Poisson model: Prussian horse kick data",
    theme = theme(plot.title = element_text(size = 15, face = "bold"))
  )

ggsave("hierarchical_model.png", width = 13, height = 5, dpi = 300)

library(ggridges)

# ------------------------------------------------------------------------------
# Ridgeline plot: posterior distributions per corps, ordered by deaths
# ------------------------------------------------------------------------------

# Sample from each corps' posterior for the ridgeline
set.seed(42)
n_samples <- 5000

ridgeline_df <- posteriors %>%
  rowwise() %>%
  reframe(
    corp        = corp,
    deaths      = deaths,
    corp_label  = paste0(corp_label, " (", deaths, " deaths)"),
    lambda      = rgamma(n_samples, shape = alpha_post, rate = beta_post)
  ) %>%
  mutate(corp_label = reorder(corp_label, deaths))

p_ridgeline <- ridgeline_df %>%
  ggplot(aes(x = lambda, y = corp_label, fill = deaths)) +
  geom_density_ridges(
    scale         = 1.8,
    rel_min_height = 0.01,
    colour        = "white",
    linewidth     = 0.4,
    alpha         = 0.85
  ) +
  geom_vline(xintercept = alpha_hat / beta_hat,
             linetype = "dashed", colour = "grey30", linewidth = 0.6) +
  annotate("text", x = alpha_hat / beta_hat + 0.02, y = 1.5,
           label = "Prior mean", colour = "grey30", size = 3.2, hjust = 0) +
  scale_fill_viridis_c(option = "plasma", name = "Deaths", direction = -1) +
  scale_x_continuous(limits = c(0, 1.8), breaks = seq(0, 1.8, by = 0.3)) +
  labs(
    x     = expression(lambda[i]),
    y     = NULL,
    subtitle = "Corps-level posterior distributions",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y     = element_text(size = 9),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------------------------------------
# Combined figure
# ------------------------------------------------------------------------------

p_shrinkage + p_ridgeline +
  plot_annotation(
    title = "Hierarchical Gamma-Poisson model: Prussian horse kick data",
    theme = theme(plot.title = element_text(size = 15, face = "bold"))
  )

ggsave("hierarchical_model.png", width = 14, height = 6, dpi = 300)

# ------------------------------------------------------------------------------
# Posterior predictive check: observed vs. predicted deaths per corps
# ------------------------------------------------------------------------------

set.seed(42)
n_samples <- 10000

ppc_df <- posteriors %>%
  rowwise() %>%
  reframe(
    corp       = corp,
    deaths     = deaths,
    corp_label = paste0(corp_label, " (", deaths, " deaths)"),
    # Posterior predictive: integrate out lambda_i
    # NegBinom with size = alpha_post, prob = beta_post / (beta_post + n_years)
    pred       = rnbinom(n_samples,
                         size = alpha_post,
                         prob = beta_post / (beta_post + n_years))
  ) %>%
  group_by(corp, deaths, corp_label) %>%
  summarise(
    pred_mean  = mean(pred),
    pred_lower = quantile(pred, 0.025),
    pred_upper = quantile(pred, 0.975),
    .groups    = "drop"
  ) %>%
  mutate(corp_label = reorder(corp_label, deaths))

p_ppc <- ppc_df %>%
  ggplot(aes(y = corp_label)) +
  geom_linerange(
    aes(xmin = pred_lower, xmax = pred_upper, colour = deaths),
    linewidth = 1.2, alpha = 0.7
  ) +
  geom_point(
    aes(x = pred_mean, colour = deaths),
    size = 3
  ) +
  geom_point(
    aes(x = deaths),
    shape = 4, size = 3, stroke = 1.2, colour = "grey20"
  ) +
  geom_vline(xintercept = mean(deaths_per_corp),
             linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  scale_colour_viridis_c(option = "plasma", name = "Deaths", direction = -1) +
  scale_x_continuous(breaks = seq(0, 35, by = 5)) +
  labs(
    x     = "Total deaths (20 years)",
    y     = NULL,
    subtitle = "Posterior predictive check",
    caption = "Dot = posterior predictive mean | Bar = 95% interval | X = observed"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y      = element_text(size = 9),
    legend.position  = "none",
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------------------------------------
# Three-panel combined figure
# ------------------------------------------------------------------------------

p_shrinkage + p_ridgeline + p_ppc +
  plot_annotation(
    title = "Hierarchical Gamma-Poisson model: Empirical Bayes",
    theme = theme(plot.title = element_text(size = 15, face = "bold"))
  ) +
  plot_layout(widths = c(1, 1.2, 1))

ggsave("figs/eb_hierarchical_model.png", width = 18, height = 6, dpi = 300)
