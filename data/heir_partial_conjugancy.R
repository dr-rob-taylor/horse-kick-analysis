library(tidyverse)
library(patchwork)

n_corps <- length(deaths_per_corp)
n_years <- length(deaths_per_year)

# ------------------------------------------------------------------------------
# Log posterior for (alpha, beta) — marginalised over lambda_i
#
# Marginal likelihood: 
#.  k_i ~ NegBinom(alpha, beta / (beta + n))
# Priors: 
#   alpha ~ Exponential(1)
#.  beta ~ Gamma(2, 1)
# ------------------------------------------------------------------------------

log_posterior <- function(log_alpha, log_beta, counts, n) {
  alpha <- exp(log_alpha)
  beta  <- exp(log_beta)
  
  # Marginal log-likelihood (summed over corps)
  mll <- sum(dnbinom(counts,
                     size = alpha,
                     prob = beta / (beta + n),
                     log  = TRUE))
  
  # Log priors (on natural scale, with Jacobian for log transformation)
  log_prior_alpha <- dexp(alpha, rate = 1, log = TRUE) + log_alpha
  log_prior_beta  <- dgamma(beta, shape = 2, rate = 1, log = TRUE) + log_beta
  
  mll + log_prior_alpha + log_prior_beta
}

# ------------------------------------------------------------------------------
# Metropolis-Hastings sampler over (log_alpha, log_beta)
#
# Sampling on the log scale enforces positivity without hard constraints.
# Proposal: random walk Normal with tunable step sizes.
# ------------------------------------------------------------------------------

run_mh_sampler <- function(n_iter    = 50000,
                           n_warmup  = 10000,
                           step_alpha = 0.3,
                           step_beta  = 0.3,
                           counts    = deaths_per_corp,
                           n         = n_years,
                           seed      = 42) {
  set.seed(seed)
  
  # Initialise at empirical Bayes estimates as a warm start
  chain <- matrix(NA_real_, nrow = n_iter, ncol = 2,
                  dimnames = list(NULL, c("log_alpha", "log_beta")))
  
  chain[1, ] <- c(log(2), log(10))  # sensible starting values
  log_post_current <- log_posterior(chain[1, 1], chain[1, 2], counts, n)
  
  n_accept <- 0L
  
  for (i in seq(2, n_iter)) {
    # Propose
    proposal <- chain[i - 1, ] + rnorm(2, sd = c(step_alpha, step_beta))
    
    # Evaluate log posterior at proposal
    log_post_proposal <- log_posterior(proposal[1], proposal[2], counts, n)
    
    # Metropolis acceptance step
    log_r <- log_post_proposal - log_post_current
    
    if (log(runif(1)) < log_r) {
      chain[i, ]       <- proposal
      log_post_current <- log_post_proposal
      n_accept         <- n_accept + 1L
    } else {
      chain[i, ] <- chain[i - 1, ]
    }
  }
  
  # Discard warmup
  posterior_samples <- chain[(n_warmup + 1):n_iter, ]
  
  cat(sprintf("Acceptance rate: %.1f%%\n",
              100 * n_accept / n_iter))
  
  as_tibble(posterior_samples) %>%
    mutate(
      alpha = exp(log_alpha),
      beta  = exp(log_beta),
      prior_mean = alpha / beta
    )
}

samples <- run_mh_sampler()

# ------------------------------------------------------------------------------
# Sampler diagnostics
# ------------------------------------------------------------------------------

p_trace_alpha <- samples %>%
  mutate(iter = row_number()) %>%
  ggplot(aes(x = iter, y = alpha)) +
  geom_line(colour = "steelblue", linewidth = 0.3, alpha = 0.7) +
  labs(x = "Iteration", y = expression(alpha), title = "Trace: alpha") +
  theme_minimal(base_size = 12)

p_trace_beta <- samples %>%
  mutate(iter = row_number()) %>%
  ggplot(aes(x = iter, y = beta)) +
  geom_line(colour = "coral", linewidth = 0.3, alpha = 0.7) +
  labs(x = "Iteration", y = expression(beta), title = "Trace: beta") +
  theme_minimal(base_size = 12)

p_joint <- samples %>%
  ggplot(aes(x = alpha, y = beta)) +
  geom_density_2d_filled(alpha = 0.85, bins = 10) +
  scale_fill_viridis_d(option = "plasma", direction = -1) +
  labs(x = expression(alpha), y = expression(beta),
       title = "Joint posterior: alpha vs beta") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

p_trace_alpha + p_trace_beta + p_joint +
  plot_annotation(title = "MCMC diagnostics",
                  theme = theme(plot.title = element_text(size = 14,
                                                          face = "bold")))

ggsave("mcmc_diagnostics.png", width = 14, height = 4, dpi = 300)

# ------------------------------------------------------------------------------
# Corps-level posteriors: analytically reconstructed for each posterior draw
#
# For each sample (alpha^(s), beta^(s)):
#   lambda_i | k_i ~ Gamma(alpha^(s) + k_i, beta^(s) + n)
#
# Summarise by averaging posterior means and quantiles across draws.
# ------------------------------------------------------------------------------

n_samples <- nrow(samples)

corp_posteriors <- map_dfr(seq_len(n_corps), function(i) {
  
  alpha_post <- samples$alpha + deaths_per_corp[i]
  beta_post  <- samples$beta  + n_years
  
  # Posterior mean for lambda_i at each MCMC draw
  post_means <- alpha_post / beta_post
  
  # Draw one lambda sample per MCMC iteration for credible intervals
  lambda_draws <- rgamma(n_samples, shape = alpha_post, rate = beta_post)
  
  tibble(
    corp       = i,
    deaths     = deaths_per_corp[i],
    mle        = deaths_per_corp[i] / n_years,
    post_mean  = mean(post_means),
    post_lower = quantile(lambda_draws, 0.025),
    post_upper = quantile(lambda_draws, 0.975)
  )
})

print(corp_posteriors)

# ------------------------------------------------------------------------------
# Compare empirical Bayes vs fully Bayesian posterior means
# ------------------------------------------------------------------------------

# Rerun empirical Bayes for comparison
marginal_log_lik <- function(params, counts, n) {
  alpha <- exp(params[1])
  beta  <- exp(params[2])
  sum(dnbinom(counts, size = alpha,
              prob = beta / (beta + n), log = TRUE))
}

eb_fit    <- optim(c(log(2), log(10)), marginal_log_lik,
                   counts = deaths_per_corp, n = n_years,
                   control = list(fnscale = -1), method = "BFGS")
alpha_eb  <- exp(eb_fit$par[1])
beta_eb   <- exp(eb_fit$par[2])

eb_posteriors <- tibble(
  corp      = seq_len(n_corps),
  deaths    = deaths_per_corp,
  mle       = deaths_per_corp / n_years,
  post_mean = (alpha_eb + deaths_per_corp) / (beta_eb + n_years)
)

# ------------------------------------------------------------------------------
# Shrinkage comparison: Empirical Bayes vs partially conjugate Bayes
# ------------------------------------------------------------------------------

comparison_df <- corp_posteriors %>%
  select(corp, deaths, mle, post_mean_bayes = post_mean,
         post_lower, post_upper) %>%
  left_join(eb_posteriors %>% select(corp, post_mean_eb = post_mean),
            by = "corp")

p_compare <- comparison_df %>%
  ggplot(aes(y = reorder(factor(corp), deaths))) +
  geom_linerange(aes(xmin = post_lower, xmax = post_upper),
                 colour = "steelblue", linewidth = 1.0, alpha = 0.5) +
  geom_point(aes(x = post_mean_bayes, colour = "Partially conjugate Bayes"),
             size = 3) +
  geom_point(aes(x = post_mean_eb, colour = "Empirical Bayes"),
             size = 3, shape = 17) +
  geom_point(aes(x = mle, colour = "MLE"),
             size = 3, shape = 4, stroke = 1.2) +
  scale_colour_manual(
    name   = "Estimate",
    values = c("Partially conjugate Bayes" = "steelblue",
               "Empirical Bayes"           = "coral",
               "MLE"                       = "grey30")
  ) +
  labs(
    x        = expression(lambda[i]),
    y        = "Corps (ordered by deaths)",
    title    = "Posterior means: MLE vs Empirical Bayes vs Partially Conjugate Bayes",
    subtitle = "Bar = 95% credible interval (partially conjugate Bayes)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())

p_compare
ggsave("comparison_plot.png", width = 11, height = 6, dpi = 300)