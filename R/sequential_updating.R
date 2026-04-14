library(tidyverse)

# ------------------------------------------------------------------------------
# Horse kick data: deaths per year, summed across 14 cavalry corps
# Source: von Bortkiewicz (1898)
# ------------------------------------------------------------------------------

data_file <- list.files( "data" )
data <- read_csv( paste0("data/", data_file), col_types = cols() ) |> 
        select(-Year)

deaths_per_year <- rowSums(data)

years <- 1875:1894
n_corps <- 14  # observations per year

# ------------------------------------------------------------------------------
# Sequential Bayesian updating
# Prior: Gamma(alpha_0, beta_0)
# Posterior after year t: Gamma(alpha + sum(k), beta + n)
# ------------------------------------------------------------------------------

alpha_0 <- 1
beta_0  <- 1

n_years <- length(years)

# Store posterior parameters and summaries after each year
params <- tibble(
  year       = years,
  alpha_post = cumsum(deaths_per_year) + alpha_0,
  beta_post  = (seq_len(n_years) * n_corps) + beta_0,
  post_mean  = alpha_post / beta_post,
  post_lower = qgamma(0.025, alpha_post, beta_post),
  post_upper = qgamma(0.975, alpha_post, beta_post)
)

# ------------------------------------------------------------------------------
# Panel 1: Evolving posterior densities
# ------------------------------------------------------------------------------

lambda_grid <- seq(0.3, 1.2, length.out = 500)

density_df <- params %>%
  rowwise() %>%
  reframe(
    year   = year,
    lambda = lambda_grid,
    density = dgamma(lambda_grid, shape = alpha_post, rate = beta_post)
  )

# Thin to every 4 years for readability, always include first and last
years_to_plot <- unique(c(
  years[seq(1, n_years, by = 4)],
  years[n_years]
))

p1 <- density_df %>%
  filter(year %in% years_to_plot) %>%
  mutate(year = factor(year)) %>%
  ggplot(aes(x = lambda, y = density, colour = year)) +
  geom_line(linewidth = 0.8) +
  scale_colour_viridis_d(option = "plasma", direction = -1, name = "Year") +
  labs(
    x = expression(lambda),
    y = "Density",
    title = "Posterior density over time"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")

# ------------------------------------------------------------------------------
# Panel 2: Posterior mean and 95% credible interval over time
# ------------------------------------------------------------------------------

p2 <- params %>%
  ggplot(aes(x = year)) +
  geom_ribbon(aes(ymin = post_lower, ymax = post_upper),
              fill = "steelblue", alpha = 0.2) +
  geom_line(aes(y = post_mean), colour = "steelblue", linewidth = 0.9) +
  geom_point(aes(y = post_mean), colour = "steelblue", size = 2) +
  geom_hline(yintercept = 196 / 280, linetype = "dashed",
             colour = "grey40", linewidth = 0.6) +
  annotate("text", x = 1876, y = 0.73,
           label = "Batch MLE (0.70)", colour = "grey40", size = 3.5, hjust = 0) +
  labs(
    x = "Year",
    y = expression(lambda),
    title = "Posterior mean and 95% credible interval"
  ) +
  scale_x_continuous(breaks = seq(1875, 1894, by = 5)) +
  theme_minimal(base_size = 13)

# ------------------------------------------------------------------------------
# Combine and save
# ------------------------------------------------------------------------------

library(patchwork)

p1 + p2 +
  plot_annotation(
    title   = "Sequential Bayesian updating: Prussian horse kick data",
    theme   = theme(plot.title = element_text(size = 15, face = "bold"))
  )

ggsave("sequential_updating.png", width = 12, height = 5, dpi = 300)
