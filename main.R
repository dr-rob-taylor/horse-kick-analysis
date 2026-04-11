rm( list = ls() )

library( tidyverse )
library( gridExtra )
library( ggpubr )

source( "src/functions.R" )

# ------------------------------------------------------------------------------
# Data load and format
# ------------------------------------------------------------------------------

data_file <- list.files( "data" )
data_long <- read_csv( paste0("data/", data_file), col_types = cols() ) %>%
             pivot_longer( -Year, names_to = "Corp", values_to = "Deaths" )

# ------------------------------------------------------------------------------
# Fit Bayesian Poisson Model 
#
# Fits Poisson likelihood with conjugate Gamma prior on lambda. 
# Produces Gamma posterior distribution over lambda.
# Model fit with  Gamma(1, 1) prior
# ------------------------------------------------------------------------------

# Data 
y <- data_long %>% pull( Deaths )

# Prior parameters
alpha_prior <- 1
beta_prior  <- 1

# Model fit and summary
bayes_fit    <- fit_poisson_gamma( y = y, alpha = alpha_prior, beta = beta_prior )
bayes_summary( fit = bayes_fit )

# Generate outputs to examine posterior and fit
plots <- bayes_plot( fit = bayes_fit )
ggsave( plots, file = "figs/bayes_plots.png", units = "cm", height = 7, width = 19 )

# ------------------------------------------------------------------------------
# Maximum Likelihood Estimation 
# ------------------------------------------------------------------------------

lambda_mle <- data_long %>% pull( Deaths ) %>% mean()

data_summary <- data_long %>%
                count( Deaths ) %>%
                mutate( freq = n / sum(n), 
                        exp = dpois(0:4, lambda_mle) )
               
plot_mle <- data_summary %>%
            ggplot( aes(x = Deaths, y = freq)) +
            geom_point( aes(col = "Data") ) +
            geom_line( aes( y = exp)) +
            geom_point( aes( y = exp, col = "Fit") ) +
            labs( y = "Density", col = NULL, subtitle = "Poisson approximation of horse kick deaths") +
            theme_bw() +
            theme( legend.title = element_blank(), 
                   legend.position = "bottom") +
            scale_color_manual( values = c("firebrick", "black"))

ggsave( plot_mle, file = "figs/plot_mle.png", units = "cm", width = 10, height = 7)

