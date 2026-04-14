rm( list = ls() )

library( tidyverse )
library( gridExtra )
#library( ggpubr )
library(patchwork)

source( "src/functions.R" )
load("data/prussian_data.RData")

# ------------------------------------------------------------------------------
# Horse kick data: deaths per year, summed across 14 cavalry corps
# Source: von Bortkiewicz (1898)
# ------------------------------------------------------------------------------

# data_path <- list.files( pattern = ".csv", recursive = T)
# data_raw  <- read.csv( data_path, row.names = 1)
# data_long <- data_raw |> 
#              pivot_longer( everything(), names_to = "Corp", values_to = "Deaths" )
# 
# deaths_per_year <- rowSums( data_raw )
# deaths_per_corp <- colSums( data_raw )

# ------------------------------------------------------------------------------
# Fit Bayesian Poisson Model 
#
# Prior: Gamma(alpha_prior, beta_prior)
# Posterior: Gamma(alpha_prior + sum(k), beta_prior + n)
#
# Change alpha_prior and beta_prior to see effect on posterior
# ------------------------------------------------------------------------------

# Data 
k <- pull( data_long[,"Deaths"] )

# Default prior parameters (change)
alpha_prior <- 1
beta_prior  <- 1

# Model fit and summary
bayes_fit <- fit_poisson_gamma( y = k, alpha = alpha_prior, beta = beta_prior )
bayes_summary( fit = bayes_fit )

# Generate outputs to examine posterior and fit
plots <- bayes_plot( fit = bayes_fit )
ggsave( plots, file = "figs/bayes_plots.png", units = "cm", height = 7, width = 19 )


if(F){
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

}