fit_poisson_gamma <- function(y, alpha, beta)
{
  
  n <- length( y )
  
  lambda_lower <- floor( qgamma( .0001, alpha, beta ) )
  lambda_upper <- ceiling( qgamma( .9999, alpha, beta ) )
  lambda_seq   <- seq(lambda_lower, lambda_upper, by = .01 )
  
  alpha_post <- alpha + sum(y)
  beta_post  <- beta + n 
  
  prior     <- dgamma( lambda_seq, alpha, beta )
  posterior <- dgamma( lambda_seq, alpha_post, beta_post )
  
  k_lower <- qpois( .0001, mean(y) )
  k_upper <- qpois( .9999, mean(y) )
  k_seq <- k_lower:k_upper
  
  posterior_pred <- dnbinom(k_seq, size = alpha_post, prob = beta_post / (1 + beta_post) )
  
  cred_int <- c( lower = qgamma(.025, alpha_post, beta_post),
                 upper = qgamma(.975, alpha_post, beta_post ) )
  
  posterior_mode   <- (alpha_post - 1) / beta_post
  posterior_median <- qgamma(.50, alpha_post, beta_post)
  posterior_mean   <- alpha_post / beta_post
  posterior_var    <- alpha_post / beta_post ^ 2
  posterior_sd     <- sqrt( posterior_var )
  
  output <- list()
  output$data <- y
  output$n_obs <- n
  output$prior <- paste0("Gamma(", alpha, ", ", beta, ")")
  output$prior_density <- data.frame( lambda = lambda_seq, density = prior )
  output$posterior <- paste0("Gamma(", alpha_post, ", ", beta_post, ")")
  output$posterior_pars <- c( alpha_post = alpha_post, beta_post = beta_post )
  output$posterior_density <- data.frame( lambda = lambda_seq, density = posterior )
  output$posterior_pred <- data.frame( k = k_seq, density = posterior_pred )
  output$posterior_summary <- c( mode   = posterior_mode, 
                                 median = posterior_median, 
                                 mean   = posterior_mean, 
                                 variance = posterior_var, 
                                 std_dev = posterior_sd )
  output$posterior_ci <- cred_int
  
  return( output )
  
}

sample_posterior <- function( fit, n_sim = 1000 )
{
  
  n <- fit$n_obs
  
  sim_output <- array(NA, dim = c(n_sim * n, 2))
  
  for( ii in 1:n_sim ){
    
    index <- ((n * (ii - 1)) + 1):(n * ii)
    
    sim_output[index,1] <- ii
    sim_output[index,2] <- rnbinom(n, size = fit$posterior_pars["alpha_post"], 
                                   prob = fit$posterior_pars["beta_post"] / (1 + fit$posterior_pars["beta_post"]) )
    
  }
  colnames( sim_output ) <- c("index", "sample")
  output <- as_tibble( sim_output )
  
  tidy_out <- output %>% 
    group_by( index ) %>% 
    count( sample ) %>% 
    mutate(freq = n / sum(n) ) %>% 
    ungroup() %>%
    complete( index, sample, fill = list(n = 0, freq = 0) )
  
  return( tidy_out )
  
}

bayes_plot <- function( fit )
{
  
  ## Compare prior and posterior distribution
  
  p1 <- fit$prior_density %>%
    ggplot( aes(x = lambda) ) +
    geom_line( aes(y = density, col = "Prior") ) +
    geom_line( data = fit$posterior_density, aes(y = density, col = "Posterior")) +
    labs( x = "Lambda", y = "Density", subtitle = "Prior and posterior density") +
    theme_bw() +
    theme( legend.title = element_blank(),
           legend.position = c(.95, .99), 
           legend.justification = c("right", "top"), 
           legend.box.just = "right",
           #legend.background = element_rect( color = "black" ),
           legend.margin = margin( rep(6, 4) ) ) +
    scale_color_manual( values = c("firebrick", "navy"))
  
  
  ## Posterior distribution and HDI
  x_lower <- floor( qgamma( .0001, fit$posterior_pars["alpha_post"], fit$posterior_pars["beta_post"] ) )
  x_upper <- ceiling( qgamma( .9999, fit$posterior_pars["alpha_post"], fit$posterior_pars["beta_post"] ) )
  
  plot_data <- fit$posterior_density %>% filter( lambda >= x_lower & lambda <= x_upper )
  
  p2 <- ggplot( data = plot_data, aes(x = lambda ) ) +
    geom_line( aes(y = density ), col = "firebrick") +
    geom_area( stat = "function", fun = dgamma, 
               args = list( shape = fit$posterior_pars["alpha_post"], rate = fit$posterior_pars["beta_post"] ), 
               fill = "red", xlim = c(fit$posterior_ci["lower"], fit$posterior_ci["upper"]), alpha = .25 ) +
    labs( x = "Lambda", y = "Density", subtitle = "Posterior density") +
    theme_bw() +
    theme( legend.position = "bottom")
  
  post_samples <- sample_posterior( fit = fit )
  y <- table(fit$data) / sum(table(fit$data)) 
  data <- tibble( x = as.numeric( names(y)), y = y )
  
  p3 <- post_samples %>%
    ggplot( aes( x = sample, y = freq)) + 
    geom_point( alpha = .05, col = "grey") +
    geom_point( data = data, aes(x = x, y = y, col = "Data")) +
    labs( x = "K", y = "Density", subtitle = "Posterior predictive") +
    theme_bw() +
    theme( legend.title = element_blank(),
           legend.position = c(.95, .99), 
           legend.justification = c("right", "top"),
           legend.box.just = "right",
           legend.margin = margin(6,6,6,6)) +
    scale_color_manual( values = c("firebrick"))
  
  ggarrange( p1, p2, p3, nrow = 1)
  
}

bayes_summary <- function( fit )
{
  
  cat( paste0(
    "--------------------------------------------\n",
    "Gamma-Poisson Model\n--------------------------------------------\n",
    "Data",
    "\nObs: ", fit$n_obs,
    "\nSum: ", sum( fit$data ),
    "\nMLE: ", mean( fit$data ), 
    "\n\nDistributions",
    "\nPrior:\t   ", fit$prior, 
    "\nPosterior: ", fit$posterior,
    "\n\nPosterior Parameters",
    "\nMode: \t", round( fit$posterior_summary["mode"], 3 ),
    "\nMedian: ", round( fit$posterior_summary["median"], 3 ),
    "\nMean: \t", round( fit$posterior_summary["mean"], 3 ),
    "\nCI: \t", "[", round( fit$posterior_ci[1], 3 ), " ", round( fit$posterior_ci[2], 3 ), "]\n"
    
  ))
  
}


