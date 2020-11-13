# Title: "A nonparametric test of constancy of a monotone function"
# Author: Avi Kenny, Marco Carone



################.
##### INFO #####
################.

# !!!!! TO DO



################.
##### SETUP #####
################.

# Set working directory
if (Sys.getenv("USERDOMAIN")=="AVI-KENNY-T460") {
  setwd("C:/Users/avike/OneDrive/Desktop/Avi/Biostats + Research/Research/Marco Carone/Project - monotone testing/z.monotest/R")
} else {
  setwd("z.monotest/R")
}

# Load packages
{
  library(simba) # devtools::install_github(repo="Avi-Kenny/simba")
  library(magrittr)
  library(ggplot2)
  library(dplyr)
}

# Load functions
{
  # source("function.R")
}

# Set code blocks to run
{
  run_main <- TRUE
}



##############################.
##### MAIN: Density test #####
##############################.

if (run_main) {
  
  sim <- new_sim()
  
  sim %<>% set_config(
    num_sim = 10000,
    parallel = "none",
    packages = c("magrittr", "dplyr")
  )
  
  sim %<>% set_levels(
    n = c(10,20,30,40,50),
    test = c("Method 1 (K=1)", "Z test", "L test", "D test", "Slope"),
    true_density = c("f(x)=1", "f(x)=2x","f(x)=ke^x")
  )
  
  # This function generates data
  #   `n` is the sample size
  #   `true_density` is the type of density used to generate the data
  sim %<>% add_creator("generate_dataset", function(n, true_density) {
    
    # !!!!! Add beta distributions
    
    # Quantile function for probability integral transform sampling
    if (true_density=="f(x)=1") {
      Q <- function(x) { x }
    } else if (true_density=="f(x)=2x") {
      Q <- function(x) { sqrt(x) }
    } else if (true_density=="f(x)=ke^x") {
      Q <- function(x) { log(x*(exp(1)-1)+1) }
    }
    
    dat = Q(runif(n));
    
  })
  
  sim %<>% add_method("Method 1 (K=1)", function(dat, n) {
    as.numeric(sqrt(n)*(2-4*mean(dat<=0.5))>2*qnorm(0.95))
  })
  
  # !!!!! Generate beta_n distribution
  beta_n_distr <- list()
  for (n in c(10,20,30,40,50)) {
    beta_ns <- c()
    for (i in 1:1000) {
      x <- runif(n)
      Theta_hat <- ecdf(x)
      mu_2n <- mean(x^2)
      mu_3n <- mean(x^3)
      beta_n <- mean((mu_2n*x^2 - mu_3n*x)*Theta_hat(x))
      beta_ns <- c(beta_ns,beta_n)
    }
    beta_n_distr[[as.character(n)]] <- beta_ns
  }
  
  sim %<>% add_constant("beta_n_distr", beta_n_distr)
  
  sim %<>% add_method("Slope", function(dat, n) {
    Theta_hat <- ecdf(dat)
    mu_2n <- mean(dat^2)
    mu_3n <- mean(dat^3)
    beta_n <- mean((mu_2n*dat^2 - mu_3n*dat)*Theta_hat(dat))
    crit_val <- quantile(beta_n_distr[[as.character(n)]],0.95)
    return(as.numeric(beta_n>crit_val))
  })
  
  # Z test from Woodroofe and Sun 1999
  sim %<>% add_method("Z test", function(dat, n) {
    as.numeric(mean(log(dat))>-1+qnorm(0.95)/sqrt(n))
  })
  
  # L test from Woodroofe and Sun 1999
  sim %<>% add_method("L test", function(dat, n) {
    as.numeric(mean(dat)>1/2+qnorm(0.95)/sqrt(12*n))
  })
  
  # D test from Woodroofe and Sun 1999 (c=0.2)
  # !!!!! Testing
  sim %<>% add_method("D test", function(dat, n) {
    t <- seq(0,1,0.001)
    # t <- seq(0,1,0.01)
    F_hat <- ecdf(dat)
    D <- sqrt(n) * max(abs(F_hat(t)-t))
    if (!(n %in% c(10,20,30,40,50))) { stop("n must be in c(10,20,30,40,50)") }
    crit_val <- case_when(
      n==10 ~ 0.960,
      n==20 ~ 0.929,
      n==30 ~ 0.943,
      n==40 ~ 0.945,
      n==50 ~ 0.956
    )
    return (as.numeric(D>crit_val))
  })
  
  # P test from Woodroofe and Sun 1999 (c=0.2)
  sim %<>% add_method("P test", function(dat, n) {
    # !!!!! TO DO
  })
  
  sim %<>% add_script("one_simulation", function() {
    
    beta_n_distr <- C$beta_n_distr
    dat <- generate_dataset(L$n, L$true_density)
    reject <- do.call(L$test, list(dat, L$n))
    return (list("reject"=reject))
    
  })
  
  sim %<>% run("one_simulation")
  
  sim$constants <- NULL # !!!!! TEMP FIX
  summ <- sim %>% summary() %>%
    rename("power"=`mean_reject`) %>%
    arrange(true_density, n, test)
  
  # Visualize results
  ggplot(summ, aes(x=n, y=power, color=factor(test))) +
    geom_line() + geom_point() +
    # facet_wrap(~true_density, ncol=3) +
    facet_wrap(~true_density, ncol=3, scales="free_y") +
    labs(
      x = "sample size", y = "Power", color = "Test type",
      title = "Power of various tests for monotone density testing (10,000 replicates per level)"
    )

}
