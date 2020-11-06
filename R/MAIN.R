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
  setwd("z.stepped.wedge/R")
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



########################.
##### MAIN: Test 1 #####
########################.

if (run_main) {
  
  sim <- new_sim()
  
  sim %<>% set_config(
    num_sim = 10000, # !!!!!
    parallel = "outer",
    packages = c("magrittr", "dplyr")
  )
  
  sim %<>% set_levels(
    n = c(20, 30),
    test = c("test_t", "test_z", "test_l"),
    true_density = c("null", "sqrt")
  )
  
  # This function generates data
  #   `n` is the sample size
  #   `true_density` is the type of density used to generate the data
  sim %<>% add_creator("generate_dataset", function(n, true_density) {
    
    if (true_density=="null") {
      Q <- function(x) { x }
    } else if (true_density=="sqrt") {
      Q <- function(x) { sqrt(x) }
    }
    
    dat = Q(runif(n));
    
  })
  
  sim %<>% add_method("test_t", function(dat, n) {
    as.numeric(sqrt(n)*(1-2*mean(dat<=0.5))>qnorm(0.95))
  })
  
  sim %<>% add_method("test_z", function(dat, n) {
    as.numeric(mean(log(dat))>-1+qnorm(0.95)/sqrt(n))
  })
  
  sim %<>% add_method("test_l", function(dat, n) {
    as.numeric(mean(dat)>1/2+qnorm(0.95)/sqrt(12*n))
  })
  
  sim %<>% add_script("one_simulation", function() {
    
    dat <- generate_dataset(L$n, L$true_density)
    reject <- do.call(L$test, list(dat, L$n))
    return (list("reject"=reject))
    
  })
  
  sim %<>% run("one_simulation")
  
  sim %>% summary() %>% rename("power"=`mean_reject`)
  
}
