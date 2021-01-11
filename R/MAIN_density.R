# Title: "A nonparametric test of constancy of a monotone function"
# Author: Avi Kenny, Marco Carone



#################.
##### SETUP #####
#################.

# Set working directory
if (Sys.getenv("USERDOMAIN")=="AVI-KENNY-T460") {
  setwd(paste0("C:/Users/avike/OneDrive/Desktop/Avi/Biostats + Research/Resear",
               "ch/Marco Carone/Project - monotone testing/z.monotest/R"))
} else {
  setwd("z.monotest/R")
}

# Load packages
{
  library(simba) # devtools::install_github(repo="Avi-Kenny/simba")
  library(ggplot2)
  library(dplyr)
  library(boot)
}

# Load functions
{
  source("fns_density.R")
}

# Set code blocks to run
{
  run_main <- TRUE
  run_testing <- FALSE
}



##########################################################.
##### MAIN: Set level sets for different simulations #####
##########################################################.

if (run_main) {
  if (Sys.getenv("run") %in% c("first", "")) {

    # Simulation 1: compare all methods
    level_set_1 <- list(
      n = c(10,20,30,40,50),
      test = list(
        "Z test" = list(type="ws_test_z"),
        "L test" = list(type="ws_test_l"),
        "slope (SS-adapted)" = list(
          type = "slope",
          params = list(subtype="SS-adapted")
        ),
        "slope (full bootstrap; quantile)" = list(
          type = "slope",
          params = list(subtype="full bootstrap", ci_type="quantile")
        ),
        "slope (mixed bootstrap; quantile)" = list(
          type = "slope",
          params = list(subtype="mixed bootstrap", ci_type="quantile")
        ),
        "cuml_incr (delta=1/2)" = list(
          type = "cuml_incr",
          params = list(delta=(1/2), wts=1)
        ),
        "cuml_incr (delta=1/3, equal weights)" = list(
          type = "cuml_incr",
          params = list(delta=(1/3), wts=c(0.5,0.5))
        )
      ),
      true_density = c("f(x)=1", "f(x)=2x", "f(x)=ke^x")
    )
    
    # Simulation 2: different weights for cuml_incr
    level_set_2 <- list(
      n = c(10,30,50),
      test = list(
        "0.0" = list(type="cuml_incr", params=list(delta=1/3,wts=c(1.0,0.0))),
        "0.1" = list(type="cuml_incr", params=list(delta=1/3,wts=c(0.9,0.1))),
        "0.2" = list(type="cuml_incr", params=list(delta=1/3,wts=c(0.8,0.2))),
        "0.3" = list(type="cuml_incr", params=list(delta=1/3,wts=c(0.7,0.3))),
        "0.4" = list(type="cuml_incr", params=list(delta=1/3,wts=c(0.6,0.4))),
        "0.5" = list(type="cuml_incr", params=list(delta=1/3,wts=c(0.5,0.5))),
        "0.6" = list(type="cuml_incr", params=list(delta=1/3,wts=c(0.4,0.6))),
        "0.7" = list(type="cuml_incr", params=list(delta=1/3,wts=c(0.3,0.7))),
        "0.8" = list(type="cuml_incr", params=list(delta=1/3,wts=c(0.2,0.8))),
        "0.9" = list(type="cuml_incr", params=list(delta=1/3,wts=c(0.1,0.9))),
        "1.0" = list(type="cuml_incr", params=list(delta=1/3,wts=c(0.0,1.0)))
      ),
      true_density = c("f(x)=1", "f(x)=2x", "f(x)=ke^x")
    )
    
    # Simulation 3: different delta values for cuml_incr
    level_set_3 <- list(
      n = c(10,30,50),
      test = list(
        "1/2" = list(type="cuml_incr", params=list(delta=1/2,wts="equal")),
        "1/3" = list(type="cuml_incr", params=list(delta=1/3,wts="equal")),
        "1/4" = list(type="cuml_incr", params=list(delta=1/4,wts="equal")),
        "1/5" = list(type="cuml_incr", params=list(delta=1/5,wts="equal")),
        "1/6" = list(type="cuml_incr", params=list(delta=1/6,wts="equal")),
        "1/7" = list(type="cuml_incr", params=list(delta=1/7,wts="equal")),
        "1/8" = list(type="cuml_incr", params=list(delta=1/8,wts="equal")),
        "1/9" = list(type="cuml_incr", params=list(delta=1/9,wts="equal")),
        "1/10" = list(type="cuml_incr", params=list(delta=1/10,wts="equal"))
      ),
      true_density = c("f(x)=1", "f(x)=2x", "f(x)=ke^x")
    )
    
    # Simulation 4: !!!!!
    level_set_4 <- list(
      # !!!!!
    )
    
  }
}



################################################.
##### MAIN: Choose which simulation to run #####
################################################.

if (run_main) {
  if (Sys.getenv("run") %in% c("first", "")) {
    
    # Set this manually
    level_set <- level_set_1
    
  }
}



##########################################.
##### MAIN: Setup and run simulation #####
##########################################.

if (run_main) {
  
  # Use these commands to run on Slurm:
  # sbatch --export=run='first',cluster='bionic',type='R',project='z.monotest' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh
  # sbatch --depend=afterok:11 --array=1-540 --export=cluster='bionic',type='R',project='z.monotest' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh
  # sbatch --depend=afterok:12 --export=run='last',cluster='bionic',type='R',project='z.monotest' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh
  
  run_on_cluster(
    
    first = {
      
      sim <- new_sim()
      
      sim %<>% set_config(
        num_sim = 2, # !!!!!
        parallel = "none",
        stop_at_error = TRUE,
        packages = c("dplyr", "boot")
      )
      
      sim <- do.call(set_levels, c(list(sim), level_set))
      
      # Add functions to simulation object
      sim %<>% add_creator(generate_data)
      sim %<>% add_method(cuml_incr)
      sim %<>% add_method(slope)
      sim %<>% add_method(ws_test_z)
      sim %<>% add_method(ws_test_l)
      
      # Generate beta_n distribution and add as a simulation constant
      # Note: the n_list values must include the n-values in the sim levels
      beta_n_distr <- list()
      n_list <- c(10,20,30,40,50)
      for (n in n_list) {
        beta_ns <- c()
        for (i in 1:1000) {
          x <- runif(n)
          Theta_hat <- ecdf(x)
          beta_n <- mean((mean(x^2)*x^2 - mean(x^3)*x)*Theta_hat(x))
          beta_ns <- c(beta_ns,beta_n)
        }
        beta_n_distr[[as.character(n)]] <- beta_ns
      }
      sim %<>% add_constants("beta_n_distr"=beta_n_distr)
      rm(beta_n_distr)
      
      sim %<>% set_script(function() {
        
        dat <- generate_data(L$n, L$true_density)
        reject <- do.call(L$test$type, list(dat, L$test$params))
        
        return (list("reject"=reject))
        
      })

    },
    
    main = {
      sim %<>% run()
    },
    
    last = {},
    
    cluster_config = list(
      js = "slurm",
      dir = "/home/akenny/z.monotest"
    )

  )

}



######################################################.
##### MAIN: Process sim #1 (compare all methods) #####
######################################################.

if (FALSE) {
  
  # sim <- readRDS("sim.simba")

  summ <- sim %>% summary(mean=list(name="power",x="reject"))
  
  # Visualize results
  ggplot(summ, aes(x=n, y=power, color=factor(test))) +
    geom_line() + geom_point() +
    facet_wrap(~true_density, ncol=3, scales="free_y") +
    labs(
      x = "sample size", y = "Power", color = "Test type",
      title = paste0("Power of tests for constant vs. monotone density ",
                     "(1,000 sims per level)")
    ) + scale_color_manual(values=c("turquoise", "salmon", "dodgerblue2",
                                    "green3", "darkorchid2", "orange"))
  
}



##################################################################.
##### MAIN: Process sim #2 (different weights for cuml_incr) #####
##################################################################.

if (FALSE) {
  
  # sim <- readRDS("sim.simba")
  
  summ <- sim %>% summary() %>% rename("power"=`mean_reject`)

  # Visualize results
  ggplot(summ, aes(x=test, y=power, group=(n), color=factor(n))) +
    geom_line() + geom_point() +
    facet_wrap(~true_density, ncol=3, scales="free_y") +
    labs(
      x = "w_2 (weight)", y = "Power", color = "n", group = "n",
      title = paste0("Power of cuml_incr test for constant vs. monotone ",
                     "density (10,000 sims per level)")
    )
  
}



######################################################################.
##### MAIN: Process sim #3: different delta values for cuml_incr #####
######################################################################.

if (FALSE) {
  
  # sim <- readRDS("sim.simba")
  
  summ <- sim_3 %>% summary() %>% rename("power"=`mean_reject`)
  summ$test <- factor(summ$test, levels=c("1/2","1/3","1/4","1/5","1/6","1/7",
                                          "1/8","1/9","1/10"))
  
  # Visualize results
  ggplot(summ, aes(x=test, y=power, group=(n), color=factor(n))) +
    geom_line() + geom_point() +
    facet_wrap(~true_density, ncol=3, scales="free_y") +
    labs(
      x = "delta", y = "Power", color = "n", group = "n",
      title = paste0("Power of cuml_incr test for constant vs. monotone ",
                     "density (10,000 sims per level)")
    )
  
}



############################################################.
##### TESTING: compare distributions of test statistic #####
############################################################.

if (run_testing) {
  
  n <- 50
  reps <- 1000
  distribution_exact <- c()
  for (i in 1:reps) {
    x <- runif(n)
    Theta_hat <- ecdf(x)
    mu_2n <- mean(x^2)
    mu_3n <- mean(x^3)
    beta_n <- mean((mu_2n*x^2 - mu_3n*x)*Theta_hat(x))
    distribution_exact <- c(distribution_exact,beta_n)
  }
  
  beta_n <- function(dat,indices) {
    x <- dat[indices]
    Theta_hat <- ecdf(x)
    mu_2n <- mean(x^2)
    mu_3n <- mean(x^3)
    return (mean((mu_2n*x^2 - mu_3n*x)*Theta_hat(x)))
  }
  boot_obj <- boot(data=runif(n), statistic=beta_n, R=reps)
  
  beta_n_mixed <- function(dat,indices) {
    x <- dat[indices]
    Theta_hat <- ecdf(x)
    mu_2n <- mean(x^2)
    mu_3n <- mean(x^3)
    return (mean((mu_2n*x^2 - mu_3n*x)*(Theta_hat(x)-x)))
  }
  boot_obj_mixed <- boot(data=runif(n), statistic=beta_n_mixed, R=reps)
  
  print(sd(distribution_exact))
  print(sd(boot_obj$t))
  print(sd(boot_obj_mixed$t))
  
  ggplot(
    data.frame(
      x = c(distribution_exact,boot_obj$t,boot_obj_mixed$t),
      type = rep(c("Exact", "Full bootstrap", "Mixed bootstrap"), each=reps)
    ),
    aes(x=x)
  ) + geom_histogram() + facet_wrap(~type)
  
}
