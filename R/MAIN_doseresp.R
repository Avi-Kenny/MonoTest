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
  library(mgcv)
}

# Load functions
{
  source("fns_doseresp.R")
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
      # n = 1000,
      n = c(500,1000),
      beta = c(0,0.7),
      mono_form = c("identity"),
      # mono_form = c("identity", "square", "sqrt", "step_0.2"),
      test = list(
        "Wald" = list(type="test_regression", params=NULL),
        "Slope: glm" = list(type="slope_dr", params=list(est="glm"))
        # "Slope: sm spline" = list(type="slope_dr", params=list(est="sm spline"))
      )
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
  # sbatch --depend=afterok:11 --array=1-48 --export=run='main',cluster='bionic',type='R',project='z.monotest' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh
  # sbatch --depend=afterok:58408889 --export=run='last',cluster='bionic',type='R',project='z.monotest' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh
  
  run_on_cluster(
    
    first = {
      
      sim <- new_sim()
      
      sim %<>% set_config(
        num_sim = 500,
        # parallel = "outer",
        # stop_at_error = TRUE, # !!!!!
        packages = c("dplyr", "boot", "mgcv")
      )
      
      sim <- do.call(set_levels, c(list(sim), level_set))
      
      # Add functions to simulation object
      sim %<>% add_creator(generate_data_dr)
      sim %<>% add_method(expit)
      sim %<>% add_method(intexpit)
      sim %<>% add_method(test_regression)
      sim %<>% add_method(slope_dr)

      sim %<>% set_script(function() {
        
        dat <- generate_data_dr(L$n, L$beta, L$mono_form)
        reject <- do.call(L$test$type, list(dat, L$test$params))
        
        return (list("reject"=reject))
        
      })
      
    },
    
    main = {
      sim %<>% run()
      # sim %<>% update()
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
  
  summ <- sim %>% summary(mean=list(name="power", x="reject"))
  
  # Visualize results
  ggplot(summ, aes(x=n, y=power, color=factor(test))) +
    geom_line() + geom_point() +
    facet_grid(rows=vars(beta), cols=vars(mono_form)) + # scales="free_y"
    labs(
      x = "sample size", y = "Power", color = "Test type",
      title = paste0("Power of tests for constant vs. monotone causal ",
                     "dose-response curve (1,000 sims per level)")
    ) + scale_color_manual(values=c("turquoise", "salmon", "dodgerblue2",
                                    "green3", "darkorchid2", "orange"))
  
}
