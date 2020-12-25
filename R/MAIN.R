# Title: "A nonparametric test of constancy of a monotone function"
# Author: Avi Kenny, Marco Carone



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
  library(ggplot2)
  library(dplyr)
  library(boot)
}

# Load functions
{
  # source("function.R")
}

# Set code blocks to run
{
  run_main <- TRUE
  run_testing <- FALSE
}



############################################.
##### MAIN: Simulation setup (density) #####
############################################.

if (run_main) {
  
  # Use these commands to run on Slurm:
  # sbatch --export=run='first',cluster='bionic',type='R',project='z.monotest' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh
  # sbatch --depend=afterok:11 --array=1-540 --export=cluster='bionic',type='R',project='z.monotest' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh
  # sbatch --depend=afterok:12 --export=run='last',cluster='bionic',type='R',project='z.monotest' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh
  
  run_on_cluster(
    
    first = {
      
      sim <- new_sim()
      
      sim %<>% set_config(
        num_sim = 1000,
        parallel = "none",
        packages = c("dplyr", "boot")
      )
      
      sim %<>% set_levels(
        n = c(25,50,100),
        # n = c(10,20,30,40,50),
        test = list(
          # "Z test" = list(type="Z test"),
          # "L test" = list(type="L test"),
          "Slope (SS-adapted)" = list(
            type = "Slope",
            params = list(subtype="SS-adapted")
          ),
          "Slope (NEW bootstrap; quantile)" = list(
            type = "Slope",
            params=list(subtype="NEW bootstrap", ci_type="quantile")
          ),
          "Slope (full bootstrap; quantile)" = list(
            type = "Slope",
            params=list(subtype="full bootstrap", ci_type="quantile")
          ),
          # "Slope (full bootstrap; normal)" = list(
          #   type = "Slope",
          #   params=list(subtype="full bootstrap", ci_type="normal")
          # ),
          "Slope (mixed bootstrap; quantile)" = list(
            type = "Slope",
            params=list(subtype="mixed bootstrap", ci_type="quantile")
          )
          # "Slope (mixed bootstrap; normal)" = list(
          #   type = "Slope",
          #   params=list(subtype="mixed bootstrap", ci_type="normal")
          # )
          # "CumlIncr (delta=1/2)" = list(
          #   type = "CumlIncr",
          #   params = list(delta=(1/2), wts=1)
          # ),
          # "CumlIncr (delta=1/3, equal weights)" = list(
          #   type = "CumlIncr",
          #   params = list(delta=(1/3), wts=c(0.5,0.5))
          # )
        ),
        true_density = c("f(x)=1", "f(x)=2x", "f(x)=ke^x")
      )
      
      # This function generates data
      #   `n` is the sample size
      #   `true_density` is the type of density used to generate the data, specified
      #       as a character string
      sim %<>% add_creator("generate_dataset", function(n, true_density) {
        
        # !!!!! Add other monotone distributions
        
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
      
      # Marco method 1
      # Currently assumes equally sized intervals of size delta
      # wts argument can either be "equal" or a vector of length (1/delta)-1
      sim %<>% add_method("CumlIncr", function(dat, n, params) {
        
        delta <- params$delta
        if (params$wts[1]=="equal") {
          wts <- rep(1/(1/delta-1),1/delta-1)
        } else {
          wts <- params$wts
        }
        
        if (params$wts[1]=="equal") {
          
          T_n <- 1 - mean(dat<=1-delta) - mean(dat<=delta)
          asy_sd <- sqrt(2*delta)
          
        } else if (delta==1/2) {
          
          # !!!!! This is technically the same test as above; combine
          T_n <- 2 - 4*mean(dat<=0.5)
          asy_sd <- 2
          
        } else if (delta==1/3) {
          
          w1 <- wts[1]
          w2 <- wts[2]
          T_n <- 3*(w1-2*w2)*mean(dat<=2/3) + 3*(w2-2*w1)*mean(dat<=1/3) + 3*w2
          asy_sd <- sqrt( 6*(w1^2-w1*w2+w2^2) )
          
        }
        
        return (as.numeric(sqrt(n)*T_n>asy_sd*qnorm(0.95)))
        
      })
      
      # Generate beta_n distribution and add as a simulation constant
      beta_n_distr <- list()
      for (n in c(25,50,100)) { # !!!!!
      # for (n in c(10,20,30,40,50)) {
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
      sim %<>% add_constants("beta_n_distr"=beta_n_distr)
      rm(beta_n_distr)
      
      # Marco method 2
      sim %<>% add_method("Slope", function(dat, n, params) {
        
        Theta_hat <- ecdf(dat)
        mu_2n <- mean(dat^2)
        mu_3n <- mean(dat^3)
        beta_n <- mean((mu_2n*dat^2 - mu_3n*dat)*Theta_hat(dat))
        
        if (params$subtype=="SS-adapted") {
          crit_val <- quantile(C$beta_n_distr[[as.character(n)]],0.95)
          return(as.numeric(beta_n>crit_val))
        }
        
        if (params$subtype=="mixed bootstrap") {
          
          # Define the statistic to bootstrap
          bootstat <- function(dat,indices) {
            d <- dat[indices]
            Theta_hat <- ecdf(d)
            mu_2n <- mean(dat^2)
            mu_3n <- mean(dat^3)
            return (mean((mu_2n*dat^2 - mu_3n*dat)*(Theta_hat(dat)-dat)))
          }
          
          # Run bootstrap
          boot_obj <- boot(data=dat, statistic=bootstat, R=1000)
          
          # Calculate critical value
          if (params$ci_type=="quantile") {
            crit_val <- as.numeric(quantile(boot_obj$t, 0.05))
          } else if (params$ci_type=="normal") {
            crit_val <- qnorm(0.05, mean=mean(boot_obj$t), sd=sd(boot_obj$t))
          }
          
          return(as.numeric(crit_val>0))
          
        }
        
        if (params$subtype=="NEW bootstrap") {
          
          # Define the statistic to bootstrap
          bootstat <- function(dat,indices) {
            d <- dat[indices]
            Theta_hat <- ecdf(d)
            mu_2n <- mean(dat^2)
            mu_3n <- mean(dat^3)
            piece_1 <- mean((mu_2n*dat^2 - mu_3n*dat)*(Theta_hat(dat)-dat))
            piece_2a <- mean(dat^2*(Theta_hat(dat)-dat))*(mean(d^2)-mean(dat^2))
            piece_2b <- mean(dat*(dat-Theta_hat(dat)))*(mean(d^3)-mean(dat^3))
            return (piece_1+piece_2a+piece_2b)
          }
          
          # Run bootstrap
          boot_obj <- boot(data=dat, statistic=bootstat, R=1000)
          
          # Calculate critical value
          if (params$ci_type=="quantile") {
            crit_val <- as.numeric(quantile(boot_obj$t, 0.05))
          } else if (params$ci_type=="normal") {
            crit_val <- qnorm(0.05, mean=mean(boot_obj$t), sd=sd(boot_obj$t))
          }
          
          return(as.numeric(crit_val>0))
          
        }
        
        if (params$subtype=="full bootstrap") {
          
          # Define the statistic to bootstrap
          bootstat <- function(dat,indices) {
            d <- dat[indices]
            Theta_hat <- ecdf(d)
            mu_2n <- mean(d^2)
            mu_3n <- mean(d^3)
            return (mean((mu_2n*d^2 - mu_3n*d)*Theta_hat(d)))
          }
          
          # Run bootstrap
          boot_obj <- boot(data=dat, statistic=bootstat, R=1000)
          
          # Calculate critical value
          if (params$ci_type=="quantile") {
            crit_val <- as.numeric(quantile(boot_obj$t, 0.05))
          } else if (params$ci_type=="normal") {
            crit_val <- qnorm(0.05, mean=mean(boot_obj$t), sd=sd(boot_obj$t))
          }
          
          return(as.numeric(crit_val>0))
          
        }
        
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
      # !!!!! Does not work yet; gives too high rejection rates
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
      
      sim %<>% set_script(function() {
        
        dat <- generate_dataset(L$n, L$true_density)
        if (is.null(L$test$params)) {
          reject <- do.call(L$test$type, list(dat, L$n))
        } else {
          reject <- do.call(L$test$type, list(dat, L$n, L$test$params))
        }
        
        return (list("reject"=reject))
        
      })

    },
    
    main = {
      sim %<>% run()
    },
    
    last = {
      
      # sim <- readRDS("sim.simba")
      # sim <- readRDS("sim_20201126.simba")

      summ <- sim %>% summary(mean=list(name="power",x="reject"))
      
      # Visualize results
      ggplot(summ, aes(x=n, y=power, color=factor(test))) +
        geom_line() + geom_point() +
        facet_wrap(~true_density, ncol=3, scales="free_y") +
        labs(
          x = "sample size", y = "Power", color = "Test type",
          title = paste0("Power of tests for constant vs. monotone density (1,000",
                         " sims per level)")
        ) + scale_color_manual(values=c("turquoise", "salmon", "dodgerblue2", "green3", "darkorchid2", "orange"))
        
      
    },
    
    cluster_config = list(
      sim_var = "sim",
      js = "slurm",  # Bionic
      dir = "/home/akenny/z.monotest"  # Bionic
    )

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



###################################################.
##### MAIN: Sim #1: compare different methods #####
###################################################.

if (FALSE) {
  
  # !!!!! Restructure
  
}



#######################################################################.
##### MAIN: Sim #2: compare different weights for CumlIncr method #####
#######################################################################.

if (FALSE) {
  
  sim_2 <- sim
  
  sim_2 %<>% set_config(num_sim=10000)
  
  sim_2 %<>% set_levels(
    n = c(10,30,50),
    test = list(
      "0.0" = list(type="CumlIncr", params=list(delta=1/3,wts=c(1.0,0.0))),
      "0.1" = list(type="CumlIncr", params=list(delta=1/3,wts=c(0.9,0.1))),
      "0.2" = list(type="CumlIncr", params=list(delta=1/3,wts=c(0.8,0.2))),
      "0.3" = list(type="CumlIncr", params=list(delta=1/3,wts=c(0.7,0.3))),
      "0.4" = list(type="CumlIncr", params=list(delta=1/3,wts=c(0.6,0.4))),
      "0.5" = list(type="CumlIncr", params=list(delta=1/3,wts=c(0.5,0.5))),
      "0.6" = list(type="CumlIncr", params=list(delta=1/3,wts=c(0.4,0.6))),
      "0.7" = list(type="CumlIncr", params=list(delta=1/3,wts=c(0.3,0.7))),
      "0.8" = list(type="CumlIncr", params=list(delta=1/3,wts=c(0.2,0.8))),
      "0.9" = list(type="CumlIncr", params=list(delta=1/3,wts=c(0.1,0.9))),
      "1.0" = list(type="CumlIncr", params=list(delta=1/3,wts=c(0.0,1.0)))
    ),
    true_density = c("f(x)=1", "f(x)=2x", "f(x)=ke^x")
  )
  
  sim_2 %<>% run("one_simulation")
  
  summ <- sim_2 %>% summary() %>% rename("power"=`mean_reject`)

  # Visualize results
  ggplot(summ, aes(x=test, y=power, group=(n), color=factor(n))) +
    geom_line() + geom_point() +
    facet_wrap(~true_density, ncol=3, scales="free_y") +
    labs(
      x = "w_2 (weight)", y = "Power", color = "n", group = "n",
      title = paste0("Power of CumlIncr test for constant vs. monotone density",
                     " (10,000 sims per level)")
    )
  
}



############################################################################.
##### MAIN: Sim #3: compare different delta values for CumlIncr method #####
############################################################################.

if (FALSE) {
  
  sim_3 <- sim
  
  sim_3 %<>% set_config(num_sim=10000)
  
  sim_3 %<>% set_levels(
    n = c(10,30,50),
    test = list(
      "1/2" = list(type="CumlIncr", params=list(delta=1/2,wts="equal")),
      "1/3" = list(type="CumlIncr", params=list(delta=1/3,wts="equal")),
      "1/4" = list(type="CumlIncr", params=list(delta=1/4,wts="equal")),
      "1/5" = list(type="CumlIncr", params=list(delta=1/5,wts="equal")),
      "1/6" = list(type="CumlIncr", params=list(delta=1/6,wts="equal")),
      "1/7" = list(type="CumlIncr", params=list(delta=1/7,wts="equal")),
      "1/8" = list(type="CumlIncr", params=list(delta=1/8,wts="equal")),
      "1/9" = list(type="CumlIncr", params=list(delta=1/9,wts="equal")),
      "1/10" = list(type="CumlIncr", params=list(delta=1/10,wts="equal"))
      # "1/4" = list(type="CumlIncr", params=list(delta=1/4,wts="equal")),
      # "1/8" = list(type="CumlIncr", params=list(delta=1/8,wts="equal")),
      # "1/16" = list(type="CumlIncr", params=list(delta=1/16,wts="equal")),
      # "1/32" = list(type="CumlIncr", params=list(delta=1/32,wts="equal")),
      # "1/64" = list(type="CumlIncr", params=list(delta=1/64,wts="equal")),
      # "1/128" = list(type="CumlIncr", params=list(delta=1/128,wts="equal"))
    ),
    true_density = c("f(x)=1", "f(x)=2x", "f(x)=ke^x")
  )
  
  sim_3 %<>% run("one_simulation")
  
  summ <- sim_3 %>% summary() %>% rename("power"=`mean_reject`)
  summ$test <- factor(summ$test, levels=c("1/2","1/3","1/4","1/5","1/6","1/7",
                                          "1/8","1/9","1/10"))
  # summ$test <- factor(summ$test, levels=c("1/2","1/4","1/8","1/16","1/32",
  #                                         "1/64","1/128"))
  
  # Visualize results
  ggplot(summ, aes(x=test, y=power, group=(n), color=factor(n))) +
    geom_line() + geom_point() +
    facet_wrap(~true_density, ncol=3, scales="free_y") +
    labs(
      x = "delta", y = "Power", color = "n", group = "n",
      title = paste0("Power of CumlIncr test for constant vs. monotone density",
                     " (10,000 sims per level)")
    )
  
}
