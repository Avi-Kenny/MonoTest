#' Generate data according to a density
#' 
#' @param n Sample size
#' @param true_density A character string representing the true density
#' @return A vector; a sample from the specified density
generate_data <- function(n, true_density) {
  
  # !!!!! Add other monotone distributions
  
  # Q is the quantile function for probability integral transform sampling
  if (true_density=="f(x)=1") {
    Q <- function(x) { x }
  } else if (true_density=="f(x)=2x") {
    Q <- function(x) { sqrt(x) }
  } else if (true_density=="f(x)=ke^x") {
    Q <- function(x) { log(x*(exp(1)-1)+1) }
  }
  
  dat = Q(runif(n));

}



#' Hypothesis testing method 1: cumulative increments
#' 
#' @param dat Data returned by generate_data()
#' @param params A list, containing the following:
#'   - `delta` is the interval size (e.g. 1/3).
#'   - `wts` can either be "equal" or a a vector of length (1/delta)-1
#' @return Binary; is null rejected (1) or not (0)
#' @notes
#'   - Currently assumes equally sized intervals of size delta
cuml_incr <- function(dat, params) {

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
  
  n <- length(dat)
  
  return (as.numeric(sqrt(n)*T_n>asy_sd*qnorm(0.95)))
  
}



#' Hypothesis testing method 2: slope
#' 
#' @param dat Data returned by generate_data()
#' @param params A list, containing the following:
#'   - `subtype` One of c("SS-adapted", "mixed bootstrap", "full bootstrap")
#'   - `ci_type` (optional); one of c("quantile", "normal"); if a bootstrap CI
#'       is used, should the cutoff points be based on the quantiles or the
#'       Normal approximation?
#' @return Binary; is null rejected (1) or not (0)
#' @notes
#'   - The "SS-adapted" subtype references the simulation constant
#'     C$beta_n_distr, which needs to be adapted for each sample size specified
#'     via the simulation levels
slope <- function(dat, params) {
  
  n <- length(dat)
  
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
  
}



#' Z test from Woodroofe and Sun 1999
#' 
#' @param dat Data returned by generate_data()
#' @param params Unused
#' @return Binary; is null rejected (1) or not (0)
ws_test_z <- function(dat, params) {
  n <- length(dat)
  return(as.numeric(mean(log(dat))>-1+qnorm(0.95)/sqrt(n)))
}



#' L test from Woodroofe and Sun 1999
#' 
#' @param dat Data returned by generate_data()
#' @param params Unused
#' @return Binary; is null rejected (1) or not (0)
ws_test_l <- function(dat, params) {
  n <- length(dat)
  return(as.numeric(mean(dat)>1/2+qnorm(0.95)/sqrt(12*n)))
}
