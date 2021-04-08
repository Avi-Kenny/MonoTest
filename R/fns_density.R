#' Generate appx distribution of beta_n via Monte Carlo sampling
#' 
#' @param n_list Sample sizes
#' @return A list of samples corresponding to each value in n_list
beta_n_distr <- function(n_list) {
  
  beta_n_distr_ls <- list()
  
  for (n in n_list) {
    beta_ns <- c()
    for (i in 1:1000) {
      x <- runif(n)
      Theta_hat <- ecdf(x)
      beta_n <- mean((mean(x^2)*x^2 - mean(x^3)*x)*Theta_hat(x))
      beta_ns <- c(beta_ns,beta_n)
    }
    beta_n_distr_ls[[as.character(n)]] <- beta_ns
  }
  
  return(beta_n_distr_ls)
  
}



#' Hypothesis testing method 1: cumulative increments
#' 
#' @param dat Data returned by generate_data()
#' @param alt_type Type of alternative hypothesis; either "incr" or "decr"
#' @param params A list, containing the following:
#'   - `delta` is the interval size (e.g. 1/3).
#'   - `wts` can either be "equal" or a a vector of length (1/delta)-1
#' @return Binary; is null rejected (1) or not (0)
#' @notes
#'   - Currently assumes equally sized intervals of size delta
cuml_incr <- function(dat, alt_type, params) {

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
#' @param alt_type Type of alternative hypothesis; either "incr" or "decr"
#' @param params A list, containing the following:
#'   - `p_star` Which distribution to use for p_star; one of c("U(0,1)", "P_0")
#'   - `subtype` One of c("SS-adapted", "mixed bootstrap", "full bootstrap")
#'   - `ci_type` (optional); one of c("quantile", "normal"); if a bootstrap CI
#'       is used, should the cutoff points be based on the quantiles or the
#'       Normal approximation?
#' @return Binary; is null rejected (1) or not (0)
#' @notes
#'   - The "SS-adapted" subtype references the simulation constant
#'     C$beta_n_distr, which needs to be adapted for each sample size specified
#'     via the simulation levels
slope <- function(dat, alt_type, params) {
  
  n <- length(dat)
  
  Theta_hat <- ecdf(dat)
  mu_2n <- mean(dat^2)
  mu_3n <- mean(dat^3)
  beta_n <- mean((mu_2n*dat^2 - mu_3n*dat)*Theta_hat(dat))
  
  if (params$subtype=="SS-adapted") {
    
    if (alt_type=="incr") {
      crit_val <- quantile(C$beta_n_distr[[as.character(n)]],0.95)
      return(as.numeric(beta_n>crit_val))
    }
    if (alt_type=="decr") {
      crit_val <- quantile(C$beta_n_distr[[as.character(n)]],0.05)
      return(as.numeric(beta_n<crit_val))
    }
    
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
    # !!!!! Adapt for alt_type=="decr"
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
    # !!!!! Adapt for alt_type=="decr"
    if (params$ci_type=="quantile") {
      crit_val <- as.numeric(quantile(boot_obj$t, 0.05))
    } else if (params$ci_type=="normal") {
      crit_val <- qnorm(0.05, mean=mean(boot_obj$t), sd=sd(boot_obj$t))
    }
    
    return(as.numeric(crit_val>0))
    
  }
  
  if (params$subtype=="asymptotic") {
    
    # Fix m at 10,000
    m <- 10000
    
    # Sample from P_0^*
    if (params$p_star=="P_0") {
      m_sample <- dat
    }
    if (params$p_star=="U(0,1)") {
      m_sample <- runif(n=m, min=0, max=1)
    }
    if (params$p_star=="B(0.4,0.4)") {
      m_sample <- rbeta(n=m, shape1=0.4, shape2=0.4)
    }
    # if (params$p_star=="N(0.5,0.1)") {         # !!!!!
    #   m_sample <- rnorm(n=m, mean=0.5, sd=0.1) # !!!!!
    # }                                          # !!!!!
    if (!(params$p_star %in% c("U(0,1)","P_0","B(0.4,0.4)"))) {
      stop("p_star must be in c('U(0,1)','P_0','B(0.4,0.4)')")
    }
    
    # Calculate test statistic
    # !!!!! Move calculation of test stat up
    lambda_2 <- mean(m_sample^2)
    lambda_3 <- mean(m_sample^3)
    theta_n <- ecdf(dat)
    beta_nm <- mean((lambda_2*m_sample^2-lambda_3*m_sample)*theta_n(m_sample))
    
    var_est <- mean(sapply(dat, function(x) {
      (mean((
        (as.numeric(x<=m_sample)-m_sample) *
          (lambda_2*m_sample^2-lambda_3*m_sample)
      )))^2
    }))
    
    # Calculate critical value
    z <- (sqrt(n)*beta_nm) / sqrt(var_est)
    if (alt_type=="incr") {
      crit_val <- qnorm(0.95)
      return(as.numeric(z>crit_val))
    }
    if (alt_type=="decr") {
      crit_val <- qnorm(0.05)
      return(as.numeric(z<crit_val))
    }
    
  }
  
}



#' Z test from Woodroofe and Sun 1999
#' 
#' @param dat Data returned by generate_data()
#' @param alt_type Currently only written for "incr"
#' @param params list(crit_val=cv); cv is either "normal appx" or "simulated"
#' @return Binary; is null rejected (1) or not (0)
ws_test_z <- function(dat, alt_type, params) {
  n <- length(dat)
  test_stat_z <- mean(log(dat))
  if (params$crit_val=="normal appx") {
    crit_val <- -1 + qnorm(0.95)/sqrt(n)
  }
  if (params$crit_val=="simulated") {
    # These were pre-calculated with a sample size of 1,000,000
    if (n %in% c(10,20,30,40,50)) {
      crit_val <- case_when(
        n==10 ~ -0.5422811,
        n==20 ~ -0.6630681,
        n==30 ~ -0.7199029,
        n==40 ~ -0.7549827,
        n==50 ~ -0.7789746
      )
    } else {
      crit_val <- quantile(sapply(c(1:10000), function(x) {
        mean(log(runif(n)))
      }),0.95)
    }
  }
  reject <- as.numeric(test_stat_z>crit_val)
  return(reject)
}



#' L test from Woodroofe and Sun 1999
#' 
#' @param dat Data returned by generate_data()
#' @param alt_type Currently only written for "incr"
#' @param params list(crit_val=cv); cv is either "normal appx" or "simulated"
#' @return Binary; is null rejected (1) or not (0)
ws_test_l <- function(dat, alt_type, params) {
  n <- length(dat)
  test_stat_l <- mean(dat)
  if (params$crit_val=="normal appx") {
    crit_val <- 0.5 + qnorm(0.95)/(sqrt(12*n))
  }
  if (params$crit_val=="simulated") {
    # These were pre-calculated with a sample size of 1,000,000
    if (n %in% c(10,20,30,40,50)) {
      crit_val <- case_when(
        n==10 ~ 0.6505658,
        n==20 ~ 0.6062897,
        n==30 ~ 0.586739,
        n==40 ~ 0.5751417,
        n==50 ~ 0.5672297
      )
    } else {
      crit_val <- quantile(sapply(c(1:10000), function(x) {
        mean(runif(n))
      }),0.95)
    }
  }
  reject <- as.numeric(test_stat_l>crit_val)
  return(reject)
}
