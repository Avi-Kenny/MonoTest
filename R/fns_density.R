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



#' Z test from Woodroofe and Sun 1999
#' 
#' @param dat Data returned by generate_data()
#' @param alt_type Currently only written for "incr"
#' @param params list(crit_val=cv); cv is either "normal appx" or "simulated"
#' @return Binary; is null rejected (1) or not (0)
test_ws_z <- function(dat, alt_type="incr", params) {
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
test_ws_l <- function(dat, alt_type="incr", params) {
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
