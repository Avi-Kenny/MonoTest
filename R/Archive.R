# The "NEW bootstrap" method used a second-order term to increase power
# This was part of the slope test
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



# D test from Woodroofe and Sun 1999 (c=0.2)
# !!!!! Does not work yet; gives too high rejection rates
sim %<>% add_method("ws_test_d", function(dat, n) {
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
