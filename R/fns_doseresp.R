#' Generate data according to a density
#' 
#' @param n Sample size
#' @param beta Height of the dose-response curve
#' @param mono_form Functional form of dose-response curve; should be a monotone
#'     increasing function with domain [0,1] such that f(0)=0 and f(1)=1;
#'     choices are c("identity", "square", "sqrt", "step_0.2", "step_0.8")
#' @return A dataframe representing a study population
generate_data_dr <- function(n, beta, mono_form) {
  
  expit <- function(x) {exp(x)/(1+exp(x))}
  
  # Fix parameters
  beta0 <- -3
  beta1 <- 0.03
  beta2 <- 0.7
  # beta0 <- 0.2
  # beta1 <- 0.003
  # beta2 <- 0.1
  
  # Sample baseline covariates
  # `antib` varies from zero to one
  bmi <- rnorm(n, mean=25, sd=4)
  sex <- rbinom(n, size=1, prob=0.5)
  shape_2 <- ifelse(sex==1, 1.5, 1.1)
  antib <- rbeta(n, shape1=0.9, shape2=shape_2)
  
  # Set transformation function
  if (mono_form=="identity") {
    mono_f <- function(x) {x}
  } else if (mono_form=="square") {
    mono_f <- function(x) {x^2}
  } else if (mono_form=="sqrt") {
    mono_f <- function(x) {sqrt(x)}
  } else if (mono_form=="step_0.2") {
    mono_f <- function(x) {as.numeric(x>0.2)}
  } else if (mono_form=="step_0.8") {
    mono_f <- function(x) {as.numeric(x>0.8)}
  } else {
    stop("mono_f incorrectly specified")
  }

  # probs <- beta0 + beta1*bmi + beta2*sex + beta*mono_f(antib)
  probs <- expit(beta0 + beta1*bmi + beta2*sex + beta*mono_f(antib))
  infected <- rbinom(n, size=1, prob=probs)
  
  dat <- data.frame(
    bmi = bmi,
    sex = sex,
    antib = antib,
    infected = infected
  )
  
  return (dat)
  
}



#' Hypothesis test based on logistic regression
#' 
#' @param dat Data returned by generate_data_dr()
#' @param params Unused
#' @return Binary; is null rejected (1) or not (0)
test_regression <- function(dat, params) {
  
  model <- glm(
    infected ~ bmi + sex + antib,
    data = dat,
    family = "binomial"
    # family = binomial(link="identity")
  )
  
  reject <- as.numeric(summary(model)$coefficients["antib",4]<0.05)
  
  return (reject)
  
}



#' Hypothesis testing method 2: slope
#' 
#' @param dat Data returned by generate_data_dr()
#' @param params A list; `est` is the estimator used for Theta_hat; one of
#'     c("glm", "sm spline")
#' @return Binary; is null rejected (1) or not (0)
slope_dr <- function(dat, params) {
  
  if (params$est=="glm") {
    
    Theta_hat_constr <- function(dat) {
      model <- glm(
        infected ~ bmi + sex + antib,
        data = dat,
        family = "binomial"
        # family = binomial(link="identity")
      )
      beta <- summary(model)$coefficients["antib",1]
      
      dat <- generate_data_dr(1000,0.6,"square")
      return (function(x) {
        0.5*beta*(x^2)
      })
    }
    
  } else if (params$est=="sm spline") {
    
    Theta_hat_constr <- function(dat) {
      
      model <- gam(
        infected ~ bmi + sex + s(antib, fx=FALSE, bs="cr", m=2, pc=0),
        data = dat,
        family = "binomial"
        # family = binomial(link="identity")
      )
      
      spline_vals <- as.numeric(predict(
        model,
        newdata = list(bmi=rep(25,101), sex=rep(0,101), antib=seq(0,1,0.01)),
        type = "terms"
      )[,3])
      
      return (function(x) {
        0.01 * sum(spline_vals[1:(round(100*round(x,2)+1,0))])
      })
    }
    
  }
  
  Theta_hat <- Vectorize(Theta_hat_constr(dat))
  
  x <- dat$antib
  mu_2n <- mean(x^2)
  mu_3n <- mean(x^3)
  beta_n <- mean((mu_2n*x^2 - mu_3n*x)*Theta_hat(x))
  
  # Define the statistic to bootstrap
  bootstat <- function(dat,indices) {
    d <- dat[indices,]
    Theta_hat <- Vectorize(Theta_hat_constr(d))
    x <- dat$antib
    mu_2n <- mean(x^2)
    mu_3n <- mean(x^3)
    return (mean((mu_2n*x^2 - mu_3n*x)*(Theta_hat(x)-x)))
  }
  
  # Run bootstrap
  boot_obj <- boot(data=dat, statistic=bootstat, R=100) # 1000
  
  # Calculate critical value
  # crit_val <- as.numeric(quantile(boot_obj$t, 0.05))
  # !!!!! This is a one-sided test; either make this two-sided or do the Wald test above as a one-sided test
  crit_val <- qnorm(0.05, mean=mean(boot_obj$t), sd=sd(boot_obj$t))
  
  return(as.numeric(crit_val>0))
  
}



#' Hypothesis test based on G-computation
#' 
#' @param dat Data returned by generate_data_dr()
#' @return Binary; is null rejected (1) or not (0)
test_gcomp <- function(dat) {
  
  model <- glm(infected~bmi+sex+antib, data=dat, family="binomial")
  
  dat %<>% mutate(
    
  )
  
  reject <- as.numeric(summary(model)$coefficients["antib",4]<0.05)
  
  return (reject)
  
}
