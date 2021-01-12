#' Expit function
#' 
#' @param x Numeric input
#' @return Numeric output
expit <- function(x) {1/(1+exp(-x))}



#' Integral of expit function
#' 
#' @param x Numeric input
#' @return Numeric output
intexpit <- function(x) {log(1+exp(x))}



#' Generate data according to a density
#' 
#' @param n Sample size
#' @param beta Height of the dose-response curve
#' @param mono_form Functional form of dose-response curve; should be a monotone
#'     increasing function with domain [0,1] such that f(0)=0 and f(1)=1;
#'     choices are c("identity", "square", "sqrt", "step_0.2", "step_0.8")
#' @return A dataframe representing a study population
generate_data_dr <- function(n, beta, mono_form) {
  
  # Fix parameters
  beta0 <- -3
  beta1 <- 0.03
  beta2 <- 0.7
  
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
  
  model <- glm(infected~bmi+sex+antib, data=dat, family="binomial")
  one_sided_p <- pnorm(summary(model)$coefficients["antib",3], lower.tail=F)
  reject <- as.numeric(one_sided_p<0.05)
  
  return (reject)
  
}



#' Hypothesis testing method 2: slope
#' 
#' @param dat Data returned by generate_data_dr()
#' @param params A list; `est` is the estimator used for Theta_hat; one of
#'     c("glm", "sm spline")
#' @return Binary; is null rejected (1) or not (0)
slope_dr <- function(dat, params) {
  
  Theta_hat_constr <- function(dat, subtype) {
    
    if (subtype=="glm") {
      
      # Run model and extract coefficients
      model <- glm(infected~bmi+sex+antib, data=dat, family="binomial")
      coeff <- summary(model)$coefficients
      beta0_hat <- coeff["(Intercept)",1]
      beta1_hat <- coeff["bmi",1]
      beta2_hat <- coeff["sex",1]
      beta3_hat <- coeff["antib",1]
      
      # Create Theta_hat function
      Theta_hat <- Vectorize(function(x) {
        
        t_i <- apply(
          X = dat,
          MARGIN = 1,
          FUN = function(r) {
            log( (1+exp(beta0_hat + beta1_hat*r[["bmi"]] +
                  beta2_hat*r[["sex"]] + beta3_hat*x)) /
                (1+exp(beta0_hat + beta1_hat*r[["bmi"]] + beta2_hat*r[["sex"]]))
            )
          }
        )
        
        return ((beta3_hat^-1)*mean(t_i))
        
      })
      
      return (Theta_hat)
      
    }
    
    if (subtype=="ss") {
      
      # Run model and extract coefficients
      model <- gam(
        infected ~ bmi + sex + s(antib, fx=FALSE, bs="cr", m=2, pc=0),
        data = dat,
        family = "binomial"
      )
      
      coeff <- model$coefficients
      beta0_hat <- coeff[["(Intercept)"]]
      beta1_hat <- coeff[["bmi"]]
      beta2_hat <- coeff[["sex"]]
      spline_vals <- as.numeric(predict(
        model,
        newdata = list(bmi=rep(25,101), sex=rep(0,101), antib=seq(0,1,0.01)),
        type = "terms"
      )[,3])
      
      # Construct theta_hat from GAM formula
      theta_hat <- Vectorize(function(x) {
        E_hat_i <- apply(
          X = dat,
          MARGIN = 1,
          FUN = function(r) {
            expit(beta0_hat + beta1_hat*r[["bmi"]] + beta2_hat*r[["sex"]] +
                  spline_vals[1:(round(100*round(x,2)+1,0))])
          }
        )
        return (mean(E_hat_i))
      })
      
      # Construct Theta_hat by integrating theta_hat
      # Approximating integral using a Riemann sum with ten intervals
      Theta_hat <- Vectorize(
        function (a) { a * mean(theta_hat(seq(a/10,a,a/10))) }
      )
      
      return (Theta_hat)
      
    }
    
  }
  
  # Construct Theta_hat function
  Theta_hat <- Theta_hat_constr(dat, params$subtype)

  # Calculate value of test statistic
  x <- dat$antib
  mu_2n <- mean(x^2)
  mu_3n <- mean(x^3)
  beta_n <- mean((mu_2n*x^2 - mu_3n*x)*Theta_hat(x))
  
  # Define the statistic to bootstrap
  bootstat <- function(dat,indices) {
    d <- dat[indices,]
    Theta_hat <- Theta_hat_constr(d, params$subtype)
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
