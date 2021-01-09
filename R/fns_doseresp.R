#' Generate data according to a density
#' 
#' @param n Sample size
#' @return A dataframe representing a study population
generate_data_dr <- function(n) {
  
  expit <- function(x) {exp(x)/(1+exp(x))}
  
  # Fix parameters
  beta0 <- -4
  beta1 <- 0.05
  beta2 <- 0.10
  beta3 <- 0.03
  
  bmi <- rnorm(n, mean=25, sd=4)
  sex <- rbinom(n, size=1, prob=0.5)
  antib <- sample(c(1:100), size=n, replace=TRUE) # !!!!! Make this a confounder
  probs <- expit(beta0 + beta1*bmi + beta2*sex + beta3*antib) # !!!!! generalize the beta3*antib bit
  # !!!!! Use identity link instead?
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
#' @return Binary; is null rejected (1) or not (0)
test_logistic <- function(dat) {
  
  model <- glm(infected~bmi+sex+antib, data=dat, family="binomial")
  reject <- as.numeric(summary(model)$coefficients["antib",4]<0.05)
  
  return (reject)
  
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
