###################.
##### Density #####
###################.

if(cfg$setting=="density") {
  
  #' Generate data according to a density
  #' 
  #' @param n Sample size
  #' @param true_density A character string representing the true density
  #' @return A vector; a sample from the specified density
  #' @notes For "Linear" and "f(x)=ke^x", inverse-CDF sampling is used. For
  #'     "Exponential" and "Spline", we approximate the continuous distribution
  #'     with a discrete distribution.
  generate_data <- function(n, true_density) {
    
    # Q is the quantile function for probability integral transform sampling
    # !!!!! Get rid of density alias
    if (true_density %in% c("f(x)=1", "Uniform")) {
      dat <- runif(n)
    }
    
    if (true_density %in% c("f(x)=2x", "Linear")) {
      Q <- function(x) { sqrt(x) }
      dat <- Q(runif(n))
    }
    
    if (true_density=="Step") {
      dat <- runif(n=n, min=0.5, max=1) - 0.5*rbinom(n=n, size=1, prob=0.05)
    }
    
    if (true_density=="Step (decr)") {
      dat <- runif(n=n, min=0, max=0.5) + 0.5*rbinom(n=n, size=1, prob=0.05)
    }
    
    if (true_density=="f(x)=ke^x") {
      Q <- function(x) { log(x*(exp(1)-1)+1) }
      dat <- Q(runif(n))
    }
    
    if (true_density %in% c("Exponential", "Exponential (decr)", "Spline")) {
      
      # Define a grid over the support
      grid <- seq(0, 1, 0.001)
      
      # Define the pdf
      if (true_density=="Exponential") {
        pdf <- function(x){
          c <- 5
          c*(exp(c)-1-c)^-1 * (exp(c*x)-1)
        }
      }
      if (true_density=="Exponential (decr)") {
        pdf <- function(x){
          c <- 5
          c*(exp(c)-1-c)^-1 * (exp(c*(1-x))-1)
        }
      }
      if (true_density=="Spline") {
        pdf <- function(x){
          (x<=0.5) * 0.1 +
            (0.5<x)*(x<=0.6) * (20*x-9.9) +
            (x>0.6) * 2.1
        }
      }
      
      # Compute probabilities over grid
      probs <- pdf(grid)
      probs <- probs/sum(probs)
      
      # Draw sample
      dat <- sample(x=grid, size=n, replace=TRUE, prob=probs)
      
    }
    
    return(dat)
    
  }
  
}



######################.
##### Regression #####
######################.

if(cfg$setting=="regression") {
  
  #' Generate data
  #' 
  #' @param n Sample size
  #' @param alpha_3 Height of the dose-response curve
  #' @param sigma Standard deviation of normal error
  #' @param mono_form Functional form of the regression curve; should be a
  #'     monotone increasing function with domain [0,1] such that f(0)=0 and
  #'     f(1)=1; choices are c("identity", "square", "sqrt", "step_0.2",
  #'     "step_0.8")
  #' @param a_distr Beta parameters for marginaldistribution of A; a list of the
  #'     form list(shape1=3, shape2=4)
  #' @return A dataframe representing a study population
  generate_data <- function(n, alpha_3, sigma, mono_form, a_distr) {
    
    # Sample baseline covariates
    a <- rbeta(n, shape1=a_distr$shape1, shape2=a_distr$shape2)
    
    # Set transformation function
    if (mono_form=="identity") {
      mono_f <- function(x) {x}
    } else if (mono_form=="square") {
      mono_f <- function(x) {x^2}
    } else if (mono_form=="sqrt") {
      mono_f <- function(x) {sqrt(x)}
    } else if (mono_form=="step_0.2") {
      mono_f <- function(x) {as.numeric(x>0.2)}
    } else if (mono_form=="step_0.5") {
      mono_f <- function(x) {as.numeric(x>0.5)}
    } else if (mono_form=="step_0.8") {
      mono_f <- function(x) {as.numeric(x>0.8)}
    } else {
      stop("mono_form incorrectly specified")
    }
    
    y <- alpha_3*mono_f(a) + rnorm(n, mean=0, sd=sigma)
    
    dat <- data.frame(a=a, y=y)
    
    return (dat)
    
  }
  
}



#########################.
##### Dose Response #####
#########################.

if(cfg$setting=="doseresp") {
  
  #' Generate data
  #' 
  #' @param n Sample size
  #' @param alpha_3 Height of the dose-response curve
  #' @param mono_form Functional form of dose-response curve; should be a monotone
  #'     increasing function with domain [0,1] such that f(0)=0 and f(1)=1;
  #'     choices are c("identity", "square", "sqrt", "step_0.2", "step_0.8")
  #' @param sampling A list. One of c("iid","two-phase")
  #' @return A dataframe representing the study population
  generate_data <- function(n, alpha_3, mono_form, sampling="iid") {
    
    # Fix parameters
    alpha_0 <- -1.5
    alpha_1 <- 0.3
    alpha_2 <- 0.7
    
    # Sample baseline covariates
    w1 <- rnorm(n)
    w2 <- rbinom(n, size=1, prob=0.5)
    shape_2 <- ifelse(w2==1, 1.5, 1.1)
    a <- rbeta(n, shape1=0.9, shape2=shape_2)
    
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
      stop("mono_form incorrectly specified")
    }
    
    # Sample outcome
    probs <- expit(alpha_0 + alpha_1*w1 + alpha_2*w2 + alpha_3*mono_f(a))
    y <- rbinom(n, size=1, prob=probs)
    
    # IID sampling
    if (sampling=="iid") {
      return (data.frame(w1=w1, w2=w2, a=a, y=y))
    }
    
    # Two-phase sampling
    if (sampling=="two-phase") {
      pi <- expit(2*y-w2)
      delta <- rbinom(n, size=1, prob=pi)
      return (data.frame(w1=w1, w2=w2, a=ifelse(delta==1,a,NA), y=y))
    }
    
  }
  
}
