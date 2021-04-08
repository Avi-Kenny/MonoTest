###################.
##### Density #####
###################.

#' Testing approach 1: cumulative increments (all functions have this structure)
#' 
#' @param dat Data returned by generate_data()
#' @param alt_type Type of alternative hypothesis; either "incr" or "decr"
#' @param params A list, specific to the test
#' @return Binary; is null rejected (1) or not (0)


if(cfg$setting=="density") {
  
  #' @param params A list, containing the following:
  #'   - `delta` is the interval size (e.g. 1/3).
  #'   - `wts` can either be "equal" or a a vector of length (1/delta)-1
  #' @notes
  #'   - Currently assumes equally sized intervals of size delta
  test1 <- function(dat, alt_type="incr", params) {
    
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
  
}



######################.
##### Regression #####
######################.

if(cfg$setting=="regression") {
  
  #' @param params A list, containing the following:
  #'   - `p1` text
  #'   - `p2` text
  #' @notes
  #'   - Note
  test1 <- function(dat, alt_type="incr", params) {
    #
  }
  
}



#########################.
##### Dose Response #####
#########################.

if(cfg$setting=="doseresp") {
  
  #' @param params A list, containing the following:
  #'   - `p1` text
  #'   - `p2` text
  #' @notes
  #'   - Note
  test1 <- function(dat, alt_type="incr", params) {
    #
  }
  
}
