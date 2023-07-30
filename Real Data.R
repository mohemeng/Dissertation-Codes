EM.case.3.alpha.sigma <- function(data, tol = 1e-2, maxiter = 500, sigma = 1, alpha = 1, verb = FALSE){
  #------------------------- First space (checking input errors) ----------------------------------------#
  if(!is.numeric(sigma) | !is.numeric(alpha) | !is.numeric(tol) | !is.numeric(data)) stop(
    "arguments 'y', 'tol', 'maxiter', 'sigma', 'alpha' must be numeric")
  if(sigma <= 0 | alpha <= 0) stop("alpha and sigma must be positive")
  #######################################################################################################
  #---------------------------- Second space (necessary functions) --------------------------------------#
  log_like <- function(iter.alpha, iter.sigma){
    sum(log(NP.pdf(data, alpha = iter.alpha, sigma = iter.sigma)))}

  # function to compute all expectations
  all.expectations <- function(iter.sigma, iter.alpha){
    a_k <- array(NA, dim = length(data)); b_k <- array(NA, dim = length(data)); c_k <- array(NA, dim = length(data))
    for(y_i in data){
      fun.1 <- function(t) (sqrt(t/2)/iter.sigma)*exp(-sqrt(2*t)*abs(y_i)/iter.sigma)*(t^(1/iter.alpha - 1)*exp(-t/iter.alpha))/((iter.alpha)^(1/iter.alpha)*gamma(1/iter.alpha))
      fun.2 <- function(t) t*(sqrt(t/2)/iter.sigma)*exp(-sqrt(2*t)*abs(y_i)/iter.sigma)*(t^(1/iter.alpha - 1)*exp(-t/iter.alpha))/((iter.alpha)^(1/iter.alpha)*gamma(1/iter.alpha))
      fun.3 <- function(t) sqrt(t)*(sqrt(t/2)/iter.sigma)*exp(-sqrt(2*t)*abs(y_i)/iter.sigma)*(t^(1/iter.alpha - 1)*exp(-t/iter.alpha))/((iter.alpha)^(1/iter.alpha)*gamma(1/iter.alpha))
      fun.4 <- function(t) log(t)*(sqrt(t/2)/iter.sigma)*exp(-sqrt(2*t)*abs(y_i)/iter.sigma)*(t^(1/iter.alpha - 1)*exp(-t/iter.alpha))/((iter.alpha)^(1/iter.alpha)*gamma(1/iter.alpha))

      dens.y <- adaptIntegrate(fun.1, lowerLimit = 0, upperLimit = Inf)$integral
      a.exp <- adaptIntegrate(fun.2, lowerLimit = 0, upperLimit = Inf)$integral/dens.y
      b.exp <- adaptIntegrate(fun.3, lowerLimit = 0, upperLimit = Inf)$integral/dens.y
      c.exp <- adaptIntegrate(fun.4, lowerLimit = 0, upperLimit = Inf)$integral/dens.y
      pos <- match(y_i, data)
      a_k[pos] <- a.exp; b_k[pos] <- b.exp; c_k[pos] <- c.exp
    }
    dat <- data.frame(a_k, b_k, c_k)
    return(dat)}
  #####################################################
  #### E step
  expect <- all.expectations(iter.alpha = alpha, iter.sigma = sigma)
  a_old <- expect$a_k; b_old <- expect$b_k; c_old <- expect$c_k
  ll_old <- log_like(iter.alpha = alpha, iter.sigma = sigma)
  k = 0
  alpha_old <- alpha; sigma_old <- sigma
  output <- c(k, alpha_old, sigma_old, ll_old)
  diff <- tol + 1

  #t <- try((diff > tol) & (k < maxiter))

  #while(k != "try-error" & t == TRUE){
  while( (diff > tol) & (k < maxiter)){
    #### M step
    alpha <- suppressWarnings(stats::nlm(function(x) length(data)*(log(x)/x + lgamma(1/x)) - sum(c_old - a_old, na.rm = TRUE)/x , p = alpha)$estimate)
    sigma <- sqrt(2)*mean(b_old*abs(data), na.rm = TRUE)
    expect <- all.expectations(iter.alpha = alpha, iter.sigma = sigma)
    a_old <- expect$a_k; b_old <- expect$b_k; c_old <- expect$c_k
    ll_new <- log_like(iter.alpha = alpha, iter.sigma = sigma)
    diff <- ll_new - ll_old
    ll_old <- ll_new
    k <- k + 1
    t <- try((diff > tol) & (k < maxiter))
    output <- rbind(output, c(k, alpha, sigma, ll_new))
    if(verb){
      cat("iteration =", k, " log.lik.diff =", diff, " log.lik =",
          ll_new, "\n")}

    if(k == maxiter){
      cat("Warning! Convergence not achieved!", "\n")
    }
    output <- data.frame(output)
    colnames(output) <- c("iteration","alpha", "sigma","log-lik")
    par <- data.frame(t(c(alpha, sigma)))
    colnames(par) <- c("alpha", "sigma")
    result <- list(par=par, log.like = ll_new, iterations=k, output=output)
  }
  return(result)
}


#### Real Data Fit
SGD <- read.csv('SGD.csv')
SGD.return <- as.numeric(SGD$Close)
SGD.return <- SGD.return[!is.na(SGD.return)]
SGD.return <- diff(log(SGD.return))
SGD.return <- SGD.return[which(SGD.return!= 0) ]
hist(SGD.return, freq = FALSE, breaks = 50, main = "", xlab = "USD/SGD return")

MME.case3.alpha.sigma(data = SGD.return)

EM.case.3.alpha.sigma(data = SGD.return, sigma = 0.5, alpha = 0.5)
