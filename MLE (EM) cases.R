library(cubature)
########## Normal Pareto Generator ###################
rNPare <- function(n = 1, alpha = 0.5, sigma = 1){
  normal <- rnorm(n); exponential <- rexp(n)
  gamm <- rgamma(n, shape = 1/alpha, scale = alpha)
  y <- sigma*sqrt(exponential/gamm)*normal
  return(y)
}

########## PDF for Normal Pareto #####################
NP.pdf <- function(y, alpha = 0.5, sigma = 1){
  bb <- array(NA, dim = length(y))
  for(i in y){
    pos <- match(i, y)
    gg <- function(t){
      uu <- (sqrt(t/2)/abs(sigma))*exp(-sqrt(2*t)*abs(i)/abs(sigma))*(t^(1/abs(alpha) - 1)*exp(-t/abs(alpha)))/(abs(alpha)^(1/abs(alpha))*gamma(1/abs(alpha)))
      return(uu)
    }
    bb[pos] <- adaptIntegrate(gg, lowerLimit = 0, upperLimit = Inf)$integral
  }
  return(bb)}

########### MLE estimation using EM Algorithm ########
# Case 1 sigma unknown, alpha known
EM.case.1.sigma <- function(data, tol = 1e-2, maxiter = 500, sigma = 1, alpha = 1, verb = FALSE){
  #------------------------- First space (checking input errors) ----------------------------------------#
  if(!is.numeric(sigma) | !is.numeric(tol) | !is.numeric(data)) stop(
    "arguments 'y', 'tol', 'maxiter', 'sigma' must be numeric")
  if(sigma <= 0) stop("alpha and sigma must be positive")
  #######################################################################################################
  #---------------------------- Second space (necessary functions) --------------------------------------#
  log_like <- function(iter.alpha = alpha, iter.sigma){
    sum(log(NP.pdf(data, alpha = iter.alpha, sigma = iter.sigma)))}

  # function to compute all expectations
  all.expectations <- function(iter.sigma, iter.alpha = alpha){
    b_k <- array(NA, dim = length(data))
    for(y_i in data){
      fun.1 <- function(t) (sqrt(t/2)/iter.sigma)*exp(-sqrt(2*t)*abs(y_i)/iter.sigma)*(t^(1/iter.alpha - 1)*exp(-t/iter.alpha))/((iter.alpha)^(1/iter.alpha)*gamma(1/iter.alpha))
      fun.3 <- function(t) sqrt(t)*(sqrt(t/2)/iter.sigma)*exp(-sqrt(2*t)*abs(y_i)/iter.sigma)*(t^(1/iter.alpha - 1)*exp(-t/iter.alpha))/((iter.alpha)^(1/iter.alpha)*gamma(1/iter.alpha))

      dens.y <- adaptIntegrate(fun.1, lowerLimit = 0, upperLimit = Inf)$integral
      b.exp <- adaptIntegrate(fun.3, lowerLimit = 0, upperLimit = Inf)$integral/dens.y

      pos <- match(y_i, data)
      b_k[pos] <- b.exp
    }
    dat <- data.frame(b_k)
    return(dat)}
  #####################################################
  #### E step
  expect <- all.expectations(iter.sigma = sigma)
  b_old <- expect$b_k
  ll_old <- log_like(iter.sigma = sigma)
  k = 0
  sigma_old <- sigma
  output <- c(k, sigma_old, ll_old)
  diff <- tol + 1

  while((diff > tol) & (k < maxiter)){
    #### M step
    sigma <- sqrt(2)*mean(b_old*abs(data))
    expect <- all.expectations(iter.sigma = sigma)
    b_old <- expect$b_k
    ll_new <- log_like(iter.sigma = sigma)
    diff <- ll_new - ll_old
    ll_old <- ll_new
    k <- k + 1
    output <- rbind(output, c(k, sigma, ll_new))
    if(verb){
      cat("iteration =", k, " log.lik.diff =", diff, " log.lik =",
          ll_new, "\n")}

    if(k == maxiter){
      cat("Warning! Convergence not achieved!", "\n")
    }
    output <- data.frame(output)
    colnames(output) <- c("iteration","sigma","log-lik")
    par <- data.frame(t(c(sigma)))
    colnames(par) <- c("sigma")
    result <- list(par=par, log.like = ll_new, iterations=k, output=output)
  }
  return(result)
}

# Iteration plots
rt1 <- rNPare(n = 1000, alpha = 1, sigma = 0.5)
rt2 <- rNPare(n = 1000, alpha = 1, sigma = 1)
rt3 <- rNPare(n = 1000, alpha = 1, sigma = 5)

l1 <- EM.case.1.sigma(data = rt1, alpha = 1, sigma = 0.2)
l2 <- EM.case.1.sigma(data = rt1, alpha = 1, sigma = MME.case1.sigma(alpha = 1, data = rt1))
l3 <- EM.case.1.sigma(data = rt1, alpha = 1, sigma = 1.5 )

plot(l1$output$iteration, l1$output$sigma, ylim = c(0, max(l1$output$sigma,l2$output$sigma,l3$output$sigma)), type = "o", pch = 19,
     xlab = "Iteration", ylab = expression(paste(sigma, " values")))
points(l2$output$iteration, l2$output$sigma, type = "o", pch = 19, col = "blue")
points(l3$output$iteration, l3$output$sigma, type = "o", col = "green", pch = 19)
abline(h = 0.5, col = "red", lty = 2)
legend(4, 1.2, legend = c("Intial values", expression(paste(sigma, " = 0.3075 (MME)" )),
                          expression(paste(sigma, " = 1.5" )), expression(paste(sigma, " = 0.2"))),
       col = c("", "blue", "green", "black"), lty = c(NA, 1, 1,1), pch = c(NA, 19, 19, 19), bty = "n",
       cex = 1.2, seg.len = 1, x.intersp = 0.1, y.intersp = 0.8)

m1 <- EM.case.1.sigma(data = rt2, alpha = 1, sigma = 0.2, tol = 1e-4)
m2 <- EM.case.1.sigma(data = rt2, alpha = 1, sigma = MME.case1.sigma(alpha = 1, data = rt2), tol = 1e-4)
m3 <- EM.case.1.sigma(data = rt2, alpha = 1, sigma = 1.5, tol = 1e-4)

plot(m1$output$iteration, m1$output$sigma, ylim = c(0, max(m1$output$sigma,m2$output$sigma,m3$output$sigma)), type = "o", pch = 19,
     xlab = "Iteration", ylab = expression(paste(sigma, " values")))
points(m2$output$iteration, m2$output$sigma, type = "o", pch = 19, col = "blue")
points(m3$output$iteration, m3$output$sigma, type = "o", col = "green", pch = 19)
abline(h = 1, col = "red", lty = 2)
legend(6, 0.7, legend = c("Intial values", expression(paste(sigma, " = 1.2481 (MME)" )),
                          expression(paste(sigma, " = 1.5" )), expression(paste(sigma, " = 0.2"))),
       col = c("", "blue", "green", "black"), lty = c(NA, 1, 1,1), pch = c(NA, 19, 19, 19), bty = "n",
       cex = 1.2, seg.len = 1, x.intersp = 0.1, y.intersp = 0.8)

n1 <- EM.case.1.sigma(data = rt3, alpha = 1, sigma = 0.2, tol = 1e-4)
n2 <- EM.case.1.sigma(data = rt3, alpha = 1, sigma = MME.case1.sigma(alpha = 1, data = rt3), tol = 1e-4)
n3 <- EM.case.1.sigma(data = rt3, alpha = 1, sigma = 10, tol = 1e-4)

plot(n1$output$iteration, n1$output$sigma, ylim = c(0, max(n1$output$sigma,n2$output$sigma,n3$output$sigma)), type = "o", pch = 19,
     xlab = "Iteration", ylab = expression(paste(sigma, " values")))
points(n2$output$iteration, n2$output$sigma, type = "o", pch = 19, col = "blue")
points(n3$output$iteration, n3$output$sigma, type = "o", col = "green", pch = 19)
abline(h = 5, col = "red", lty = 2)
legend(6, 8, legend = c("Intial values", expression(paste(sigma, " = 3.3465 (MME)")),
                          expression(paste(sigma, " = 10" )), expression(paste(sigma, " = 0.2"))),
       col = c("", "blue", "green", "black"), lty = c(NA, 1, 1,1), pch = c(NA, 19, 19, 19), bty = "n",
       cex = 1.2, seg.len = 1, x.intersp = 0.1, y.intersp = 0.8)

sim.MLE.case1.sigma <- function(k = 1000, alpha = 0.5, sigma = 1, num = 10){
  num.vec <- matrix(NA, nrow = k, ncol = 5)
  colnames(num.vec) <- c("Sample", "alpha", "sigma", "Estimate", "RelativeError")
  num.vec[ , 1] <- rep(num, k)
  num.vec[ , 2] <- rep(alpha, k)
  num.vec[ , 3] <- rep(sigma, k)

  # Generate independent Y's
  for(i in 1:k){
    cc <- rNPare(n = num, alpha = alpha, sigma = sigma)
    num.vec[i, 4] <- as.numeric(EM.case.1.sigma(alpha = alpha, data = cc, sigma = MME.case1.sigma(alpha = alpha, data = cc))$par)
    num.vec[i, 5] <- (num.vec[i, 4] - num.vec[i, 3])/num.vec[i, 3]
  }
  return(data.frame(num.vec))
}

l <- rbind(sim.MLE.case1.sigma(alpha = 0.7, sigma = 0.5, num = 50), sim.MLE.case1.sigma(alpha = 0.7, sigma = 0.5, num = 150), sim.MLE.case1.sigma(alpha = 0.7, sigma = 0.5, num = 300),
           sim.MLE.case1.sigma(alpha = 1, sigma = 1, num = 50), sim.MLE.case1.sigma(alpha = 1, sigma = 1, num = 150), sim.MLE.case1.sigma(alpha = 1, sigma = 1, num = 300),
           sim.MLE.case1.sigma(alpha = 2, sigma = 5, num = 50), sim.MLE.case1.sigma(alpha = 2, sigma = 5, num = 150), sim.MLE.case1.sigma(alpha = 2, sigma = 5, num = 300),
           sim.MLE.case1.sigma(alpha = 1, sigma = 0.5, num = 50), sim.MLE.case1.sigma(alpha = 2, sigma = 0.5, num = 50), sim.MLE.case1.sigma(alpha = 1, sigma = 0.5, num = 150),
           sim.MLE.case1.sigma(alpha = 2, sigma = 0.5, num = 150), sim.MLE.case1.sigma(alpha = 1, sigma = 0.5, num = 300), sim.MLE.case1.sigma(alpha = 2, sigma = 0.5, num = 300),
           sim.MLE.case1.sigma(alpha = 0.7, sigma = 1, num = 50), sim.MLE.case1.sigma(alpha = 2, sigma = 1, num = 50), sim.MLE.case1.sigma(alpha = 0.7, sigma = 1, num = 150),
           sim.MLE.case1.sigma(alpha = 2, sigma = 1, num = 150), sim.MLE.case1.sigma(alpha = 0.7, sigma = 1, num = 300), sim.MLE.case1.sigma(alpha = 2, sigma = 1, num = 300),
           sim.MLE.case1.sigma(alpha = 0.7, sigma = 5, num = 50), sim.MLE.case1.sigma(alpha = 1, sigma = 5, num = 50), sim.MLE.case1.sigma(alpha = 0.7, sigma = 5, num = 150),
           sim.MLE.case1.sigma(alpha = 1, sigma = 5, num = 150), sim.MLE.case1.sigma(alpha = 0.7, sigma = 5, num = 300), sim.MLE.case1.sigma(alpha = 1, sigma = 5, num = 300))

load(file = "C:/Users/Esther/Desktop/Normal Pareto/MLE.case1.sigma.1000.simulations1.RData")
load(file = "C:/Users/Esther/Desktop/Normal Pareto/MLE.case1.sigma.1000.simulations2.RData")
load(file = "C:/Users/Esther/Desktop/Normal Pareto/MLE.case1.sigma.1000.simulations3.RData")

save(l, file = "MLE.case1.sigma.1000.simulations.RData")
l1 <- l
l <- rbind(l1, l2, l3)

library(ggplot2)
# Relative error graphs
gh <- ggplot(l, aes(x = RelativeError)) + geom_histogram(aes(y = ..density..)) + facet_grid(sigma ~ Sample + alpha, scales = "free", labeller = label_both) + theme_classic()
gh <- gh + labs(x = "Relative Error")
gh

# Boxplots
th <- ggplot(l, aes(x = as.factor(Sample), y = Estimate)) + geom_boxplot() + geom_hline(aes(yintercept = sigma), color = "red") +
  facet_wrap(sigma ~ alpha, scales = "free", labeller = label_both) + labs(x = "Sample size") + theme_classic()
th

# Generating table
tab1 <- l[which(l$Sample == 300 & l$sigma == 5 & l$alpha == 2), ]
(esti <- round(mean(tab1$Estimate), 4))
(MSE <- round(mean( (tab1$sigma - tab1$Estimate)^2 ), 4 ))

# Case 2 alpha unknown, sigma known
EM.case.2.alpha <- function(data, tol = 1e-2, maxiter = 500, sigma = 1, alpha = 1, verb = FALSE){
  #------------------------- First space (checking input errors) ----------------------------------------#
  if(!is.numeric(sigma) | !is.numeric(tol) | !is.numeric(data)) stop(
    "arguments 'y', 'tol', 'maxiter', 'sigma' must be numeric")
  if(sigma <= 0 | alpha <= 0) stop("alpha and sigma must be positive")
  #######################################################################################################
  #---------------------------- Second space (necessary functions) --------------------------------------#
  log_like <- function(iter.alpha, iter.sigma = sigma){
    sum(log(NP.pdf(data, alpha = iter.alpha, sigma = iter.sigma)))}

  # function to compute all expectations
  all.expectations <- function(iter.sigma = sigma, iter.alpha){
    a_k <- array(NA, dim = length(data)); c_k <- array(NA, dim = length(data))
    for(y_i in data){
      fun.1 <- function(t) (sqrt(t/2)/iter.sigma)*exp(-sqrt(2*t)*abs(y_i)/iter.sigma)*(t^(1/iter.alpha - 1)*exp(-t/iter.alpha))/((iter.alpha)^(1/iter.alpha)*gamma(1/iter.alpha))
      fun.2 <- function(t) t*(sqrt(t/2)/iter.sigma)*exp(-sqrt(2*t)*abs(y_i)/iter.sigma)*(t^(1/iter.alpha - 1)*exp(-t/iter.alpha))/((iter.alpha)^(1/iter.alpha)*gamma(1/iter.alpha))
      fun.4 <- function(t) log(t)*(sqrt(t/2)/iter.sigma)*exp(-sqrt(2*t)*abs(y_i)/iter.sigma)*(t^(1/iter.alpha - 1)*exp(-t/iter.alpha))/((iter.alpha)^(1/iter.alpha)*gamma(1/iter.alpha))

      dens.y <- adaptIntegrate(fun.1, lowerLimit = 0, upperLimit = Inf)$integral
      a.exp <- adaptIntegrate(fun.2, lowerLimit = 0, upperLimit = Inf)$integral/dens.y
      c.exp <- adaptIntegrate(fun.4, lowerLimit = 0, upperLimit = Inf)$integral/dens.y
      pos <- match(y_i, data)
      a_k[pos] <- a.exp; c_k[pos] <- c.exp
    }
    dat <- data.frame(a_k, c_k)
    return(dat)}
  #####################################################
  #### E step
  expect <- all.expectations(iter.alpha = alpha)
  a_old <- expect$a_k; c_old <- expect$c_k
  ll_old <- log_like(iter.alpha = alpha)
  k = 0
  alpha_old <- alpha
  output <- c(k, alpha_old, ll_old)
  diff <- tol + 1

  while((diff > tol) & (k < maxiter)){
    #### M step
    alpha <- suppressWarnings(stats::nlm(function(x) length(data)*(log(x)/x + lgamma(1/x)) - sum(c_old - a_old)/x , p = alpha)$estimate)
    expect <- all.expectations(iter.alpha = alpha)
    a_old <- expect$a_k; c_old <- expect$c_k
    ll_new <- log_like(iter.alpha = alpha)
    diff <- ll_new - ll_old
    ll_old <- ll_new
    k <- k + 1
    output <- rbind(output, c(k, alpha, ll_new))
    if(verb){
      cat("iteration =", k, " log.lik.diff =", diff, " log.lik =",
          ll_new, "\n")}

    if(k == maxiter){
      cat("Warning! Convergence not achieved!", "\n")
    }
    output <- data.frame(output)
    colnames(output) <- c("iteration","alpha","log-lik")
    par <- data.frame(t(c(alpha)))
    colnames(par) <- c("alpha")
    result <- list(par=par, log.like = ll_new, iterations=k, output=output)
  }
  return(result)
}

# Iteration plots
rt1 <- rNPare(n = 1000, alpha = 0.7, sigma = 1)
rt2 <- rNPare(n = 1000, alpha = 1, sigma = 1)
rt3 <- rNPare(n = 1000, alpha = 2, sigma = 1)

l1 <- EM.case.2.alpha(data = rt1, alpha = 1.5, sigma = 1)
l2 <- EM.case.2.alpha(data = rt1, alpha = MME.case2.alpha(sigma = 1, data = rt1), sigma = 1)
l3 <- EM.case.2.alpha(data = rt1, alpha = 0.3, sigma = 1 )

plot(l1$output$iteration, l1$output$alpha, ylim = c(0, max(l1$output$alpha,l2$output$alpha,l3$output$alpha)), type = "o", pch = 19,
     xlab = "Iteration", ylab = expression(paste(alpha, " values")), xlim = c(0, max(l1$output$iteration, l2$output$iteration, l3$output$iteration)))
points(l2$output$iteration, l2$output$alpha, type = "o", pch = 19, col = "blue")
points(l3$output$iteration, l3$output$alpha, type = "o", col = "green", pch = 19)
abline(h = 0.7, col = "red", lty = 2)
legend(15, 1.2, legend = c("Intial values", expression(paste(alpha, " = 1.0559 (MME)" )),
                          expression(paste(alpha, " = 1.5" )), expression(paste(alpha, " = 0.3"))),
       col = c("", "blue", "green", "black"), lty = c(NA, 1, 1,1), pch = c(NA, 19, 19, 19), bty = "n",
       cex = 1.2, seg.len = 1, x.intersp = 0.1, y.intersp = 0.8)

m1 <- EM.case.2.alpha(data = rt2, alpha = 0.5, sigma = 1)
m2 <- EM.case.2.alpha(data = rt2, alpha = MME.case2.alpha(sigma = 1, data = rt2), sigma = 1)
m3 <- EM.case.2.alpha(data = rt2, alpha = 1.5, sigma = 1)

plot(m1$output$iteration, m1$output$alpha, ylim = c(0, max(m1$output$alpha,m2$output$alpha,m3$output$alpha)), type = "o", pch = 19,
     xlab = "Iteration", ylab = expression(paste(alpha, " values")), xlim = c(0, max(m1$output$iteration, m2$output$iteration, m3$output$iteration)))
points(m2$output$iteration, m2$output$alpha, type = "o", pch = 19, col = "blue")
points(m3$output$iteration, m3$output$alpha, type = "o", col = "green", pch = 19)
abline(h = 1, col = "red", lty = 2)
legend(8, 0.5, legend = c("Intial values", expression(paste(alpha, " = 0.8969 (MME)" )),
                          expression(paste(alpha, " = 1.5" )), expression(paste(alpha, " = 0.5"))),
       col = c("", "blue", "green", "black"), lty = c(NA, 1, 1,1), pch = c(NA, 19, 19, 19), bty = "n",
       cex = 1.2, seg.len = 1, x.intersp = 0.1, y.intersp = 0.8)

n1 <- EM.case.2.alpha(data = rt3, alpha = 0.7, sigma = 1)
n2 <- EM.case.2.alpha(data = rt3, alpha = MME.case2.alpha(sigma = 1, data = rt3), sigma = 1)
n3 <- EM.case.2.alpha(data = rt3, alpha = 5, sigma = 1)

plot(n1$output$iteration, n1$output$alpha, ylim = c(0, max(n1$output$alpha,n2$output$alpha,n3$output$alpha)), type = "o", pch = 19,
     xlab = "Iteration", ylab = expression(paste(alpha, " values")), xlim = c(0, max(n1$output$iteration, n2$output$iteration, n3$output$iteration)))
points(n2$output$iteration, n2$output$alpha, type = "o", pch = 19, col = "blue")
points(n3$output$iteration, n3$output$alpha, type = "o", col = "green", pch = 19)
abline(h = 2, col = "red", lty = 2)
legend(5, 4, legend = c("Intial values", expression(paste(alpha, " = 1.7960 (MME)")),
                        expression(paste(alpha, " = 5" )), expression(paste(alpha, " = 0.7"))),
       col = c("", "blue", "green", "black"), lty = c(NA, 1, 1,1), pch = c(NA, 19, 19, 19), bty = "n",
       cex = 1.2, seg.len = 1, x.intersp = 0.1, y.intersp = 0.8)

first <- list(rt1,l1,l2,l3)
second <- list(rt2,m1,m2,m3)
third <-list(rt3,n1,n2,n3)

save(first, file = "EM.data.1000.alpha.unknown.0.7.sigma.known.1.RData")
save(second, file = "EM.data.1000.alpha.unknown.1.sigma.known.1.RData")
save(third, file = "EM.data.1000.alpha.unknown.2.sigma.known.1.RData")

sim.MLE.case2.alpha <- function(k = 1000, alpha = 0.5, sigma = 1, num = 10){
  num.vec <- matrix(NA, nrow = k, ncol = 5)
  colnames(num.vec) <- c("Sample", "alpha", "sigma", "Estimate", "RelativeError")
  num.vec[ , 1] <- rep(num, k)
  num.vec[ , 2] <- rep(alpha, k)
  num.vec[ , 3] <- rep(sigma, k)

  # Generate independent Y's
  for(i in 1:k){
    cc <- rNPare(n = num, alpha = alpha, sigma = sigma)
    num.vec[i, 4] <- as.numeric(EM.case.2.alpha(alpha = MME.case2.alpha(sigma = sigma, data = cc), data = cc, sigma = sigma)$par)
    num.vec[i, 5] <- (num.vec[i, 4] - num.vec[i, 2])/num.vec[i, 2]
  }
  return(data.frame(num.vec))
}

l <- rbind(sim.MLE.case2.alpha(alpha = 0.7, sigma = 0.5, num = 50), sim.MLE.case2.alpha(alpha = 0.7, sigma = 0.5, num = 150), sim.MLE.case2.alpha(alpha = 0.7, sigma = 0.5, num = 300),
           sim.MLE.case2.alpha(alpha = 1, sigma = 1, num = 50), sim.MLE.case2.alpha(alpha = 1, sigma = 1, num = 150), sim.MLE.case2.alpha(alpha = 1, sigma = 1, num = 300),
           sim.MLE.case2.alpha(alpha = 2, sigma = 5, num = 50), sim.MLE.case2.alpha(alpha = 2, sigma = 5, num = 150), sim.MLE.case2.alpha(alpha = 2, sigma = 5, num = 300),
           sim.MLE.case2.alpha(alpha = 1, sigma = 0.5, num = 50), sim.MLE.case2.alpha(alpha = 2, sigma = 0.5, num = 50), sim.MLE.case2.alpha(alpha = 1, sigma = 0.5, num = 150),
           sim.MLE.case2.alpha(alpha = 2, sigma = 0.5, num = 150), sim.MLE.case2.alpha(alpha = 1, sigma = 0.5, num = 300), sim.MLE.case2.alpha(alpha = 2, sigma = 0.5, num = 300),
           sim.MLE.case2.alpha(alpha = 0.7, sigma = 1, num = 50), sim.MLE.case2.alpha(alpha = 2, sigma = 1, num = 50), sim.MLE.case2.alpha(alpha = 0.7, sigma = 1, num = 150),
           sim.MLE.case2.alpha(alpha = 2, sigma = 1, num = 150), sim.MLE.case2.alpha(alpha = 0.7, sigma = 1, num = 300), sim.MLE.case2.alpha(alpha = 2, sigma = 1, num = 300),
           sim.MLE.case2.alpha(alpha = 0.7, sigma = 5, num = 50), sim.MLE.case2.alpha(alpha = 1, sigma = 5, num = 50), sim.MLE.case2.alpha(alpha = 0.7, sigma = 5, num = 150),
           sim.MLE.case2.alpha(alpha = 1, sigma = 5, num = 150), sim.MLE.case2.alpha(alpha = 0.7, sigma = 5, num = 300), sim.MLE.case2.alpha(alpha = 1, sigma = 5, num = 300))

save(l, file = "MLE.case2.alpha.1000.simulations.RData")

library(ggplot2)
# Relative error graphs
gh <- ggplot(l, aes(x = RelativeError)) + geom_histogram(aes(y = ..density..)) + facet_grid(sigma ~ Sample + alpha, scales = "free", labeller = label_both) + theme_classic()
gh <- gh + labs(x = "Relative Error")
gh

# Boxplots
th <- ggplot(l, aes(x = as.factor(Sample), y = Estimate)) + geom_boxplot() + geom_hline(aes(yintercept = alpha), color = "red") +
  facet_wrap(sigma ~ alpha, scales = "free", labeller = label_both) + labs(x = "Sample size") + theme_classic()
th

# Generating table
tab1 <- l[which(l$Sample == 50 & l$sigma == 0.5 & l$alpha == 0.7), ]
(esti <- round(mean(tab1$Estimate), 4))
(MSE <- round(mean( (tab1$sigma - tab1$Estimate)^2 ), 4 ))

tab1 <- l[which(l$Sample == 300 & l$sigma == 5 & l$alpha == 2), ]
(esti <- round(mean(tab1$Estimate), 4 ))
(MSE <- round(mean( (tab1$sigma - tab1$Estimate)^2 ), 4 ))

# Case 3 alpha and sigma unknown
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

  while((diff > tol) & (k < maxiter)){
    #### M step
    alpha <- suppressWarnings(stats::nlm(function(x) length(data)*(log(x)/x + lgamma(1/x)) - sum(c_old - a_old)/x , p = alpha)$estimate)
    sigma <- sqrt(2)*mean(b_old*abs(data))
    expect <- all.expectations(iter.alpha = alpha, iter.sigma = sigma)
    a_old <- expect$a_k; b_old <- expect$b_k; c_old <- expect$c_k
    ll_new <- log_like(iter.alpha = alpha, iter.sigma = sigma)
    diff <- ll_new - ll_old
    ll_old <- ll_new
    k <- k + 1
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

library(rootSolve); library(cubature)
# MME
MME.case3.alpha.sigma <- function(data){
  R <- mean(log(abs(data))) - 0.5*(log(2) + digamma(1/2) + digamma(1)) - 2*log(gamma(3/4)*gamma(1/4)/(sqrt(pi)*(2^0.25)*mean(abs(data)^(-0.5))))

  fin <- function(alpha){
    h <- 0.5*(log(1/alpha) - digamma(1/alpha)) + 2*log((alpha^0.25)*gamma(1/alpha + 0.25)/gamma(1/alpha)) - R
    return(h)
  }
  alpha <- uniroot(fin, lower = 0.01, upper = 100, extendInt = "yes")$root
  sigma <- ( (alpha^0.25 * gamma(1/alpha + 0.25)*gamma(3/4)*gamma(1/4)  )/(mean(abs(data)^(-0.5))*gamma(1/2)*(2^0.25)*gamma(1/alpha)))^2
  return(c(alpha, sigma))
}

# Iteration plots
rt1 <- rNPare(n = 1000, alpha = 0.7, sigma = 0.5)
rt2 <- rNPare(n = 1000, alpha = 0.7, sigma = 1)
rt3 <- rNPare(n = 1000, alpha = 0.7, sigma = 5)
rt4 <- rNPare(n = 1000, alpha = 1, sigma = 0.5)
rt5 <- rNPare(n = 1000, alpha = 1, sigma = 1)
rt6 <- rNPare(n = 1000, alpha = 1, sigma = 5)
rt7 <- rNPare(n = 1000, alpha = 2, sigma = 0.5)
rt8 <- rNPare(n = 1000, alpha = 2, sigma = 1)
rt9 <- rNPare(n = 1000, alpha = 2, sigma = 5)

l1 <- EM.case.3.alpha.sigma(data = rt1, alpha =  1.5, sigma = 1)
l2 <- EM.case.3.alpha.sigma(data = rt1, alpha = MME.case3.alpha.sigma(data = rt1)[1], sigma = MME.case3.alpha.sigma(data = rt1)[2])
l3 <- EM.case.3.alpha.sigma(data = rt1, alpha = 0.3, sigma = 0.2)

plot(l1$output$iteration, l1$output$alpha, ylim = c(0, max(l1$output$alpha,l2$output$alpha,l3$output$alpha)), type = "o", pch = 19,
     xlab = "Iteration", ylab = expression(paste(alpha, " values")), xlim = c(0, max(l1$output$iteration, l2$output$iteration, l3$output$iteration)))
points(l2$output$iteration, l2$output$alpha, type = "o", pch = 19, col = "blue")
points(l3$output$iteration, l3$output$alpha, type = "o", col = "green", pch = 19)
abline(h = 0.7, col = "red", lty = 2)
legend(15, 1.45, legend = c("Intial values", expression(paste(alpha, " = 0.3014, ", sigma, " = 0.6055 (MME)")),
                           expression(paste(alpha, " = 1.5, ", sigma, " = 1" )), expression(paste(alpha, " = 0.3, ", sigma, " = 0.2"))),
       col = c("", "blue", "black", "green"), lty = c(NA, 1, 1,1), pch = c(NA, 19, 19, 19), bty = "n",
       cex = 1.2, seg.len = 1, x.intersp = 0.1, y.intersp = 0.8)

plot(l1$output$iteration, l1$output$sigma, ylim = c(0, max(l1$output$sigma,l2$output$sigma,l3$output$sigma)), type = "o", pch = 19,
     xlab = "Iteration", ylab = expression(paste(sigma, " values")), xlim = c(0, max(l1$output$iteration, l2$output$iteration, l3$output$iteration)))
points(l2$output$iteration, l2$output$sigma, type = "o", pch = 19, col = "blue")
points(l3$output$iteration, l3$output$sigma, type = "o", col = "green", pch = 19)
abline(h = 0.5, col = "red", lty = 2)
legend(15, 0.4, legend = c("Intial values", expression(paste(alpha, " = 0.3014; ", sigma, " = 0.6055 (MME)")),
                           expression(paste(alpha, " = 1.5, ", sigma, " = 1" )), expression(paste(alpha, " = 0.3, ", sigma, " = 0.2"))),
       col = c("", "blue", "black", "green"), lty = c(NA, 1, 1,1), pch = c(NA, 19, 19, 19), bty = "n",
       cex = 1.2, seg.len = 1, x.intersp = 0.1, y.intersp = 0.8)

m1 <- EM.case.3.alpha.sigma(data = rt2, alpha = 1.5, sigma = 3)
m2 <- EM.case.3.alpha.sigma(data = rt2, alpha = MME.case3.alpha.sigma(data = rt2)[1], sigma = MME.case3.alpha.sigma(data = rt2)[2])
m3 <- EM.case.3.alpha.sigma(data = rt2, alpha = 0.3, sigma = 0.5)

plot(m1$output$iteration, m1$output$alpha, ylim = c(0, max(m1$output$alpha,m2$output$alpha,m3$output$alpha)), type = "o", pch = 19,
     xlab = "Iteration", ylab = expression(paste(alpha, " values")), xlim = c(0, max(m1$output$iteration, m2$output$iteration, m3$output$iteration)))
points(m2$output$iteration, m2$output$alpha, type = "o", pch = 19, col = "blue")
points(m3$output$iteration, m3$output$alpha, type = "o", col = "green", pch = 19)
abline(h = 0.7, col = "red", lty = 2)
legend(8, 0.5, legend = c("Intial values", expression(paste(alpha, " = 0.8969 (MME)" )),
                          expression(paste(alpha, " = 1.5" )), expression(paste(alpha, " = 0.5"))),
       col = c("", "blue", "green", "black"), lty = c(NA, 1, 1,1), pch = c(NA, 19, 19, 19), bty = "n",
       cex = 1.2, seg.len = 1, x.intersp = 0.1, y.intersp = 0.8)

n1 <- EM.case.2.alpha(data = rt3, alpha = 0.7, sigma = 1)
n2 <- EM.case.2.alpha(data = rt3, alpha = MME.case2.alpha(sigma = 1, data = rt3), sigma = 1)
n3 <- EM.case.2.alpha(data = rt3, alpha = 5, sigma = 1)

plot(n1$output$iteration, n1$output$alpha, ylim = c(0, max(n1$output$alpha,n2$output$alpha,n3$output$alpha)), type = "o", pch = 19,
     xlab = "Iteration", ylab = expression(paste(alpha, " values")), xlim = c(0, max(n1$output$iteration, n2$output$iteration, n3$output$iteration)))
points(n2$output$iteration, n2$output$alpha, type = "o", pch = 19, col = "blue")
points(n3$output$iteration, n3$output$alpha, type = "o", col = "green", pch = 19)
abline(h = 2, col = "red", lty = 2)
legend(5, 4, legend = c("Intial values", expression(paste(alpha, " = 1.7960 (MME)")),
                        expression(paste(alpha, " = 5" )), expression(paste(alpha, " = 0.7"))),
       col = c("", "blue", "green", "black"), lty = c(NA, 1, 1,1), pch = c(NA, 19, 19, 19), bty = "n",
       cex = 1.2, seg.len = 1, x.intersp = 0.1, y.intersp = 0.8)

o1 <- EM.case.2.alpha(data = rt3, alpha = 0.7, sigma = 1)
o2 <- EM.case.2.alpha(data = rt3, alpha = MME.case2.alpha(sigma = 1, data = rt3), sigma = 1)
o3 <- EM.case.2.alpha(data = rt3, alpha = 5, sigma = 1)

plot(n1$output$iteration, n1$output$alpha, ylim = c(0, max(n1$output$alpha,n2$output$alpha,n3$output$alpha)), type = "o", pch = 19,
     xlab = "Iteration", ylab = expression(paste(alpha, " values")), xlim = c(0, max(n1$output$iteration, n2$output$iteration, n3$output$iteration)))
points(n2$output$iteration, n2$output$alpha, type = "o", pch = 19, col = "blue")
points(n3$output$iteration, n3$output$alpha, type = "o", col = "green", pch = 19)
abline(h = 2, col = "red", lty = 2)
legend(5, 4, legend = c("Intial values", expression(paste(alpha, " = 1.7960 (MME)")),
                        expression(paste(alpha, " = 5" )), expression(paste(alpha, " = 0.7"))),
       col = c("", "blue", "green", "black"), lty = c(NA, 1, 1,1), pch = c(NA, 19, 19, 19), bty = "n",
       cex = 1.2, seg.len = 1, x.intersp = 0.1, y.intersp = 0.8)

p1 <- EM.case.2.alpha(data = rt3, alpha = 0.7, sigma = 1)
p2 <- EM.case.2.alpha(data = rt3, alpha = MME.case2.alpha(sigma = 1, data = rt3), sigma = 1)
p3 <- EM.case.2.alpha(data = rt3, alpha = 5, sigma = 1)

plot(n1$output$iteration, n1$output$alpha, ylim = c(0, max(n1$output$alpha,n2$output$alpha,n3$output$alpha)), type = "o", pch = 19,
     xlab = "Iteration", ylab = expression(paste(alpha, " values")), xlim = c(0, max(n1$output$iteration, n2$output$iteration, n3$output$iteration)))
points(n2$output$iteration, n2$output$alpha, type = "o", pch = 19, col = "blue")
points(n3$output$iteration, n3$output$alpha, type = "o", col = "green", pch = 19)
abline(h = 2, col = "red", lty = 2)
legend(5, 4, legend = c("Intial values", expression(paste(alpha, " = 1.7960 (MME)")),
                        expression(paste(alpha, " = 5" )), expression(paste(alpha, " = 0.7"))),
       col = c("", "blue", "green", "black"), lty = c(NA, 1, 1,1), pch = c(NA, 19, 19, 19), bty = "n",
       cex = 1.2, seg.len = 1, x.intersp = 0.1, y.intersp = 0.8)

q1 <- EM.case.2.alpha(data = rt3, alpha = 0.7, sigma = 1)
q2 <- EM.case.2.alpha(data = rt3, alpha = MME.case2.alpha(sigma = 1, data = rt3), sigma = 1)
q3 <- EM.case.2.alpha(data = rt3, alpha = 5, sigma = 1)

plot(n1$output$iteration, n1$output$alpha, ylim = c(0, max(n1$output$alpha,n2$output$alpha,n3$output$alpha)), type = "o", pch = 19,
     xlab = "Iteration", ylab = expression(paste(alpha, " values")), xlim = c(0, max(n1$output$iteration, n2$output$iteration, n3$output$iteration)))
points(n2$output$iteration, n2$output$alpha, type = "o", pch = 19, col = "blue")
points(n3$output$iteration, n3$output$alpha, type = "o", col = "green", pch = 19)
abline(h = 2, col = "red", lty = 2)
legend(5, 4, legend = c("Intial values", expression(paste(alpha, " = 1.7960 (MME)")),
                        expression(paste(alpha, " = 5" )), expression(paste(alpha, " = 0.7"))),
       col = c("", "blue", "green", "black"), lty = c(NA, 1, 1,1), pch = c(NA, 19, 19, 19), bty = "n",
       cex = 1.2, seg.len = 1, x.intersp = 0.1, y.intersp = 0.8)

r1 <- EM.case.2.alpha(data = rt3, alpha = 0.7, sigma = 1)
r2 <- EM.case.2.alpha(data = rt3, alpha = MME.case2.alpha(sigma = 1, data = rt3), sigma = 1)
r3 <- EM.case.2.alpha(data = rt3, alpha = 5, sigma = 1)

plot(n1$output$iteration, n1$output$alpha, ylim = c(0, max(n1$output$alpha,n2$output$alpha,n3$output$alpha)), type = "o", pch = 19,
     xlab = "Iteration", ylab = expression(paste(alpha, " values")), xlim = c(0, max(n1$output$iteration, n2$output$iteration, n3$output$iteration)))
points(n2$output$iteration, n2$output$alpha, type = "o", pch = 19, col = "blue")
points(n3$output$iteration, n3$output$alpha, type = "o", col = "green", pch = 19)
abline(h = 2, col = "red", lty = 2)
legend(5, 4, legend = c("Intial values", expression(paste(alpha, " = 1.7960 (MME)")),
                        expression(paste(alpha, " = 5" )), expression(paste(alpha, " = 0.7"))),
       col = c("", "blue", "green", "black"), lty = c(NA, 1, 1,1), pch = c(NA, 19, 19, 19), bty = "n",
       cex = 1.2, seg.len = 1, x.intersp = 0.1, y.intersp = 0.8)

s1 <- EM.case.2.alpha(data = rt3, alpha = 0.7, sigma = 1)
s2 <- EM.case.2.alpha(data = rt3, alpha = MME.case2.alpha(sigma = 1, data = rt3), sigma = 1)
s3 <- EM.case.2.alpha(data = rt3, alpha = 5, sigma = 1)

plot(n1$output$iteration, n1$output$alpha, ylim = c(0, max(n1$output$alpha,n2$output$alpha,n3$output$alpha)), type = "o", pch = 19,
     xlab = "Iteration", ylab = expression(paste(alpha, " values")), xlim = c(0, max(n1$output$iteration, n2$output$iteration, n3$output$iteration)))
points(n2$output$iteration, n2$output$alpha, type = "o", pch = 19, col = "blue")
points(n3$output$iteration, n3$output$alpha, type = "o", col = "green", pch = 19)
abline(h = 2, col = "red", lty = 2)
legend(5, 4, legend = c("Intial values", expression(paste(alpha, " = 1.7960 (MME)")),
                        expression(paste(alpha, " = 5" )), expression(paste(alpha, " = 0.7"))),
       col = c("", "blue", "green", "black"), lty = c(NA, 1, 1,1), pch = c(NA, 19, 19, 19), bty = "n",
       cex = 1.2, seg.len = 1, x.intersp = 0.1, y.intersp = 0.8)

t1 <- EM.case.2.alpha(data = rt3, alpha = 0.7, sigma = 1)
t2 <- EM.case.2.alpha(data = rt3, alpha = MME.case2.alpha(sigma = 1, data = rt3), sigma = 1)
t3 <- EM.case.2.alpha(data = rt3, alpha = 5, sigma = 1)

plot(n1$output$iteration, n1$output$alpha, ylim = c(0, max(n1$output$alpha,n2$output$alpha,n3$output$alpha)), type = "o", pch = 19,
     xlab = "Iteration", ylab = expression(paste(alpha, " values")), xlim = c(0, max(n1$output$iteration, n2$output$iteration, n3$output$iteration)))
points(n2$output$iteration, n2$output$alpha, type = "o", pch = 19, col = "blue")
points(n3$output$iteration, n3$output$alpha, type = "o", col = "green", pch = 19)
abline(h = 2, col = "red", lty = 2)
legend(5, 4, legend = c("Intial values", expression(paste(alpha, " = 1.7960 (MME)")),
                        expression(paste(alpha, " = 5" )), expression(paste(alpha, " = 0.7"))),
       col = c("", "blue", "green", "black"), lty = c(NA, 1, 1,1), pch = c(NA, 19, 19, 19), bty = "n",
       cex = 1.2, seg.len = 1, x.intersp = 0.1, y.intersp = 0.8)


first <- list(rt1,l1,l2,l3)
second <- list(rt2,m1,m2,m3)
third <-list(rt3,n1,n2,n3)

save(first, file = "EM.data.1000.alpha.unknown.0.7.sigma.known.1.RData")
save(second, file = "EM.data.1000.alpha.unknown.1.sigma.known.1.RData")
save(third, file = "EM.data.1000.alpha.unknown.2.sigma.known.1.RData")

sim.MLE.case2.alpha <- function(k = 1000, alpha = 0.5, sigma = 1, num = 10){
  num.vec <- matrix(NA, nrow = k, ncol = 5)
  colnames(num.vec) <- c("Sample", "alpha", "sigma", "Estimate", "RelativeError")
  num.vec[ , 1] <- rep(num, k)
  num.vec[ , 2] <- rep(alpha, k)
  num.vec[ , 3] <- rep(sigma, k)

  # Generate independent Y's
  for(i in 1:k){
    cc <- rNPare(n = num, alpha = alpha, sigma = sigma)
    num.vec[i, 4] <- as.numeric(EM.case.2.alpha(alpha = MME.case2.alpha(sigma = sigma, data = cc), data = cc, sigma = sigma)$par)
    num.vec[i, 5] <- (num.vec[i, 4] - num.vec[i, 2])/num.vec[i, 2]
  }
  return(data.frame(num.vec))
}

l <- rbind(sim.MLE.case2.alpha(alpha = 0.7, sigma = 0.5, num = 50), sim.MLE.case2.alpha(alpha = 0.7, sigma = 0.5, num = 150), sim.MLE.case2.alpha(alpha = 0.7, sigma = 0.5, num = 300),
           sim.MLE.case2.alpha(alpha = 1, sigma = 1, num = 50), sim.MLE.case2.alpha(alpha = 1, sigma = 1, num = 150), sim.MLE.case2.alpha(alpha = 1, sigma = 1, num = 300),
           sim.MLE.case2.alpha(alpha = 2, sigma = 5, num = 50), sim.MLE.case2.alpha(alpha = 2, sigma = 5, num = 150), sim.MLE.case2.alpha(alpha = 2, sigma = 5, num = 300),
           sim.MLE.case2.alpha(alpha = 1, sigma = 0.5, num = 50), sim.MLE.case2.alpha(alpha = 2, sigma = 0.5, num = 50), sim.MLE.case2.alpha(alpha = 1, sigma = 0.5, num = 150),
           sim.MLE.case2.alpha(alpha = 2, sigma = 0.5, num = 150), sim.MLE.case2.alpha(alpha = 1, sigma = 0.5, num = 300), sim.MLE.case2.alpha(alpha = 2, sigma = 0.5, num = 300),
           sim.MLE.case2.alpha(alpha = 0.7, sigma = 1, num = 50), sim.MLE.case2.alpha(alpha = 2, sigma = 1, num = 50), sim.MLE.case2.alpha(alpha = 0.7, sigma = 1, num = 150),
           sim.MLE.case2.alpha(alpha = 2, sigma = 1, num = 150), sim.MLE.case2.alpha(alpha = 0.7, sigma = 1, num = 300), sim.MLE.case2.alpha(alpha = 2, sigma = 1, num = 300),
           sim.MLE.case2.alpha(alpha = 0.7, sigma = 5, num = 50), sim.MLE.case2.alpha(alpha = 1, sigma = 5, num = 50), sim.MLE.case2.alpha(alpha = 0.7, sigma = 5, num = 150),
           sim.MLE.case2.alpha(alpha = 1, sigma = 5, num = 150), sim.MLE.case2.alpha(alpha = 0.7, sigma = 5, num = 300), sim.MLE.case2.alpha(alpha = 1, sigma = 5, num = 300))

save(l, file = "MLE.case2.alpha.1000.simulations.RData")

library(ggplot2)
# Relative error graphs
gh <- ggplot(l, aes(x = RelativeError)) + geom_histogram(aes(y = ..density..)) + facet_grid(sigma ~ Sample + alpha, scales = "free", labeller = label_both) + theme_classic()
gh <- gh + labs(x = "Relative Error")
gh

# Boxplots
th <- ggplot(l, aes(x = as.factor(Sample), y = Estimate)) + geom_boxplot() + geom_hline(aes(yintercept = alpha), color = "red") +
  facet_wrap(sigma ~ alpha, scales = "free", labeller = label_both) + labs(x = "Sample size") + theme_classic()
th

# Generating table
tab1 <- l[which(l$Sample == 50 & l$sigma == 0.5 & l$alpha == 0.7), ]
(esti <- round(mean(tab1$Estimate), 4))
(MSE <- round(mean( (tab1$sigma - tab1$Estimate)^2 ), 4 ))

tab1 <- l[which(l$Sample == 300 & l$sigma == 5 & l$alpha == 2), ]
(esti <- round(mean(tab1$Estimate), 4 ))
(MSE <- round(mean( (tab1$sigma - tab1$Estimate)^2 ), 4 ))

# This is for Case 2
load(file = "MLE.case2.alpha.1000.simulations.RData")
tab1 <- l[which(l$Sample == 300 & l$sigma == 5 & l$alpha == 2), ]
(esti <- round(mean(tab1$Estimate), 4 ))
(MSE <- round(mean( (tab1$alpha - tab1$Estimate)^2 ), 4 ))

# Relative error graphs
gh <- ggplot(l, aes(x = RelativeError)) + geom_histogram(aes(y = ..density..)) + facet_grid( alpha ~ sigma + Sample, scales = "free", labeller = label_both) + theme_classic()
gh <- gh + labs(x = "Relative Error")
gh


# Boxplots
th <- ggplot(l, aes(x = as.factor(Sample), y = Estimate)) + geom_boxplot() + geom_hline(aes(yintercept = alpha), color = "red") +
  facet_wrap(alpha ~ sigma, scales = "free", labeller = label_both) + labs(x = "Sample size") + theme_classic()
th



# This is for case 3
load(file = "MLE.case3.alpha.sigma.1.simulations.RData")
load(file = "MLE.case3.alpha.sigma.4.simulations.RData")
load(file = "MLE.case3.alpha.sigma.5.simulations.RData")
load(file = "MLE.case3.alpha.sigma.7.simulations.RData")
load(file = "MLE.case3.alpha.sigma.8.simulations.RData")
load(file = "MLE.case3.alpha.sigma.10.simulations.RData")
load(file = "MLE.case3.alpha.sigma.11.simulations.RData")
load(file = "MLE.case3.alpha.sigma.13.simulations.RData")
load(file = "MLE.case3.alpha.sigma.14.simulations.RData")
load(file = "MLE.case3.alpha.sigma.16.simulations.RData")
load(file = "MLE.case3.alpha.sigma.17.simulations.RData")
load(file = "MLE.case3.alpha.sigma.19.simulations.RData")
load(file = "MLE.case3.alpha.sigma.20.simulations.RData")
load(file = "MLE.case3.alpha.sigma.22.simulations.RData")
load(file = "MLE.case3.alpha.sigma.23.simulations.RData")
load(file = "MLE.case3.alpha.sigma.25.simulations.RData")
load(file = "MLE.case3.alpha.sigma.26.simulations.RData")

l <- rbind(l1[complete.cases(l1), ], l4[complete.cases(l4), ], l5[complete.cases(l5), ], l7[complete.cases(l7), ], l8[complete.cases(l8), ], l10[complete.cases(l10), ], l11[complete.cases(l11), ],
           l13[complete.cases(l13), ], l14[complete.cases(l14), ], l16[complete.cases(l16), ], l17[complete.cases(l17), ], l19[complete.cases(l19), ], l20[complete.cases(l20), ], l22[complete.cases(l22), ],
           l23[complete.cases(l23), ], l25[complete.cases(l25), ], l26[complete.cases(l26), ])

# Plots for alpha
# Relative error graphs
gh <- ggplot(l, aes(x = RelativeErrorAlpha)) + geom_histogram(aes(y = ..density..)) + facet_grid( alpha ~ sigma + Sample, scales = "free", labeller = label_both) + theme_classic()
gh <- gh + labs(x = "Relative Error")
gh


gh <- ggplot(l, aes(x = RelativeErrorAlpha)) + geom_histogram(aes(y = ..density..)) + facet_grid( alpha ~ sigma + Sample, scales = "free", labeller = label_both) + theme_classic()
gh <- gh + labs(x = "Relative Error")
gh

# Boxplots
th <- ggplot(l, aes(x = as.factor(Sample), y = EstimateAlpha)) + geom_boxplot() + geom_hline(aes(yintercept = alpha), color = "red") +
  facet_wrap(alpha ~ sigma, scales = "free", labeller = label_both) + labs(x = "Sample size", y = expression(paste("Estimate ", alpha))) + theme_classic()
th

# Plots for sigma
gh <- ggplot(l, aes(x = RelativeErrorSigma)) + geom_histogram(aes(y = ..density..)) + facet_grid(alpha ~ sigma + Sample, scales = "free", labeller = label_both) + theme_classic()
gh <- gh + labs(x = "Relative Error")
gh

# Boxplots
th <- ggplot(l, aes(x = as.factor(Sample), y = EstimateSigma)) + geom_boxplot() + geom_hline(aes(yintercept = sigma), color = "red") +
  facet_wrap(alpha ~ sigma, scales = "free", labeller = label_both) + labs(x = "Sample size", y = expression(paste("Estimate ", sigma))) + theme_classic()
th
