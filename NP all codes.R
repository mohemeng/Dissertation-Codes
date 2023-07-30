library(cubature);library(rootSolve)
####---- The PDF -------#####
NP.pdf <- function(y, alpha = 0.5, sigma = 1){
  bb <- array(NA, dim = length(y))
  for(i in y){
    pos <- match(i, y)
    gg <- function(t){
      uu <- (sqrt(t/2)/sigma)*exp(-sqrt(2*t)*abs(i)/sigma)*(t^(1/alpha - 1)*exp(-t/alpha))/(alpha^(1/alpha)*gamma(1/alpha))
      return(uu)
    }
  bb[pos] <- adaptIntegrate(gg, lowerLimit = 0, upperLimit = Inf)$integral
  #bb[pos] <- integrate(gg, lower = 0, upper = Inf)$value
  }
return(bb)
}

###----- Curve for PDF -----#####
graph.f.x.NP <- function(alpha.par = 0.5, sigma.par = 1, start = -10, end = 10){
  NP.pdf <- function(y, alpha= alpha.par, sigma = sigma.par){
    bb <- array(NA, dim = length(y))
    for(i in y){
      pos <- match(i, y)
      gg <- function(t){
        uu <- (sqrt(t/2)/sigma)*exp(-sqrt(2*t)*abs(i)/sigma)*(t^(1/alpha - 1)*exp(-t/alpha))/(alpha^(1/alpha)*gamma(1/alpha))
        return(uu)
      }
      bb[pos] <- integrate(gg, lower = 0, upper = Inf)$value
    }
    return(bb)
  }
  curve(NP.pdf, from = start, to = end, ylab = "f(y)", xlab = "y", col = "blue")
  legend(x = 3, y = 0.3, legend = paste(c("alpha =", alpha.par, ";", "sigma =", sigma.par),
                                        collapse = " "), lty = 1, col =c("blue"), cex = 0.7, box.col = "white", text.width = 3)
}

graph.f.x.NP(alpha.par = 0.5, sigma.par = 1)

### PDF plots
val <- seq(-10, 10, by = 0.01)
f.1 <- NP.pdf(y = val, alpha = 0.6, sigma = 3)
f.2 <- NP.pdf(y = val, alpha = 1, sigma = 1)
f.3 <- NP.pdf(y = val, alpha = 1, sigma = 1.5)

plot(val, f.1, type = "l", col = "black", xlab = "y", ylab = "f(y)", ylim = c(0, max(f.1, f.2, f.3)))
points(val, f.2, type = "l", lty = "dotted", col = "red")
points(val, f.3, type = "l", lty = "twodash", col = "blue")
legend(x = "right", legend = c("alpha = 0.6, sigma = 3", "alpha = 1, sigma = 1",
                                    "alpha = 2, sigma = 1.5"), lty = c(1, 3, 5), col = c("black", "red", "blue"), cex = 0.7, text.width = 2.5)

###---- The CDF ------#####
NP.cdf <- function(y, alpha = 1, sigma = 1){
 bb <- array(NA, dim = length(y))
  for(i in y){
   pos <- match(i, y)
    gg <- function(t){(0.5*exp(-sqrt(2*t)*abs(i)/sigma)*t^((1/alpha) - 1)*exp(-t/alpha))/(alpha^(1/alpha)*gamma(1/alpha))}
  if(i < 0){
  bb[pos] <- integrate(gg, lower = 0, upper = Inf)$value
    }
  else if(i >= 0){
  bb[pos] <- 1 - integrate(gg, lower = 0, upper = Inf)$value
  }
 }
return(bb)
}

####----- Curve for CDF -----####
curve(NP.cdf, from = -10, to = 10, ylab = "F(y)", xlab = "y")

### Other CDFs with different parameter combinations
val <- seq(-13, 13, by = 0.01)
f.1 <- NP.cdf(y = val, alpha = 0.6, sigma = 3)
f.2 <- NP.cdf(y = val, alpha = 1, sigma = 1)
f.3 <- NP.cdf(y = val, alpha = 2, sigma = 1.5)

plot(val, f.1, type = "l", col = "black", xlab = "y", ylab = "F(y)")
points(val, f.2, type = "l", lty = "dotted", col = "red")
points(val, f.3, type = "l", lty = "twodash", col = "blue")
legend(x = 0.5, y = 0.4, legend = c("alpha = 0.6, sigma = 3", "alpha = 1, sigma = 1",
                                    "alpha = 2, sigma = 1.5"), lty = c(1, 3, 5), col = c("black", "red", "blue"), cex = 0.7, text.width = 8)


# Obtain n random samples from NP
rNPare <- function(n = 1, alpha = 0.5, sigma = 1){
    normal <- rnorm(n); exponential <- rexp(n)
    gamm <- rgamma(n, shape = 1/alpha, scale = alpha)
    y <- sigma*sqrt(exponential/gamm)*normal
    return(y)
}

# Graph overlay of numbers generated versus the density
par(mfrow = c(2, 2))

alpha.par <- 0.7; sigma.par <- 1
X <- rNPare(n = 1000, alpha = alpha.par, sigma = sigma.par)
jj <- NP.pdf(y = seq(min(X), max(X), by = 0.01), alpha = alpha.par, sigma = sigma.par)
hist(X, freq = FALSE, breaks = 60, ylim = c(0, max(jj)), main = "", xlab = "y")
points(x = seq(min(X), max(X), by = 0.01), y = jj, type = "l", col = "red")
legend(x = 3, y = 0.4, legend = expression(paste(alpha, " = ", 0.7, "; ", sigma, " = ", 1)), bty = "n", cex = 1.2, text.width = 6.5 )

alpha.par <- 1; sigma.par <- 1
X <- rNPare(n = 1000, alpha = alpha.par, sigma = sigma.par)
jj <- NP.pdf(y = seq(min(X), max(X), by = 0.01), alpha = alpha.par, sigma = sigma.par)
hist(X, freq = FALSE, breaks = 60, ylim = c(0, max(jj)), main = "", xlab = "y")
points(x = seq(min(X), max(X), by = 0.01), y = jj, type = "l", col = "red")
legend(x = 3, y = 0.4, legend = expression(paste(alpha, " = ", 1, "; ", sigma, " = ", 1)), bty = "n", cex = 1.2, text.width = 6.5 )

#alpha.par <- 2; sigma.par <- 1
#X <- rNPare(n = 300, alpha = alpha.par, sigma = sigma.par)
#jj <- NP.pdf(y = seq(min(X), max(X), by = 0.01), alpha = alpha.par, sigma = sigma.par)
#hist(X, freq = FALSE, breaks = 120, ylim = c(0, max(jj)), main = "", xlab = "y")
#points(x = seq(min(X), max(X), by = 0.01), y = jj, type = "l", col = "red")
#legend(x = -400, y = 0.4, legend = expression(paste(alpha, " = ", 2, "; ", sigma, " = ", 1)), bty = "n", cex = 1.2, text.width = 6.5 )

alpha.par <- 0.7; sigma.par <- 5
X <- rNPare(n = 1000, alpha = alpha.par, sigma = sigma.par)
jj <- NP.pdf(y = seq(min(X), max(X), by = 0.01), alpha = alpha.par, sigma = sigma.par)
hist(X, freq = FALSE, breaks = 60, ylim = c(0, max(jj)), main = "", xlab = "y")
points(x = seq(min(X), max(X), by = 0.01), y = jj, type = "l", col = "red")
legend(x = -100, y = 0.06, legend = expression(paste(alpha, " = ", 0.7, "; ", sigma, " = ", 5)), bty = "n", cex = 1.2, text.width = 6.5 )

alpha.par <- 1; sigma.par <- 5
X <- rNPare(n = 1000, alpha = alpha.par, sigma = sigma.par)
jj <- NP.pdf(y = seq(min(X), max(X), by = 0.01), alpha = alpha.par, sigma = sigma.par)
hist(X, freq = FALSE, breaks = 60, ylim = c(0, max(jj)), main = "", xlab = "y")
points(x = seq(min(X), max(X), by = 0.01), y = jj, type = "l", col = "red")
legend(x = -100, y = 0.06, legend = expression(paste(alpha, " = ", 1, "; ", sigma, " = ", 5)), bty = "n", cex = 1.2, text.width = 6.5 )

#alpha.par <- 2; sigma.par <- 5
#X <- rNPare(n = 300, alpha = alpha.par, sigma = sigma.par)
#jj <- NP.pdf(y = seq(min(X), max(X), by = 0.01), alpha = alpha.par, sigma = sigma.par)
#hist(X, freq = FALSE, breaks = 60, ylim = c(0, max(jj)), main = "", xlab = "y")
#points(x = seq(min(X), max(X), by = 0.01), y = jj, type = "l", col = "red")
#legend(x = 3, y = 0.4, legend = expression(paste(alpha, " = ", 2, "; ", sigma, " = ", 5)), bty = "n", cex = 1.2, text.width = 6.5 )





#alpha.par <- 1; sigma.par <- 1
#X <- rNPare(n = 1000, alpha = alpha.par, sigma = sigma.par)
#jj <- NP.pdf(y = seq(min(X), max(X), by = 0.01), alpha = alpha.par, sigma = sigma.par)
#hist(X, freq = FALSE, ylim = c(0, max(jj)+0.12), main = "", col = "gray", xlab = "y")
#points(x = seq(min(X), max(X), by = 0.01), y = jj, type = "l", col = "red")
#legend(x = "top", legend = paste(c("alpha =", alpha.par, ";", "sigma =", sigma.par), collapse = " "),box.col = "white", cex = 0.8, text.width = 6.5)

#alpha.par <- 1.5; sigma.par <- 1.5
#X <- rNPare(n = 1000, alpha = alpha.par, sigma = sigma.par)
#jj <- NP.pdf(y = seq(min(X), max(X), by = 0.01), alpha = alpha.par, sigma = sigma.par)
#hist(X, freq = FALSE, ylim = c(0, max(jj)+0.12), main = "", col = "gray", xlab = "y")
#points(x = seq(min(X), max(X), by = 0.01), y = jj, type = "l", col = "red")
#legend(x = "top", legend = paste(c("alpha =", alpha.par, ";", "sigma =", sigma.par), collapse = " "),box.col = "white", cex = 0.8, text.width = 6.5)

#alpha.par <- 0.1; sigma.par <- 10
#X <- rNPare(n = 1000, alpha = alpha.par, sigma = sigma.par)
#jj <- NP.pdf(y = seq(min(X), max(X), by = 0.01), alpha = alpha.par, sigma = sigma.par)
#hist(X, freq = FALSE, ylim = c(0, max(jj)+0.02), main = "", col = "gray", xlab = "y")
#points(x = seq(min(X), max(X), by = 0.01), y = jj, type = "l", col = "red")
#legend(x = "top", legend = paste(c("alpha =", alpha.par, ";", "sigma =", sigma.par), collapse = " "),box.col = "white", cex = 0.8, text.width = 6.5)

par(mfrow = c(1, 1))


#####------------ Quantile Function -----------------######
#NP.quantile <- function(alpha = 0.5, sigma = 1, u = 0.1){
#  sam <- rgamma(10000, shape = 1/alpha, scale = alpha)

#  g.y.1 <- function(y){
#    fin <- mean( 0.5*exp(-sqrt(2*sam)*abs(y))) - u
#    return(fin)}
#  g.y.2 <- function(y){
#    fin <- 1 - u - mean(0.5*exp(-sqrt(2*sam)*abs(y)))
#    return(fin)}
#  require(rootSolve)

#  ahy <- array(NA, dim = length(u))

#  for(i in u){
#    pos <- match(i, u)

#    if(i == 0.5){
#      ahy[pos] <- 0
#    }
#    else if(i < 0.5){
#      ahy[pos] <- return(uniroot(g.y.1, lower = -15, upper = -0.01,
#                                 extendInt = "yes", maxiter = 100 )$root)
#    }
#    else if(i > 0.5){
#      ahy[pos] <- return(uniroot(g.y.2, lower = 0.001, upper = 15,
#                                 extendInt = "yes", maxiter = 100 )$root)
#    }
#  }
#  return(ahy)
#}


#NP.quantile <- function(alpha = 0.5, sigma = 1, u = 0.1){
#  sam <- rgamma(10000, shape = 1/alpha, scale = alpha)

#  g.y.1 <- function(y){
#    fin <- mean( 0.5*exp(-sqrt(2*sam)*abs(y))) - u
#    return(fin)}
#  g.y.2 <- function(y){
#    fin <- 1 - u - mean(0.5*exp(-sqrt(2*sam)*abs(y)))
#    return(fin)}
#  require(rootSolve)

#  ahy <- array(NA, dim = length(u))

#  for(i in u){
#    pos <- match(i, u)

#    if(i == 0.5){
#      ahy[pos] <- 0
#    }
#    if(i < 0.5){
#      ahy[pos] <- return(uniroot(g.y.1, lower = -15, upper = -0.01,
#                                 extendInt = "yes", maxiter = 100 )$root)
#    }
#    if(i > 0.5){
#      ahy[pos] <- return(uniroot(g.y.2, lower = 0.001, upper = 15,
#                                 extendInt = "yes", maxiter = 100 )$root)
#    }
#  }
#  return(ahy)
#}


##################################################
#Y.quantile <- function(alpha, sigma, u = 0.1){
 # sam <- rgamma(10000, shape = 1/alpha, scale = alpha)

#  require(rootSolve)

 # ahy <- array(NA, dim = length(u))

#  for(i in u){
  #  pos <- match(i, u)

  #  g.y.1 <- function(y){
  #    fin <- mean( 0.5*exp(-sqrt(2*sam)*abs(y))) - i
  #    return(fin)}
 #   g.y.2 <- function(y){
#      fin <- 1 - i - mean(0.5*exp(-sqrt(2*sam)*abs(y)))
 #     return(fin)}

#    if(i == 0.5){
#      ahy[pos] <- 0
 #   }
#    else if(i < 0.5){
 #     ahy[pos] <- return(uniroot(g.y.1, lower = -15, upper = -0.00000001,
#                                 extendInt = "yes", maxiter = 10000 )$root)
 #   }
#    else if(i > 0.5){
 #     ahy[pos] <- return(uniroot(g.y.2, lower = 0.00000001, upper = 15,
#                                 extendInt = "yes", maxiter = 10000 )$root)
#    }
#  }
#  return(ahy)
#}


NP.quantile <- function(alpha = 0.5, sigma = 1, u = 0.01){
  ahy <- array(NA, dim = length(u))
  for(i in u){
  pos <- match(i, u)

  g.y.1 <- function(y){
    fin <- NP.cdf(y = y, alpha = alpha, sigma = sigma) - i
    return(fin)
  }
  ahy[pos] <- uniroot(g.y.1, lower = -0.001, upper = 100, extendInt = "yes")$root
  }
  return(ahy)
}


# A graph of the quantile function
plot(x = seq(0.01, 0.99, 0.01), y = NP.quantile(u = seq(0.01, 0.99, 0.01)), type = "l", xlab = "u", ylab = "Quantile")

#emp.cdf <- function(data, cor = 0.25){
#  s <- sort(data); n <- length(data)
#  uni <- unique(sort(data))
#  cdf <- array(NA, dim = length(uni))
#  for(i in uni){
#    pos <- match(i, uni)
#    cdf[pos] <- length(which(s <= i))/n
#  }
#  mon <- n/(n + cor)  # using correction for final (tail) value
#  cdf[length(uni)] <- mon
#  uni[length(uni)] <- quantile(s, probs = mon)
#  return(data.frame(X = uni, Empirical_CDF = cdf))
#}

emp.cdf <- function(data, cor = 0.25){
  s <- sort(data); n <- length(data)
  uni <- unique(sort(data))
  cdf <- array(NA, dim = length(uni))
  for(i in uni){
   pos <- match(i, uni)
   cdf[pos] <- length(which(s <= i))/(n + cor )
  }
  return(data.frame(X = uni, Empirical_CDF = cdf))
}

QQ.plot.y <- function(data, alpha.est, sigma.est, corr = 0.25){
  jj <- emp.cdf(data, cor = corr)
  mm <- array(NA, dim = length(jj$X))
  for(i in jj$Empirical_CDF){
    pos <- match(i, jj$Empirical_CDF)
    mm[pos] <- NP.quantile(u = i, alpha = alpha.est, sigma = sigma.est)
  }
  plot(mm, jj$X, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
   #    main = "Normal-Pareto Q-Q Plot")
  abline(a = 0, b = 1, col = "red")
}
















CNBD_quantile <- function(r.par = 1, beta.par = 2, u = 0.1){
  require(rootSolve)
  ans.roots <- array(NA, dim = length(u))
  for(i in u){
    pos <- match(i, u)
    func <- function(x, t){
      h <- function(t) ((beta.par^r.par)/gamma(r.par))*pgamma(t, x)*t^(r.par - 1)*exp(- beta.par*t)
      f <-  1 - integrate(h, lower = 0, upper = Inf)$value - i
      return(f)
    }
    ans.roots[pos] <- uniroot(func, lower =  0.0002, upper = 100, extendInt = "yes")$root
  }
  return(ans.roots)
}

###---------------- QQ Plot for CNBD -------------------#####
QQPlot_CNBD <- function(data, r, beta, corr = 0.25){
  emp.cdf <- function(dat = data, cor = corr){
    s <- sort(dat); n <- length(dat)
    uni <- unique(sort(dat))
    cdf <- array(NA, dim = length(uni))
    for(i in uni){
      pos <- match(i, uni)
      cdf[pos] <- length(which(s <= i))/(n + cor)
    }
    return(data.frame(X = uni, Empirical_CDF = cdf))
  }

  jj <- emp.cdf( )

  mm <- array(NA, dim = length(jj$X))
  for(i in jj$Empirical_CDF){
    pos <- match(i, jj$Empirical_CDF)
    mm[pos] <- CNBD_quantile(r.par = r, beta.par = beta, u = i)
  }
  plot(mm, jj$X, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
       main = "CNBD Q-Q Plot")
  abline(a = 0, b = 1, col = "red")
}

# Example showing QQ plot works
jj <- CNBD_rand(n = 100, r = 2, beta = 1)
QQPlot_CNBD(data = jj, r = 2, beta = 1)


###------ Estimation procedures -------###
##### MME
lmix.MME <- function(obser.y){
  require(rootSolve)
  f.alpha <- function(alpha){
    C.1 <- 0.5*log(2) + 0.5*digamma(0.5) + 0.5*digamma(1)
    C.2 <- (gamma(3/4)*gamma(1/4))/(sqrt(pi)*(2^0.25))
    res <- 2*lgamma((1/alpha) + 0.25) - 2*lgamma(1/alpha) -
      0.5*digamma(1/alpha) -
      mean(log(abs(obser.y))) + C.1 +
      2*log(C.2/mean(abs(obser.y)^(-0.5)))
    return(res)
  }
  R.fun <- function(cc){
    C.1 <- 0.5*(log(2) + digamma(1/2) + digamma(1))
    C.2 <- (gamma(0.75)*gamma(0.25))/(sqrt(pi)*(2^0.25))
    ves <- mean(log(abs(cc))) - C.1 -
      2*log( C.2/mean(abs(cc)^(-0.5)))
    return(ves)
  }

  if(R.fun(obser.y) > 0){
    out.alpha <- uniroot(f.alpha, lower = 0.01,
                         upper = 100, extendInt = 'yes')$root
    C.2 <- (gamma(3/4)*gamma(1/4))/(sqrt(pi)*(2^0.25))
    out.sigma <- exp(2*(log(((out.alpha^(0.25))/mean(abs(obser.y)^(-0.5)))*C.2) +
                          lgamma( (1/out.alpha) + 0.25) - lgamma(1/out.alpha)))
    out <- c(alpha = out.alpha, sigma = out.sigma, R.value = R.fun(obser.y))
    return(out)
  }
  else{
    out.alpha <- 0
    C.2 <- (gamma(3/4)*gamma(1/4))/(sqrt(pi)*(2^0.25))
    out.sigma <-  (C.2/mean((abs(obser.y)^(-0.5))))^2
    out <- c(alpha = out.alpha, sigma = out.sigma, R.value = R.fun(obser.y))
    return(out)
  }
}


### EM Algorithm
