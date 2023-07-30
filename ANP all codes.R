library(cubature); library(rootSolve)
####---- The PDF -------#####
ANP.pdf <- function(y, alpha = 0.5, sigma = 1, mu = 1){
 bb <- array(NA, dim = length(y))
  for(i in y){
    pos <- match(i, y)
    if(i >= 0){
    gg <- function(t){
      yy <- t*exp(-2*abs(i)*t/(sqrt(2*sigma^2*t + mu^2) + mu))/sqrt(2*sigma^2*t + mu^2)
      uu <- yy*(t^(1/alpha - 1)*exp(-t/alpha))/(alpha^(1/alpha)*gamma(1/alpha))
      return(uu)
    }
    #bb[pos] <- integrate(gg, lower = 0, upper = Inf)$value
    bb[pos] <- adaptIntegrate(gg, lowerLimit = 0, upperLimit = Inf)$integral
    }
  else{
    gg <- function(t){
      yy <- t*exp(-2*abs(i)*t/(sqrt(2*sigma^2*t + mu^2) - mu))/sqrt(2*sigma^2*t + mu^2)
      uu <- yy*(t^(1/alpha - 1)*exp(-t/alpha))/(alpha^(1/alpha)*gamma(1/alpha))
      return(uu)
    }
    bb[pos] <- integrate(gg, lower = 0, upper = Inf)$value
  }
  }
return(bb)
}

###----- Curve for PDF -----#####
graph.f.x.ANP <- function(alpha.par = 0.5, sigma.par = 1, mu.par = 1, start = -10, end = 10){
  ANP.pdf <- function(y, alpha = alpha.par, sigma = sigma.par, mu = mu.par){
    bb <- array(NA, dim = length(y))
    for(i in y){
      pos <- match(i, y)
      if(i >= 0){
        gg <- function(t){
          yy <- t*exp(-2*abs(i)*t/(sqrt(2*sigma^2*t + mu^2) + mu))/sqrt(2*sigma^2*t + mu^2)
          uu <- yy*(t^(1/alpha - 1)*exp(-t/alpha))/(alpha^(1/alpha)*gamma(1/alpha))
          return(uu)
        }

        bb[pos] <- integrate(gg, lower = 0, upper = Inf)$value
      }
      else{
        gg <- function(t){
          yy <- t*exp(-2*abs(i)*t/(sqrt(2*sigma^2*t + mu^2) - mu))/sqrt(2*sigma^2*t + mu^2)
          uu <- yy*(t^(1/alpha - 1)*exp(-t/alpha))/(alpha^(1/alpha)*gamma(1/alpha))
          return(uu)
        }

        bb[pos] <- integrate(gg, lower = 0, upper = Inf)$value
      }

    }
    return(bb)
  }
  curve(ANP.pdf, from = start, to = end, ylab = "f(y)", xlab = "y", col = "blue")
  legend(x = 3, y = 0.3, legend = paste(c("alpha =", alpha.par, "; sigma =", sigma.par, "; mu =", mu.par),
                                        collapse = " "), lty = 1, col =c("blue"), cex = 0.7, box.col = "white", text.width = 3)
}

graph.f.x.ANP(alpha.par = 0.5, sigma.par = 1)

### PDF plots
val <- seq(-10, 10, by = 0.01)
f.1 <- ANP.pdf(y = val, alpha = 0.7, sigma = 1, mu = 1)
f.2 <- ANP.pdf(y = val, alpha = 2, sigma = 1, mu = 1)
f.3 <- ANP.pdf(y = val, alpha = 5, sigma = 1, mu = 1)

plot(val, f.1, type = "l", col = "black", xlab = "y", ylab = "Density", ylim = c(0, max(f.1, f.2, f.3)))
points(val, f.2, type = "l", lty = "dotted", col = "red")
points(val, f.3, type = "l", lty = "twodash", col = "blue")
legend(x = -9, y = 0.3, legend = c(expression(paste(mu, " = 1; ", alpha, " = 0.7; ", sigma, " = 1")), expression(paste(mu, " = 1; ", alpha, " = 2; ", sigma, " = 1")), expression(paste(mu, " = 1; ", alpha, " = 5; ", sigma, " = 1"))),
       lty = c(1, 3, 5), col = c("black", "red", "blue"), cex = 1.2, text.width = 2.5, bty = "n")

g.1 <- ANP.pdf(y = val, alpha = 1.5, sigma = 1, mu = -1)
g.2 <- ANP.pdf(y = val, alpha = 1.5, sigma = 1, mu = 0)
g.3 <- ANP.pdf(y = val, alpha = 1.5, sigma = 1, mu = 2)

par(mfrow = c(1, 3))
plot(val, g.1, type = "l", xlab = "y", ylab = "Density")
legend(x = -2, y = 0.3, legend = c(expression(paste(mu, " = -1; ", alpha, " = 1.5; ", sigma, " = 1"))),
       cex = 1.2, text.width = 2.5, bty = "n")

plot(val, g.2, type = "l", xlab = "y", ylab = "Density")
legend(x = 0, y = 0.3, legend = c(expression(paste(mu, " = 0; ", alpha, " = 1.5; ", sigma, " = 1"))),
       cex = 1.2, text.width = 2.5, bty = "n")

plot(val, g.3, type = "l", xlab = "y", ylab = "Density")
legend(x = 0, y = 0.2, legend = c(expression(paste(mu, " = 2; ", alpha, " = 1.5; ", sigma, " = 1"))),
       cex = 1.2, text.width = 2.5, bty = "n")

###---- The CDF ------#####
ANP.cdf <- function(y, alpha = 1, sigma = 1, mu = 1){
 bb <- array(NA, dim = length(y))
  for(i in y){
   pos <- match(i, y)
#    gg <- function(t){(0.5*exp(-sqrt(2*t)*abs(i)/sigma)*t^((1/alpha) - 1)*exp(-t/alpha))/(alpha^(1/alpha)*gamma(1/alpha))}
  gg <- function(t){0.5*((1 + (mu/(sqrt(2*(sigma^2)*t + mu^2))))*exp(-2*t*abs(i)/( sqrt(2*(sigma^2)*t + mu^2 ) + mu))*(t^(1/alpha - 1)*exp(-t/alpha)/(alpha^(1/alpha)*gamma(1/alpha))))}
  gg.1 <- function(t){0.5*((1 - (mu/(sqrt(2*(sigma^2)*t + mu^2))))*exp(-2*t*abs(i)/( sqrt(2*(sigma^2)*t + mu^2 ) - mu))*(t^(1/alpha - 1)*exp(-t/alpha)/(alpha^(1/alpha)*gamma(1/alpha))))}

  if(i < 0){
      #bb[pos] <- integrate(gg.1, lower = 0, upper = Inf)$value
    bb[pos] <- adaptIntegrate(gg.1, lowerLimit = 0, upperLimit = Inf)$integral
      }
    else if(i >= 0){
    #  bb[pos] <- 1 - integrate(gg, lower = 0, upper = Inf)$value
      bb[pos] <- 1 - adaptIntegrate(gg, lowerLimit = 0, upperLimit = Inf)$integral
      }
    }
 return(bb)
}

####----- Curve for CDF -----####
curve(ANP.cdf, from = -30, to = 30, ylab = "F(y)", xlab = "y")

### Other CDFs with different parameter combinations
val <- seq(-30, 30, by = 0.01)
f.1 <- ANP.cdf(y = val, alpha = 0.7, sigma = 1, mu = 1)
f.2 <- ANP.cdf(y = val, alpha = 1, sigma = 1, mu = 1)
f.3 <- ANP.cdf(y = val, alpha = 2, sigma = 1, mu = 1)

plot(val, f.1, type = "l", col = "black", xlab = "y", ylab = "CDF", ylim = c(0, max(f.1, f.2, f.3)))
points(val, f.2, type = "l", lty = "dotted", col = "red")
points(val, f.3, type = "l", lty = "twodash", col = "blue")
legend(x = -35, y = 0.4, legend = c(expression(paste(mu, " = 1; ", alpha, " = 0.7; ", sigma, " = 1")), expression(paste(mu, " = 1; ", alpha, " = 1; ", sigma, " = 1")), expression(paste(mu, " = 1; ", alpha, " = 2; ", sigma, " = 1"))),
       lty = c(1, 3, 5), col = c("black", "red", "blue"), cex = 1.2, text.width = 1, bty = "n", seg.len = 0.5)

val <- seq(-30, 30, by = 0.01)
g.1 <- ANP.cdf(y = val, alpha = 1, sigma = 1, mu = -1)
g.2 <- ANP.cdf(y = val, alpha = 1, sigma = 1, mu = 0)
g.3 <- ANP.cdf(y = val, alpha = 1, sigma = 1, mu = 2)

plot(val, g.1, type = "l", col = "black", xlab = "y", ylab = "CDF", ylim = c(0, max(g.1, g.2, g.3)))
points(val, g.2, type = "l", lty = "dotted", col = "red")
points(val, g.3, type = "l", lty = "twodash", col = "blue")
legend(x = -35, y = 0.6, legend = c(expression(paste(mu, " = -1; ", alpha, " = 1; ", sigma, " = 1")), expression(paste(mu, " = 0; ", alpha, " = 1; ", sigma, " = 1")), expression(paste(mu, " = 2; ", alpha, " = 1; ", sigma, " = 1"))),
       lty = c(1, 3, 5), col = c("black", "red", "blue"), cex = 1.2, text.width = 1, bty = "n", seg.len = 0.5)

#plot(val, f.1, type = "l", col = "black", xlab = "y", ylab = "F(y)", ylim = c(min(f.1,f.2,f.3), 1))
#points(val, f.2, type = "l", lty = "dotted", col = "red")
#points(val, f.3, type = "l", lty = "twodash", col = "blue")
#legend("bottomright", legend = c("alpha = 0.6, sigma = 1, mu = -1.2", "alpha = 0.6, sigma = 1, mu = 0",
#                                    "alpha = 0.6, sigma = 1, mu = 3"), lty = c(1, 3, 5), col = c("black", "red", "blue"), cex = 0.7, text.width = 8)
#legend(x = 0.5, y = 0.4, legend = c("alpha = 0.6, sigma = 3", "alpha = 1, sigma = 1",
#                                    "alpha = 2, sigma = 1.5"), lty = c(1, 3, 5), col = c("black", "red", "blue"), cex = 0.7, text.width = 8)

# Obtain n random samples from NP
rANPare <- function(n = 1, alpha = 0.5, sigma = 1, mu = 1){
    normal <- rnorm(n); exponential <- rexp(n)
    gamm <- rgamma(n, shape = 1/alpha, scale = alpha)
    y <- mu*(exponential/gamm) + sigma*sqrt(exponential/gamm)*normal
    return(y)
}

# Graph overlay of numbers generated versus the density
par(mfrow = c(2, 2))

alpha.par <- 0.5; sigma.par <- 1; mu.par = 1
X <- rANPare(n = 1000, alpha = alpha.par, sigma = sigma.par, mu = mu.par)
jj <- ANP.pdf(y = seq(min(X), max(X), by = 0.001), alpha = alpha.par, sigma = sigma.par, mu = mu.par)
hist(X, freq = FALSE, ylim = c(0, max(jj)), main = "", col = "gray", xlab = "y", breaks = 80)
points(x = seq(min(X), max(X), by = 0.001), y = jj, type = "l", col = "red")
legend(x = 7, y = 0.3, legend = c(expression(paste(alpha, " = 0.5; ", sigma, " = 1; ",mu, " = 1" ))), bty = "n", cex = 1.2, text.width = 6.5)

#legend(x = 9, y = 0.3, legend = paste(c("alpha =", alpha.par, ";", "sigma =", sigma.par,";","mu =", mu.par), collapse = " "),box.col = "white", cex = 0.8, text.width = 6.5)
#alpha.par <- 0.7; sigma.par <- 1; mu.par = -1
#X <- rANPare(n = 1000, alpha = alpha.par, sigma = sigma.par, mu = mu.par)
#jj <- ANP.pdf(y = seq(min(X), max(X), by = 0.01), alpha = alpha.par, sigma = sigma.par, mu = mu.par)
#hist(X, freq = FALSE, ylim = c(0, max(jj)), main = "", col = "gray", xlab = "y", breaks = 80)
#points(x = seq(min(X), max(X), by = 0.01), y = jj, type = "l", col = "red")
#legend(x = "center", legend = paste(c("alpha =", alpha.par, ";", "sigma =", sigma.par,";","mu =", mu.par), collapse = " "),box.col = "white", cex = 0.8, text.width = 6.5)
#legend(x = "top", legend = paste(c("alpha =", alpha.par, ";", "sigma =", sigma.par), collapse = " "),box.col = "white", cex = 0.8, text.width = 6.5)

alpha.par <- 0.5; sigma.par <- 1; mu.par = - 1
X <- rANPare(n = 1000, alpha = alpha.par, sigma = sigma.par, mu =  mu.par)
jj <- ANP.pdf(y = seq(min(X), max(X), by = 0.01), alpha = alpha.par, sigma = sigma.par, mu = mu.par)
hist(X, freq = FALSE, ylim = c(0, max(jj)), main = "", col = "gray", xlab = "y", breaks = 80)
points(x = seq(min(X), max(X), by = 0.01), y = jj, type = "l", col = "red")
legend(x = -60, y = 0.3, legend = c(expression(paste(alpha, " = 0.5; ", sigma, " = 1; ",mu, " = - 1" ))), bty = "n", cex = 1.2, text.width = 6.5)
#legend(x = "topleft", legend = paste(c("alpha =", alpha.par, ";", "sigma =", sigma.par,"; mu =", mu.par), collapse = " "),box.col = "white", cex = 0.8, text.width = 6.5)

alpha.par <- 0.5; sigma.par <- 1; mu.par = 0
X <- rANPare(n = 1000, alpha = alpha.par, sigma = sigma.par, mu =  mu.par)
jj <- ANP.pdf(y = seq(min(X), max(X), by = 0.01), alpha = alpha.par, sigma = sigma.par, mu = mu.par)
hist(X, freq = FALSE, ylim = c(0, max(jj)), main = "", col = "gray", xlab = "y", breaks = 80)
points(x = seq(min(X), max(X), by = 0.01), y = jj, type = "l", col = "red")
legend(x = -2, y = 0.4, legend = c(expression(paste(alpha, " = 0.5; ", sigma, " = 1; ",mu, " = 0" ))), bty = "n", cex = 1.2, text.width = 6.5)

#alpha.par <- 0.1; sigma.par <- 10; mu.par = 5
#X <- rANPare(n = 1000, alpha = alpha.par, sigma = sigma.par, mu = mu.par)
#jj <- ANP.pdf(y = seq(min(X), max(X), by = 0.01), alpha = alpha.par, sigma = sigma.par, mu = mu.par)
#hist(X, freq = FALSE, ylim = c(0, max(jj)+0.02), main = "", col = "gray", xlab = "y", breaks = 45)
#points(x = seq(min(X), max(X), by = 0.01), y = jj, type = "l", col = "red")
#legend(x = "top", legend = paste(c("alpha =", alpha.par, ";", "sigma =", sigma.par, "; mu =", mu.par), collapse = " "),box.col = "white", cex = 0.8, text.width = 6.5)

par(mfrow = c(1, 1))

#####------------ Quantile Function -----------------######
AY.quantile <- function(alpha = 0.5, sigma = 1, mu = 1, u = 0.1){
  require(rootSolve)
  xxx <- ANP.cdf(alpha = alpha, sigma = sigma, mu = mu, y = 0)
  ans.roots <- array(NA, dim = length(u))

  for(i in u){
  pos <- match(i, u)
  #i = 0.4177055
 func <- function(x, t){
  #h <- function(t) ((beta.par^r.par)/gamma(r.par))*pgamma(t, x)*t^(r.par - 1)*exp(- beta.par*t)
  h1 <- function(t)(0.5*((1 + mu/sqrt(2*sigma^2*t + mu^2))*exp(-2*t*abs(x)/(sqrt(2*sigma^2*t + mu^2) + mu))*(t^(1/alpha - 1)*exp(-t/alpha)/(alpha^(1/alpha)*gamma(1/alpha)))))
  f <-  1 - integrate(h1, lower = 0, upper = Inf)$value - i
  return(f)
  }

 func2 <- function(x, t){
   #h <- function(t) ((beta.par^r.par)/gamma(r.par))*pgamma(t, x)*t^(r.par - 1)*exp(- beta.par*t)
   h2 <- function(t)(0.5*((1 - mu/sqrt(2*sigma^2*t + mu^2))*exp(-2*t*abs(x)/(sqrt(2*sigma^2*t + mu^2) - mu))*(t^(1/alpha - 1)*exp(-t/alpha)/(alpha^(1/alpha)*gamma(1/alpha)))))
   f <-  integrate(h2, lower = 0, upper = Inf)$value - i
   return(f)
 }

 if(i == xxx){
   ans.roots[pos] <- 0
 }
 else if(i < xxx){
   ans.roots[pos] <- uniroot(func2, lower = -100, upper = -0.0002, extendInt = "yes")$root
 }
 else{
   ans.roots[pos] <- uniroot(func, lower =  0.0002, upper = 100, extendInt = "yes")$root
 }
  }
 return(ans.roots)
}

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

QQ.plot.Ay <- function(data, alpha.est, sigma.est, mu.est, corr = 0.25){
  jj <- emp.cdf(data, cor = corr)
  mm <- array(NA, dim = length(jj$X))
  for(i in jj$Empirical_CDF){
    pos <- match(i, jj$Empirical_CDF)
    mm[pos] <- AY.quantile(i, alpha = alpha.est, sigma = sigma.est, mu = mu.est)
  }
  #plot(mm, jj$X, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
   #    main = "Asymmetric Normal-Pareto Q-Q Plot")
  plot(mm, jj$X, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
  abline(a = 0, b = 1, col = "red")
}

# Example showing QQ plot works
jj <- rANPare(n = 1000, alpha = 0.5, sigma = 1, mu = 2)
QQ.plot.Ay(data = jj, alpha.est = 0.5, sigma.est = 1, mu.est = 2)


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









#ANP.quantile <- function(alpha = 0.5, sigma = 1, mu = 1, u = 0.1){
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
#   pos <- match(i, u)

#    if(i == 0.5){
#      ahy[pos] <- 0
#    }
#   else if(i < 0.5){
#      ahy[pos] <- return(uniroot(g.y.1, lower = -15, upper = -0.00000001,
#                                extendInt = "yes", maxiter = 10000 )$root)
#    }
#   else if(i > 0.5){
#      ahy[pos] <- return(uniroot(g.y.2, lower = 0.00000001, upper = 15,
#                               extendInt = "yes", maxiter = 10000 )$root)
#   }
#  }
#  return(ahy)
#}








## testing whether inverse function is working well
#rrr <- seq(0.1, 0.9, by = 0.001)
#sss <- AY.quantile(u = rrr)
#qqq <- ANP.cdf(y = sss, alpha = 0.5)












#### Other coding stuff

#for(i in u){
##  pos <- match(i, u)
#  func <- function(x, t){
#    h <- function(t) ((beta.par^r.par)/gamma(r.par))*pgamma(t, x)*t^(r.par - 1)*exp(- beta.par*t)
#    f <-  1 - integrate(h, lower = 0, upper = Inf)$value - i
#    return(f)
#  }

#}
#sam <- rgamma(10000, shape = 1/alpha, scale = alpha)

#g.y.1 <- function(y){
#  fin <- mean( 0.5*exp(-sqrt(2*sam)*abs(y))) - u
#  return(fin)}
#g.y.2 <- function(y){
#  fin <- 1 - u - mean(0.5*exp(-sqrt(2*sam)*abs(y)))
#  return(fin)}
#require(rootSolve)

#ahy <- array(NA, dim = length(u))

#for(i in u){
#  pos <- match(i, u)

#  if(i == xxx){
 #   ahy[pos] <- 0
#  }
#  else if(i < xxx){
#    ahy[pos] <- return(uniroot(g.y.1, lower = -15, upper = -0.00000001,
#                               extendInt = "yes", maxiter = 10000 )$root)
#  }
#  else if(i > xxx){
#    ahy[pos] <- return(uniroot(g.y.2, lower = 0.00000001, upper = 15,
#                               extendInt = "yes", maxiter = 10000 )$root)
#  }
#}
#return(ahy)
#}

#require(rootSolve)
#ans.roots <- array(NA, dim = length(u))
#for(i in u){
#  pos <- match(i, u)
#  func <- function(x, t){
#    h <- function(t) ((beta.par^r.par)/gamma(r.par))*pgamma(t, x)*t^(r.par - 1)*exp(- beta.par*t)
#    f <-  1 - integrate(h, lower = 0, upper = Inf)$value - i
#    return(f)
#  }
#  ans.roots[pos] <- uniroot(func, lower =  0.0002, upper = 100, extendInt = "yes")$root
#}
#return(ans.roots)
#}


