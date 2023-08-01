pcpois <- function(q, lambda, lower.tail = TRUE){
  if(!is.numeric(q) | !is.numeric(lambda)) stop("inputs x and lambda must be numeric")
  if(lambda <= 0) stop("lambda must be positive")

  suppressWarnings(df <- pgamma(lambda, shape = q, lower.tail = !lower.tail))
  df[which(is.na(df) == TRUE)] = 0
  return(df)
}

# Poisson discrete
disP <- function(k, lambda){
  return(pcpois(k, lambda = lambda) - pcpois(k-1, lambda = lambda))
}

# Inverse Poisson discrete
disIP <- function(k, lambda){
  #if(k == 1){
#  return(1 - pcpois(1/k, lambda = lambda))
#  }
#  else{
    return(pcpois(1/(k - 1), lambda = lambda) - pcpois(1/k, lambda = lambda))
#  }
}

disP <- function(k, lambda){
  return(pcpois(k+1, lambda = lambda) - pcpois(k, lambda = lambda))
}

disPcum <- function(x, lambda){
  ar <- array(NA, dim = length(x))
  for(i in x){
  ar[match(i, x)] <- sum(disP(1:i, lambda = lambda))
  }
  return(ar)
 # return(sum(disP(1:x, lambda = lambda)))
}

disIPcum <- function(x, lambda){
  ar <- array(NA, dim = length(x))
  for(i in x){
    ar[match(i, x)] <- sum(disIP(1:i, lambda = lambda))
  }
  return(ar)
  # return(sum(disP(1:x, lambda = lambda)))
}

plot(seq(1, 10), disP(1:10, 2), type = "h", col = "red")
points(seq(1, 10), disP(1:10, 2), type = "o", col = "red")
points(seq(1, 10), disIP(1:10, 2), type = "o", col = "blue")

plot(seq(1, 10), disP(1:10, 2), pch = 19, col = "red", ylim = c(0, max(disP(1:10, 2), disIP(1:10, 2))), xlab = "n", ylab = "P(n)")
points(seq(1, 10), disP(1:10, 2), col = "black", type = "h", lty = 3)
points(seq(1, 10), disIP(1:10, 2), col = "blue")
points(seq(1, 10), disIP(1:10, 2), col = "black", type = "h", lty = 3)
axis(1, at = seq(0:10))
#abline(v = seq(0, 10), col = "black", lty = 3)
legend(5, 0.85, pch = c(19, 1), col = c("red", "blue"),
       legend = c(expression(paste("Poisson (", lambda, " = ", 2, ")")), expression(paste("Inverse Poisson (", lambda, " = ",2, ")"))), bty ="n", cex = 1.2)

plot(seq(1, 15), disPcum(1:15, 2), type = "s", col = "red", xlab = "n", ylab = "F(n)")
points(seq(1, 15), disPcum(1:15, 2), type = "s", col = "red")
points(seq(1, 15), disIPcum(1:15, 2), type = "s", col = "blue")
axis(1, at = seq(0:15))
abline(v = seq(0, 15), col = "grey", lty = 3)
abline(h = c(0, 0.2, 0.4, 0.6, 0.8, 1), col = "grey", lty = 3)
legend(7, 0.6, col = c("red", "blue"), legend = c("Poisson", "Inverse Poisson"), lty = c(1, 1), bty = "n", cex = 1.2, seg.len = 0.5)
par(mfrow = c(1, 2))
grid(ny = NA, col = "black", lty = 3)

points(seq(1, 10), disP(1:10, 2), type = "o", col = "red")
axis(1, at = seq(0:10))
grid()
points(seq(1, 10), disIP(1:10, 2), type = "o", col = "blue")

# Quick check of distribution for random mean
yy <- function(p, sims ){
  rr <- rgeom(n = sims, prob = p) + 1
  jj <- array(NA, dim = sims)

  for(i in rr){
  pos <- match(i, rr)
  jj[pos] <- (sum(runif(n = i, min = -1, max = 1))/i)/( (1/3)/sqrt(i) )
  }
  hist(jj, freq = FALSE, ylim = c(0 , dnorm(0)));points(seq(-3, 3, 0.01), dnorm(seq(-3, 3, 0.01)), type = "l" )
}

hist(rnorm(100))
points(seq(-2, 2, 0.01),seq(-2, 2, 0.01), type = "l" )

# Geometric discrete
disG <- function(k, p){
  lambda = -log(1 - p)
  return(pexp(k+1, rate = lambda) - pexp(k, rate = lambda))
}

# Inverse Geometric discrete
disIG <- function(k, p){
  lambda = -log(1 - p)
  jj <- array(NA, dim = length(k))
  for(i in k){
  pos <- match(i, k)
  jj[pos] <- exp(-(1*lambda)/(i + 1)) - exp(-(1*lambda)/i)
  }
  return(jj)
}

disGcum <- function(x, p){
  ar <- array(NA, dim = length(x))
  for(i in x){
    ar[match(i, x)] <- sum(disG(0:i, p = p))
  }
  return(ar)
}

disIGcum <- function(x, p){
  ar <- array(NA, dim = length(x))
  for(i in x){
    ar[match(i, x)] <- sum(disIG(0:i, p = p))
  }
  return(ar)
}

plot(seq(1, 10), disG(1:10, 1/3), type = "h", col = "red")
points(seq(1, 10), disG(1:10, 1/3), type = "o", col = "red")
points(seq(1, 10), disIG(1:10, 1/3), type = "o", col = "blue")

plot(seq(0, 10), disG(0:10, p=1/3), pch = 19, col = "red", ylim = c(0, max(disG(0:10, 1/3), disIG(0:10, 1/3))), xlab = "n", ylab = "P(n)")
points(seq(0, 10), disG(0:10, 1/3), col = "black", type = "h", lty = 3)
points(seq(0, 10), disIG(0:10, 1/3), col = "blue", pch = 19)
points(seq(0, 10), disIG(0:10, 1/3), col = "black", type = "h", lty = 3)
axis(1, at = seq(0:10))
#abline(v = seq(0, 10), col = "black", lty = 3)
legend(3, 0.5, pch = c(19, 19), col = c("red", "blue"),
       legend = c(expression(paste("Geometric (", p, " = ", 1/3, ")")), expression(paste("Inverse Geometric (", p, " = ",1/3, ")"))), bty ="n", cex = 1.2)

plot(seq(0, 15), disGcum(0:15, 1/3), type = "s", col = "red", xlab = "n", ylab = "F(n)")
points(seq(0, 15), disGcum(0:15, 1/3), type = "s", col = "red")
points(seq(0, 15), disIGcum(0:15, 1/3), type = "s", col = "blue")
axis(1, at = seq(0:15))
abline(v = seq(0, 15), col = "grey", lty = 3)
abline(h = c(0, 0.2, 0.4, 0.6, 0.8, 1), col = "grey", lty = 3)
legend(4, 0.5, col = c("red", "blue"), legend = c("Geometric", "Inverse Geometric"), lty = c(1, 1), bty = "n", cex = 1.2, seg.len = 0.5)

# Geometric discrete
disG <- function(k){
  ar <- array(NA, dim = length(k))
  for(i in k){
    ar[match(i, k)] <- (exp(-1))^(i - 1)*(1 - exp(-1))
  }
  return(ar)}

# Inverse Geometric discrete
disIG <- function(k){
  jj <- array(NA, dim = length(k))
  for(i in k){
    pos <- match(i, k)
    jj[pos] <- exp(-1/i) - exp(-1/(i - 1))
  }
  return(jj)}

disGcum <- function(x){
  ar <- array(NA, dim = length(x))
  for(i in x){
    ar[match(i, x)] <- sum(disG(1:i))
  }
  return(ar)}

disIGcum <- function(x){
  ar <- array(NA, dim = length(x))
  for(i in x){
    ar[match(i, x)] <- sum(disIG(1:i))
  }
  return(ar)}

#plot(seq(1, 10), disG(1:10), type = "h", col = "red")
#points(seq(1, 10), disG(1:10), type = "o", col = "red")
#points(seq(1, 10), disIG(1:10), type = "o", col = "blue")

plot(seq(1, 10), disG(1:10), pch = 19, col = "red", ylim = c(0, max(disG(1:10), disIG(1:10))), xlab = "n", ylab = "P(n)")
points(seq(1, 10), disG(1:10), col = "black", type = "h", lty = 3)
points(seq(1, 10), disIG(1:10), col = "blue", pch = 19)
points(seq(1, 10), disIG(1:10), col = "black", type = "h", lty = 3)
axis(1, at = seq(1:10))
#abline(v = seq(0, 10), col = "black", lty = 3)
legend(5, 0.5, pch = c(19, 19), col = c("red", "blue"),
       legend = c("Geometric", "Inverse Geometric"), bty ="n", cex = 1.2)

plot(seq(1, 10), disGcum(1:10), ylim = c(0, max(disGcum(1:10), disIGcum(1:10))), type = "s", col = "red", xlab = "n", ylab = "F(n)")
#points(seq(1, 15), disGcum(1:15), type = "s", col = "red")
points(seq(1, 10), disIGcum(1:10), type = "s", col = "blue")
axis(1, at = seq(1:10))
abline(v = seq(0, 10), col = "grey", lty = 3)
#abline(h = c(0, 0.2, 0.4, 0.6, 0.8, 1), col = "grey", lty = 3)
legend(4, 0.5, col = c("red", "blue"), legend = c("Geometric", "Inverse Geometric"), lty = c(1, 1), bty = "n", cex = 1.2, seg.len = 0.5)

library(expint)

####### discrete gamma and inverse gamma
disGam <- function(k, beta, alpha){
  ar <- array(NA, dim = length(k))
  for(i in k){
#   ar[match(i, k)] <- pgamma(beta*(i-i), shape = alpha, lower.tail = F) - pgamma(beta*i, shape = alpha, lower.tail = F)
    ar[match(i, k)] <- pgamma(i, shape = alpha, rate = beta) - pgamma(i-1, shape = alpha, rate = beta)
  }
  return(ar)}

# Inverse Geometric discrete
disIGam <- function(k, beta, alpha){
  jj <- array(NA, dim = length(k))
  for(i in k){
    pos <- match(i, k)
    jj[pos] <- gammainc(alpha, beta/i)/gamma(alpha) - (gammainc(alpha, beta/(i - 1))/gamma(alpha))
  }
  return(jj)}

disGamcum <- function(x, beta, alpha){
  ar <- array(NA, dim = length(x))
  for(i in x){
    ar[match(i, x)] <- suppressWarnings(sum(disGam(k = 1:i, beta = beta, alpha = alpha)))
  }
  return(ar)}

disIGamcum <- function(x, beta, alpha){
  ar <- array(NA, dim = length(x))
  for(i in x){
    ar[match(i, x)] <-suppressWarnings(sum(disIGam(1:i, beta = beta, alpha = alpha)))
  }
  return(ar)}

#plot(seq(1, 10), disG(1:10), type = "h", col = "red")
#points(seq(1, 10), disG(1:10), type = "o", col = "red")
#points(seq(1, 10), disIG(1:10), type = "o", col = "blue")

plot(seq(1, 10), disGam(1:10, alpha = 2, beta = 1), pch = 19, col = "red", ylim = c(0, max(disGam(1:10, alpha = 2, beta = 1), disIGam(1:10, alpha = 2, beta = 1))), xlab = "n", ylab = "P(n)")
points(seq(1, 10), disGam(1:10, alpha = 2, beta = 1), col = "black", type = "h", lty = 3)
points(seq(1, 10), disIGam(1:10, alpha = 2, beta = 1), col = "blue", pch = 19)
points(seq(1, 10), disIGam(1:10, alpha = 2, beta = 1), col = "black", type = "h", lty = 3)
axis(1, at = seq(1:10))
#abline(v = seq(0, 10), col = "black", lty = 3)
legend(3, 0.6, pch = c(19, 19), col = c("red", "blue"),
       legend = c(expression(paste("Discrete Gamma(", alpha, " = 2, ",beta," = 1)")), expression(paste("Inverse discrete Gamma(", alpha, " = 2, ",beta," = 1)" ))), bty = "n", cex = 1.2)
#       legend = c("Geometric", "Inverse Geometric"), bty ="n", cex = 1.2)

plot(seq(1, 15), disGamcum(1:15, alpha = 2, beta = 1), ylim = c(0, max(disGamcum(1:15, alpha = 2, beta = 1), disIGamcum(1:15, alpha = 2, beta = 1))), type = "s", col = "red", xlab = "n", ylab = "F(n)")
#points(seq(1, 15), disGcum(1:15), type = "s", col = "red")
points(seq(1, 15), disIGamcum(1:15, alpha = 2, beta = 1), type = "s", col = "blue")
axis(1, at = seq(1:15))
abline(v = seq(0, 10), col = "grey", lty = 3)
#abline(h = c(0, 0.2, 0.4, 0.6, 0.8, 1), col = "grey", lty = 3)
legend(5, 0.5, col = c("red", "blue"), legend = c(expression(paste("Discrete Gamma(", alpha, " = 2, ",beta," = 1)")), expression(paste("Inverse discrete Gamma(", alpha, " = 2, ",beta," = 1)" ))), lty = c(1, 1), bty = "n", cex = 1.2, seg.len = 0.3)

# Weibull and inverse Weibull
disWei <- function(k, alpha){
  ar <- array(NA, dim = length(k))
  for(i in k){
    #   ar[match(i, k)] <- pgamma(beta*(i-i), shape = alpha, lower.tail = F) - pgamma(beta*i, shape = alpha, lower.tail = F)
    ar[match(i, k)] <- exp(-(i - 1)^alpha) - exp(-i^alpha)
  }
  return(ar)}

# Inverse Geometric discrete
disIWei <- function(k, alpha){
  jj <- array(NA, dim = length(k))
  for(i in k){
    pos <- match(i, k)
    jj[pos] <- exp(-1/(i^alpha)) - exp(-1/(i - 1)^alpha)
  }
  return(jj)}

disWeicum <- function(x, alpha){
  ar <- array(NA, dim = length(x))
  for(i in x){
    ar[match(i, x)] <- suppressWarnings(sum(disWei(k = 1:i, alpha = alpha)))
  }
  return(ar)}

disIWeicum <- function(x, alpha){
  ar <- array(NA, dim = length(x))
  for(i in x){
    ar[match(i, x)] <-suppressWarnings(sum(disIWei(1:i, alpha = alpha)))
  }
  return(ar)}

#plot(seq(1, 10), disG(1:10), type = "h", col = "red")
#points(seq(1, 10), disG(1:10), type = "o", col = "red")
#points(seq(1, 10), disIG(1:10), type = "o", col = "blue")

plot(seq(1, 10), disWei(1:10, alpha = 0.5), pch = 19, col = "red", ylim = c(0, max(disWei(1:10, alpha = 0.5), disIWei(1:10, alpha = 0.5))), xlab = "n", ylab = "P(n)")
points(seq(1, 10), disWei(1:10, alpha = 0.5), col = "black", type = "h", lty = 3)
points(seq(1, 10), disIWei(1:10, alpha = 0.5), col = "blue", pch = 19)
points(seq(1, 10), disIWei(1:10, alpha = 0.5), col = "black", type = "h", lty = 3)
axis(1, at = seq(1:10))
#abline(v = seq(0, 10), col = "black", lty = 3)
legend(4, 0.4, pch = c(19, 19), col = c("red", "blue"),
       legend = c(expression(paste("Discrete Weibull(", alpha, " = 0.5)")), expression(paste("Inverse discrete Weibull(", alpha, " = 0.5)" ))), bty = "n", cex = 1.2)
#       legend = c("Geometric", "Inverse Geometric"), bty ="n", cex = 1.2)

plot(seq(1, 15), disWeicum(1:15, alpha = 0.5), ylim = c(0, max(disWeicum(1:15, alpha = 0.5), disIWeicum(1:15, alpha = 0.5))), type = "s", col = "red", xlab = "n", ylab = "F(n)")
#points(seq(1, 15), disGcum(1:15), type = "s", col = "red")
points(seq(1, 15), disIWeicum(1:15, alpha = 0.5), type = "s", col = "blue")
axis(1, at = seq(1:15))
abline(v = seq(0, 10), col = "grey", lty = 3)
#abline(h = c(0, 0.2, 0.4, 0.6, 0.8, 1), col = "grey", lty = 3)
legend(5, 0.5, col = c("red", "blue"), legend = c(expression(paste("Discrete Weibull(", alpha, " = 0.5)")), expression(paste("Inverse discrete Weibull(", alpha, " = 0.5)"))), lty = c(1, 1), bty = "n", cex = 1.2, seg.len = 0.3)

## Half-logistic
disHL <- function(k){
  ar <- array(NA, dim = length(k))
  for(i in k){
    ar[match(i, k)] <- ((1 - exp(-i))/(1 + exp(-i))) - ((1 - exp(-(i - 1)))/(1 + exp(-(i - 1))))
  }
  return(ar)}

disIHL <- function(k){
  jj <- array(NA, dim = length(k))
  for(i in k){
    pos <- match(i, k)
    jj[pos] <-  2*( ( exp(-1/i)/(1 + exp(-1/i)) ) - ( exp(-1/(i - 1))/(1 + exp(-1/(i - 1)))))
  }
  return(jj)}

disHLcum <- function(x, alpha){
  ar <- array(NA, dim = length(x))
  for(i in x){
    ar[match(i, x)] <- suppressWarnings(sum(disHL(k = 1:i)))
  }
  return(ar)}

disIHLcum <- function(x, alpha){
  ar <- array(NA, dim = length(x))
  for(i in x){
    ar[match(i, x)] <-suppressWarnings(sum(disIHL(1:i)))
  }
  return(ar)}

#plot(seq(1, 10), disG(1:10), type = "h", col = "red")
#points(seq(1, 10), disG(1:10), type = "o", col = "red")
#points(seq(1, 10), disIG(1:10), type = "o", col = "blue")

plot(seq(1, 10), disHL(1:10), pch = 19, col = "red", ylim = c(0, max(disHL(1:10), disIHL(1:10))), xlab = "n", ylab = "P(n)")
points(seq(1, 10), disHL(1:10), col = "black", type = "h", lty = 3)
points(seq(1, 10), disIHL(1:10), col = "blue", pch = 19)
points(seq(1, 10), disIHL(1:10), col = "black", type = "h", lty = 3)
axis(1, at = seq(1:10))
#abline(v = seq(0, 10), col = "black", lty = 3)
legend(4, 0.4, pch = c(19, 19), col = c("red", "blue"),
       legend = c("Discrete Half-Logistic", "Discrete inverse Half-Logistic"), bty = "n", cex = 1.2)
#       legend = c("Geometric", "Inverse Geometric"), bty ="n", cex = 1.2)

plot(seq(1, 15), disHLcum(1:15), ylim = c(0, max(disHLcum(1:15), disIHLcum(1:15))), type = "s", col = "red", xlab = "n", ylab = "F(n)")
#points(seq(1, 15), disGcum(1:15), type = "s", col = "red")
points(seq(1, 15), disIHLcum(1:15), type = "s", col = "blue")
axis(1, at = seq(1:15))
abline(v = seq(0, 10), col = "grey", lty = 3)
#abline(h = c(0, 0.2, 0.4, 0.6, 0.8, 1), col = "grey", lty = 3)
legend(5, 0.4, col = c("red", "blue"), legend = c("Discrete Half-logistic", "Discrete inverse Half-logistic"), lty = c(1, 1), bty = "n", cex = 1.2, seg.len = 0.3)

# Discrete Pareto
disPare <- function(k, alpha){
  ar <- array(NA, dim = length(k))
  for(i in k){
    ar[match(i, k)] <- ( 1/(1 + alpha*(i -1)))^(1/alpha) - ( 1/(1 + alpha*i))^(1/alpha)
  }
  return(ar)}

disIPare <- function(k, alpha){
  jj <- array(NA, dim = length(k))
  for(i in k){
    pos <- match(i, k)
    jj[pos] <-  (1/(1 + (alpha/i)))^(1/alpha) - (1/(1 + alpha/(i - 1)))^(1/alpha)
  }
  return(jj)}

disParecum <- function(x, alpha){
  ar <- array(NA, dim = length(x))
  for(i in x){
    ar[match(i, x)] <- suppressWarnings(sum(disPare(k = 1:i, alpha = alpha)))
  }
  return(ar)}

disIParecum <- function(x, alpha){
  ar <- array(NA, dim = length(x))
  for(i in x){
    ar[match(i, x)] <-suppressWarnings(sum(disIPare(1:i, alpha = alpha)))
  }
  return(ar)}

#plot(seq(1, 10), disG(1:10), type = "h", col = "red")
#points(seq(1, 10), disG(1:10), type = "o", col = "red")
#points(seq(1, 10), disIG(1:10), type = "o", col = "blue")

plot(seq(1, 10), disPare(1:10, alpha = 0.5), pch = 19, col = "red", ylim = c(0, max(disPare(1:10, alpha = 0.5), disIPare(1:10, alpha = 0.5))), xlab = "n", ylab = "P(n)")
points(seq(1, 10), disPare(1:10, alpha = 0.5), col = "black", type = "h", lty = 3)
points(seq(1, 10), disIPare(1:10, alpha = 0.5), col = "blue", pch = 19)
points(seq(1, 10), disIPare(1:10, alpha = 0.5), col = "black", type = "h", lty = 3)
axis(1, at = seq(1:10))
#abline(v = seq(0, 10), col = "black", lty = 3)
legend(3, 0.5, pch = c(19, 19), col = c("red", "blue"),
       legend = c(expression(paste("Discrete Pareto (", alpha, " = 0.5)")), expression(paste("Discrete Inverse Pareto (", alpha, " = 0.5)"))), bty = "n", cex = 1.2)
#       legend = c("Geometric", "Inverse Geometric"), bty ="n", cex = 1.2)

plot(seq(1, 15), disParecum(1:15, alpha = 0.5), ylim = c(0, max(disParecum(1:15, alpha = 0.5), disIParecum(1:15, alpha = 0.5))), type = "s", col = "red", xlab = "n", ylab = "F(n)")
#points(seq(1, 15), disGcum(1:15), type = "s", col = "red")
points(seq(1, 15), disIParecum(1:15, alpha = 0.5), type = "s", col = "blue")
axis(1, at = seq(1:15))
abline(v = seq(0, 10), col = "grey", lty = 3)
#abline(h = c(0, 0.2, 0.4, 0.6, 0.8, 1), col = "grey", lty = 3)
legend(5, 0.4, col = c("red", "blue"), legend = c(expression(paste("Discrete Pareto (", alpha, " = 0.5)")), expression(paste("Discrete Inverse Pareto (", alpha, " = 0.5)"))), lty = c(1, 1), bty = "n", cex = 1.2, seg.len = 0.3)

### Poisson
disIPois <- function(k, lambda){
  jj <- array(NA, dim = length(k))
  for(i in k){
    pos <- match(i, k)
     jj[pos] <- pgamma(lambda, shape = 1/i, lower.tail = FALSE) - pgamma(lambda, shape = 1/(i + 1), lower.tail = FALSE)
  }
  return(jj)}

disIPoiscum <- function(x, lambda){
  ar <- array(NA, dim = length(x))
  for(i in x){
    ar[match(i, x)] <-suppressWarnings(sum(disIPois(0:i, lambda = lambda)))
  }
  return(ar)}

#plot(seq(1, 10), disG(1:10), type = "h", col = "red")
#points(seq(1, 10), disG(1:10), type = "o", col = "red")
#points(seq(1, 10), disIG(1:10), type = "o", col = "blue")

plot(seq(0, 10), dpois(0:10, lambda = 2), pch = 19, col = "red", ylim = c(0, max(dpois(0:10, lambda = 2), disIPois(0:10, lambda = 2))), xlab = "n", ylab = "P(n)")
points(seq(0, 10), dpois(0:10, lambda = 2), col = "black", type = "h", lty = 3)
points(seq(0, 10), disIPois(0:10, lambda = 2), col = "blue", pch = 19)
points(seq(0, 10), disIPois(0:10, lambda = 2), col = "black", type = "h", lty = 3)
axis(1, at = seq(1:10))
#abline(v = seq(0, 10), col = "black", lty = 3)
legend(2, 0.6, pch = c(19, 19), col = c("red", "blue"),
       legend = c(expression(paste("Poisson (", lambda, " = 2)")), expression(paste("Discrete Inverse Poisson (", lambda, " = 2)"))), bty = "n", cex = 1.2)
#       legend = c("Geometric", "Inverse Geometric"), bty ="n", cex = 1.2)

plot(seq(0, 15), ppois(0:15, lambda = 2), ylim = c(0, max(ppois(0:15, lambda = 2), disIPoiscum(0:15, lambda = 2))), type = "s", col = "red", xlab = "n", ylab = "F(n)")
#points(seq(1, 15), disGcum(1:15), type = "s", col = "red")
points(seq(0, 15), disIPoiscum(0:15, lambda = 2), type = "s", col = "blue")
axis(1, at = seq(1:15))
abline(v = seq(0, 10), col = "grey", lty = 3)
#abline(h = c(0, 0.2, 0.4, 0.6, 0.8, 1), col = "grey", lty = 3)
legend(5, 0.4, col = c("red", "blue"), legend = c(expression(paste("Poisson (", lambda, " = 2)")), expression(paste("Discrete Inverse Poisson (", lambda, " = 2)"))), lty = c(1, 1), bty = "n", cex = 1.2, seg.len = 0.4)

# Negative binomial distribution and its inverse
disINB <- function(k, r, beta){
  jj <- array(NA, dim = length(k))
  for(i in k){
    if(i == 0){
      jj[match(i, k)] <- (beta^r/(gamma(r)))*(integrate(f = function(x) pgamma(x, shape = 1/(i + 1))*x^(r - 1)*exp(-beta*x), lower = 0, upper = Inf)$value)
    }
    else{
      pos <- match(i, k)
      jj[pos] <- (beta^r/(gamma(r)))*(integrate(f = function(x) pgamma(x, shape = 1/(i + 1))*x^(r - 1)*exp(-beta*x), lower = 0, upper = Inf)$value -
        integrate(f = function(x) pgamma(x, shape = 1/i)*x^(r - 1)*exp(-beta*x), lower = 0, upper = Inf)$value )
    }
  }
  return(jj)}

disINBcum <- function(x, r, beta){
  ar <- array(NA, dim = length(x))
  for(i in x){
    ar[match(i, x)] <-suppressWarnings(sum(disINB(0:i, r = r, beta = beta)))
  }
  return(ar)}

#plot(seq(1, 10), disG(1:10), type = "h", col = "red")
#points(seq(1, 10), disG(1:10), type = "o", col = "red")
#points(seq(1, 10), disIG(1:10), type = "o", col = "blue")

plot(seq(0, 10), dnbinom(0:10, size = 3, prob = 2/3), pch = 19, col = "red", ylim = c(0, max(dnbinom(0:10, size = 3, prob = 2/3), disINB(0:10, r = 3, beta = 2))), xlab = "n", ylab = "P(n)")
points(seq(0, 10), dnbinom(0:10, size = 3, prob = 2/3), col = "black", type = "h", lty = 3)
points(seq(0, 10), disINB(0:10, r = 3, beta = 2), col = "blue", pch = 19)
points(seq(0, 10), disINB(0:10, r=3, beta = 2), col = "black", type = "h", lty = 3)
axis(1, at = seq(1:10))
#abline(v = seq(0, 10), col = "black", lty = 3)
legend(1, 0.5, pch = c(19, 19), col = c("red", "blue"),
       legend = c(expression(paste("Negative Binomial (", r, " = 3, ", beta, " = 2)")), expression(paste("Inverse Negative Binomial (", r, " = 3, ", beta, " = 2)"))), bty = "n", cex = 1.2)
#       legend = c("Geometric", "Inverse Geometric"), bty ="n", cex = 1.2)

plot(seq(0, 15), pnbinom(0:15, size = 3, prob = 2/3), ylim = c(0, max(pnbinom(0:15, size = 3, prob = 2/3), disINBcum(0:15, r = 3, beta = 2))), type = "s", col = "red", xlab = "n", ylab = "F(n)")
#points(seq(1, 15), disGcum(1:15), type = "s", col = "red")
points(seq(0, 15), disINBcum(0:15, r = 3, beta = 2), type = "s", col = "blue")
axis(1, at = seq(1:15))
abline(v = seq(0, 10), col = "grey", lty = 3)
#abline(h = c(0, 0.2, 0.4, 0.6, 0.8, 1), col = "grey", lty = 3)
legend(4, 0.4, col = c("red", "blue"), legend = c(expression(paste("Negative Binomial (", r, " = 3, ", beta, " = 2)")), expression(paste("Inverse Negative Binomial (", r, " = 3, ", beta, " = 2)"))), lty = c(1, 1), bty = "n", cex = 1.2, seg.len = 0.4)

# Generalized Sibuya distribution
disGS <- function(k, nu, alpha){
  ar <- array(NA, dim = length(k))
  for(i in k){
    ar[match(i, k)] <- (gamma(nu + 1)*gamma(nu + i - alpha))/(gamma(nu + 1 - alpha)*gamma(nu + i)) -
                        (gamma(nu + 1)*gamma(nu + 1+i - alpha))/(gamma(nu + 1 - alpha)*gamma(nu + 1 + i))
  }
  return(ar)}

disIGS <- function(k, nu, alpha){
  jj <- array(NA, dim = length(k))
  for(i in k){
    pos <- match(i, k)
    if(i == 1){
     jj[pos] <- (gamma(nu + 1)*gamma(nu + 1 + 1/i - alpha))/(gamma(nu + 1 - alpha)*gamma(nu + 1 + 1/i))
    }
    else{
    jj[pos] <-  -( (gamma(nu + 1)*gamma(nu + 1 + 1/(i-1) - alpha ))/(gamma(nu + 1 - alpha)*gamma(nu + 1 + 1/(i-1))) -
      (gamma(nu + 1)*gamma(nu + 1 + 1/i - alpha))/(gamma(nu + 1 - alpha)*gamma(nu + 1 + 1/i)) )
    }
  }
  return(jj)}

disGScum <- function(x, nu, alpha){
  ar <- array(NA, dim = length(x))
  for(i in x){
    ar[match(i, x)] <- suppressWarnings(sum(disGS(k = 1:i, nu = nu, alpha = alpha)))
  }
  return(ar)}

disIGScum <- function(x, nu, alpha){
  ar <- array(NA, dim = length(x))
  for(i in x){
    ar[match(i, x)] <-suppressWarnings(sum(disIGS(1:i, nu = nu, alpha = alpha)))
  }
  return(ar)}

plot(seq(1, 10), disGS(1:10, nu = 2, alpha = 1.5), pch = 19, col = "red", ylim = c(0, max(disGS(1:10, nu = 2, alpha = 1.5),
                                                                                          disIGS(1:10, nu = 2, alpha = 1.5))), xlab = "n", ylab = "P(n)")
points(seq(1, 10), disGS(1:10, nu = 2, alpha = 1.5), col = "black", type = "h", lty = 3)
points(seq(1, 10), disIGS(1:10, nu = 2, alpha = 1.5), col = "blue", pch = 19)
points(seq(1, 10), disIGS(1:10, nu = 2, alpha = 1.5), col = "black", type = "h", lty = 3)
axis(1, at = seq(1:10))
#abline(v = seq(0, 10), col = "black", lty = 3)
legend(2.5, 0.4, pch = c(19, 19), col = c("red", "blue"),
       legend = c(expression(paste("Generalized Sibuya(", nu, " = 2, ", alpha, " = 1.5)" )),
                  expression(paste("Generalized Inverse Sibuya(", nu, " = 2, ", alpha, " = 1.5)" ))), bty = "n", cex = 1.2)

plot(seq(1, 15), disGScum(1:15, nu = 2, alpha = 1.5), ylim = c(0, max(disGScum(1:15, alpha = 1.5, nu = 2), disIGScum(1:15, alpha = 1.5, nu = 2))), type = "s", col = "red", xlab = "n", ylab = "F(n)")
#points(seq(1, 15), disGcum(1:15), type = "s", col = "red")
points(seq(1, 15), disIGScum(1:15, alpha = 1.5, nu =2), type = "s", col = "blue")
axis(1, at = seq(1:15))
abline(v = seq(0, 15), col = "grey", lty = 3)
#abline(h = c(0, 0.2, 0.4, 0.6, 0.8, 1), col = "grey", lty = 3)
legend(3, 0.4, col = c("red", "blue"), legend = c(expression(paste("Generalized Sibuya(", nu, " = 2, ", alpha, " = 1.5)" )),
                                                  expression(paste("Generalized Inverse Sibuya(", nu, " = 2, ", alpha, " = 1.5)" ))), lty = c(1, 1), bty = "n", cex = 1.2, seg.len = 0.3)

# Zeta- function
library(VGAM)

# logarithmic distribution
dislog <- function(k, q){
  ar <- array(NA, dim = length(k))
  for(i in k){
    ar[match(i, k)] <- - (q^i)/(log(1 - q)*i)
  }
  return(ar)}

disIlog <- function(k, q){
  jj <- array(NA, dim = length(k))
  for(i in k){
    pos <- match(i, k)
    jj[pos] <- (-1/log(1 - q))*integrate(function(x){ ( x^(1/i) - x^(1/(i - 1)))/(1 - x) }, lower = 0, upper = q)$value
  }
  return(jj)}

dislogcum <- function(x, q){
  ar <- array(NA, dim = length(x))
  for(i in x){
    ar[match(i, x)] <- suppressWarnings(sum(dislog(k = 1:i, q = q)))
  }
  return(ar)}

disIlogcum <- function(x, q){
  ar <- array(NA, dim = length(x))
  for(i in x){
    ar[match(i, x)] <-suppressWarnings(sum(disIlog(1:i, q = q)))
  }
  return(ar)}

plot(seq(1, 10), dislog(1:10, q = 0.7), pch = 19, col = "red", ylim = c(0, max(dislog(1:10, q = 0.7),
                                                                                          dislog(1:10, q = 0.7))), xlab = "n", ylab = "P(n)")
points(seq(1, 10), dislog(1:10, q = 0.7), col = "black", type = "h", lty = 3)
points(seq(1, 10), disIlog(1:10, q = 0.7), col = "blue", pch = 19)
points(seq(1, 10), disIlog(1:10, q = 0.7), col = "black", type = "h", lty = 3)
axis(1, at = seq(1:10))
#abline(v = seq(0, 10), col = "black", lty = 3)
legend(4, 0.4, pch = c(19, 19), col = c("red", "blue"), legend = c(expression(paste("logarithmic (", q," = 0.7)")),
                  expression(paste("Inverse logarithmic(", q, " = 0.7)" ))), bty = "n", cex = 1.2)

plot(seq(1, 15), dislogcum(1:15, q = 0.7), ylim = c(0, max(dislogcum(1:15, q = 0.7), disIlogcum(1:15, q = 0.7))), type = "s", col = "red", xlab = "n", ylab = "F(n)")
#points(seq(1, 15), disGcum(1:15), type = "s", col = "red")
points(seq(1, 15), disIlogcum(1:15, q = 0.7), type = "s", col = "blue")
axis(1, at = seq(1:15))
abline(v = seq(0, 15), col = "grey", lty = 3)
#abline(h = c(0, 0.2, 0.4, 0.6, 0.8, 1), col = "grey", lty = 3)
legend(3, 0.3, col = c("red", "blue"), legend = c(expression(paste("logarithmic(", q, " = 0.7)" )),
                                                  expression(paste("Inverse logarithmic(", q, " = 0.7)" ))), lty = c(1, 1), bty = "n", cex = 1.2, seg.len = 0.5)

# zipf distribution
diszipf <- function(k, alpha){
  ar <- array(NA, dim = length(k))
  for(i in k){
    ar[match(i, k)] <- (1/zeta(x = alpha + 1))*((1/(i + 1))^(alpha + 1) )
  }
  return(ar)}

disIzipf <- function(k, alpha){
  jj <- array(NA, dim = length(k))
  for(i in k){
    pos <- match(i, k)
    jj[pos] <- (1/zeta(x = alpha + 1))*( zeta(x = alpha + 1, shift = ((i + 2)/(i + 1))) - zeta(x = alpha + 1, shift = ((i + 1)/(i))))
  }
  return(jj)}

diszipfcum <- function(x, alpha){
  ar <- array(NA, dim = length(x))
  for(i in x){
    ar[match(i, x)] <- suppressWarnings(sum(diszipf(k = 0:i, alpha = alpha)))
  }
  return(ar)}

disIzipfcum <- function(x, alpha){
  ar <- array(NA, dim = length(x))
  for(i in x){
    ar[match(i, x)] <-suppressWarnings(sum(disIzipf(0:i, alpha = alpha)))
  }
  return(ar)}

plot(seq(0, 10), diszipf(0:10, alpha = 2), pch = 19, col = "red", ylim = c(0, max(diszipf(0:10, alpha = 2),
                                                                               diszipf(0:10, alpha = 2))), xlab = "n", ylab = "P(n)")
points(seq(0, 10), diszipf(0:10, alpha = 2), col = "black", type = "h", lty = 3)
points(seq(0, 10), disIzipf(0:10, alpha = 2), col = "blue", pch = 19)
points(seq(0, 10), disIzipf(0:10, alpha = 2), col = "black", type = "h", lty = 3)
axis(1, at = seq(1:10))
#abline(v = seq(0, 10), col = "black", lty = 3)
legend(4, 0.4, pch = c(19, 19), col = c("red", "blue"), legend = c(expression(paste("zipf (", alpha," = 2)")),
                                                                   expression(paste("Inverse zipf(", alpha, " = 2)" ))), bty = "n", cex = 1.2)

plot(seq(0, 15), diszipfcum(0:15, alpha = 2), ylim = c(0, max(diszipfcum(0:15, alpha = 2), disIzipfcum(0:15, alpha = 2))), type = "s", col = "red", xlab = "n", ylab = "F(n)")
#points(seq(1, 15), disGcum(1:15), type = "s", col = "red")
points(seq(0, 15), disIzipfcum(0:15, alpha = 2), type = "s", col = "blue")
axis(1, at = seq(0:15))
abline(v = seq(0, 15), col = "grey", lty = 3)
#abline(h = c(0, 0.2, 0.4, 0.6, 0.8, 1), col = "grey", lty = 3)
legend(3, 0.3, col = c("red", "blue"), legend = c(expression(paste("zipf(", alpha, " = 2)" )),
                                                  expression(paste("Inverse zipf(", alpha, " = 2)" ))), lty = c(1, 1), bty = "n", cex = 1.2, seg.len = 0.5)
