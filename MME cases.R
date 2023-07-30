0########## Normal Pareto Generator ###################
rNPare <- function(n = 1, alpha = 0.5, sigma = 1){
  normal <- rnorm(n); exponential <- rexp(n)
  gamm <- rgamma(n, shape = 1/alpha, scale = alpha)
  y <- sigma*sqrt(exponential/gamm)*normal
  return(y)
}

########################## MME case 1 sigma unknown but alpha known ##########################
MME.case1.sigma <- function(alpha, data){
  fin <- (gamma(3/4)*gamma(1/4)*gamma(1/alpha + 0.25)*(alpha^0.25)/(2^0.25 * sqrt(pi)*gamma(1/alpha)*mean( abs(data)^(-0.5))))^2
  return(fin)
}

MME.case1.sigma(data = rNPare(5000), alpha = 0.5) # An example

sim.MME.case1.sigma <- function(k = 1000, alpha = 0.5, sigma = 1, num = 10){
  num.vec <- matrix(NA, nrow = k, ncol = 5)
  colnames(num.vec) <- c("Sample", "alpha", "sigma", "Estimate", "RelativeError")
  num.vec[ , 1] <- rep(num, k)
  num.vec[ , 2] <- rep(alpha, k)
  num.vec[ , 3] <- rep(sigma, k)

  # Generate independent Y's
  for(i in 1:k){
    cc <- rNPare(n = num, alpha = alpha, sigma = sigma)
    num.vec[i, 4] <- MME.case1.sigma(alpha = alpha, data = cc)
    num.vec[i, 5] <- (num.vec[i, 4] - num.vec[i, 3])/num.vec[i, 3]
  }
  #  return(data.frame(num.vec))
  return(data.frame(num.vec))
}

l <- rbind(sim.MME.case1.sigma(alpha = 0.7, sigma = 0.5, num = 50), sim.MME.case1.sigma(alpha = 0.7, sigma = 0.5, num = 150), sim.MME.case1.sigma(alpha = 0.7, sigma = 0.5, num = 300),
           sim.MME.case1.sigma(alpha = 1, sigma = 1, num = 50), sim.MME.case1.sigma(alpha = 1, sigma = 1, num = 150), sim.MME.case1.sigma(alpha = 1, sigma = 1, num = 300),
           sim.MME.case1.sigma(alpha = 2, sigma = 5, num = 50), sim.MME.case1.sigma(alpha = 2, sigma = 5, num = 150), sim.MME.case1.sigma(alpha = 2, sigma = 5, num = 300),
           sim.MME.case1.sigma(alpha = 1, sigma = 0.5, num = 50), sim.MME.case1.sigma(alpha = 2, sigma = 0.5, num = 50), sim.MME.case1.sigma(alpha = 1, sigma = 0.5, num = 150),
           sim.MME.case1.sigma(alpha = 2, sigma = 0.5, num = 150), sim.MME.case1.sigma(alpha = 1, sigma = 0.5, num = 300), sim.MME.case1.sigma(alpha = 2, sigma = 0.5, num = 300),
           sim.MME.case1.sigma(alpha = 0.7, sigma = 1, num = 50), sim.MME.case1.sigma(alpha = 2, sigma = 1, num = 50), sim.MME.case1.sigma(alpha = 0.7, sigma = 1, num = 150),
           sim.MME.case1.sigma(alpha = 2, sigma = 1, num = 150), sim.MME.case1.sigma(alpha = 0.7, sigma = 1, num = 300), sim.MME.case1.sigma(alpha = 2, sigma = 1, num = 300),
           sim.MME.case1.sigma(alpha = 0.7, sigma = 5, num = 50), sim.MME.case1.sigma(alpha = 1, sigma = 5, num = 50), sim.MME.case1.sigma(alpha = 0.7, sigma = 5, num = 150),
           sim.MME.case1.sigma(alpha = 1, sigma = 5, num = 150), sim.MME.case1.sigma(alpha = 0.7, sigma = 5, num = 300), sim.MME.case1.sigma(alpha = 1, sigma = 5, num = 300))

save(l, file = "MME.case1.sigma.1000.simulations.RData")

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
tab1 <- l[which(l$Sample == 50 & l$sigma == 0.5 & l$alpha == 0.7), ]
(esti <- round(mean(tab1$Estimate), 4))
(MSE <- round(mean( (tab1$sigma - tab1$Estimate)^2 ), 4 ))

tab1 <- l[which(l$Sample == 300 & l$sigma == 5 & l$alpha == 2), ]
(esti <- round(mean(tab1$Estimate), 4 ))
(MSE <- round(mean( (tab1$sigma - tab1$Estimate)^2 ), 4 ))

f <- function(alpha){
  h <- log(alpha) + digamma(1/alpha)
  return(h)
}

curve(f, from = 0.01, to = 70)

g <- function(z){
  h <- (gamma(z + 0.25)/gamma(z))/(z^0.25)
  return(h)
}

curve(g, from = 0.001, to = 10000, xlab = expression(paste(alpha)), ylab = expression(paste("G(",alpha,")")))

##################### MME case 2 sigma unknown ########################
library(rootSolve)
MME.case2.alpha <- function(sigma, data){
  fin <- function(alpha){
    h <- log(alpha) + digamma(1/alpha) - 2*log(sigma) - digamma(1) - log(2) - digamma(1/2) + 2*mean(log(abs(data)))
    return(h)
  }
  return(uniroot(fin, lower = 0.01, upper = 100, extendInt = "yes")$root)
}

MME.case2.alpha(data = rNPare(500, alpha = 2), sigma = 1) # An example

sim.MME.case2.alpha <- function(k = 8000, alpha = 0.5, sigma = 1, num = 10){
  num.vec <- matrix(NA, nrow = k, ncol = 5)
  colnames(num.vec) <- c("Sample", "alpha", "sigma", "Estimate", "RelativeError")
  num.vec[ , 1] <- rep(num, k)
  num.vec[ , 2] <- rep(alpha, k)
  num.vec[ , 3] <- rep(sigma, k)

  # Generate independent Y's
  for(i in 1:k){
    cc <- rNPare(n = num, alpha = alpha, sigma = sigma)
    y <- suppressWarnings(try(MME.case2.alpha(alpha = alpha, data = cc)))
    if(y == "try-error"){
      num.vec[i, 4] <- NA
    }
    else{
      q <- suppressWarnings(try(x <- MME.case2.alpha(sigma = sigma, data = cc)))
        if(class(q) != "try-error"){
          num.vec[i, 4] <- q
        }
      }
    num.vec[i, 5] <- (num.vec[i, 4] - num.vec[i, 2])/num.vec[i, 2]
  }
  v <- data.frame(num.vec)
  v <- v[complete.cases(v), ]
  v <- v[1:1000, ]
  return(v)
}

l <- rbind(sim.MME.case2.alpha(alpha = 0.7, sigma = 0.5, num = 50), sim.MME.case2.alpha(alpha = 0.7, sigma = 0.5, num = 150), sim.MME.case2.alpha(alpha = 0.7, sigma = 0.5, num = 300),
           sim.MME.case2.alpha(alpha = 1, sigma = 1, num = 50), sim.MME.case2.alpha(alpha = 1, sigma = 1, num = 150), sim.MME.case2.alpha(alpha = 1, sigma = 1, num = 300),
           sim.MME.case2.alpha(alpha = 2, sigma = 5, num = 50), sim.MME.case2.alpha(alpha = 2, sigma = 5, num = 150), sim.MME.case2.alpha(alpha = 2, sigma = 5, num = 300),
           sim.MME.case2.alpha(alpha = 1, sigma = 0.5, num = 50), sim.MME.case2.alpha(alpha = 2, sigma = 0.5, num = 50), sim.MME.case2.alpha(alpha = 1, sigma = 0.5, num = 150),
           sim.MME.case2.alpha(alpha = 2, sigma = 0.5, num = 150), sim.MME.case2.alpha(alpha = 1, sigma = 0.5, num = 300), sim.MME.case2.alpha(alpha = 2, sigma = 0.5, num = 300),
           sim.MME.case2.alpha(alpha = 0.7, sigma = 1, num = 50), sim.MME.case2.alpha(alpha = 2, sigma = 1, num = 50), sim.MME.case2.alpha(alpha = 0.7, sigma = 1, num = 150),
           sim.MME.case2.alpha(alpha = 2, sigma = 1, num = 150), sim.MME.case2.alpha(alpha = 0.7, sigma = 1, num = 300), sim.MME.case2.alpha(alpha = 2, sigma = 1, num = 300),
           sim.MME.case2.alpha(alpha = 0.7, sigma = 5, num = 50), sim.MME.case2.alpha(alpha = 1, sigma = 5, num = 50), sim.MME.case2.alpha(alpha = 0.7, sigma = 5, num = 150),
           sim.MME.case2.alpha(alpha = 1, sigma = 5, num = 150), sim.MME.case2.alpha(alpha = 0.7, sigma = 5, num = 300), sim.MME.case2.alpha(alpha = 1, sigma = 5, num = 300))

save(l, file = "MME.case2.alpha.unknown.RData")

library(ggplot2)
# Relative error graphs
gh <- ggplot(l, aes(x = RelativeError)) + geom_histogram(aes(y = ..density..)) + facet_grid( alpha ~ sigma + Sample, scales = "free", labeller = label_both) + theme_classic()
gh <- gh + labs(x = "Relative Error")
gh


# Boxplots
th <- ggplot(l, aes(x = as.factor(Sample), y = Estimate)) + geom_boxplot() + geom_hline(aes(yintercept = alpha), color = "red") +
  facet_wrap(alpha ~ sigma, scales = "free", labeller = label_both) + labs(x = "Sample size") + theme_classic()
th

# Generating table
tab1 <- l[which(l$Sample == 50 & l$sigma == 0.5 & l$alpha == 0.7), ]
(esti <- round(mean(tab1$Estimate), 4 ))
(MSE <- round(mean( (tab1$alpha - tab1$Estimate)^2 ), 4 ))

tab1 <- l[which(l$Sample == 300 & l$sigma == 5 & l$alpha == 2), ]
(esti <- round(mean(tab1$Estimate), 4 ))
(MSE <- round(mean( (tab1$alpha - tab1$Estimate)^2 ), 4 ))


# Constant is negative addressed proportion of zeros
delta <- function(sigma, alpha, size = 10){
  ikr <- log(2*(sigma^2)) + digamma(1) + digamma(1/2) - 2*mean(log(abs(rNPare(n=size, alpha = alpha, sigma = sigma))))
  return(ikr)
}

ty.delta <- function(alpha, sigma, n = 100, k = 10000){
  num.vec <- matrix(NA, nrow = k, ncol = 4)
  colnames(num.vec) <- c("sigma", "alpha", "sample size", "delta < 0")
  num.vec[ , 1] <- rep(sigma, k)
  num.vec[ , 2] <- rep(alpha, k)
  num.vec[ , 3] <- rep(n, k)

  for(i in 1:k){
    if(delta(alpha = alpha, sigma = sigma, size = n) < 0){num.vec[i, 4] <- "Yes" }
    else{num.vec[i, 4] <- "No"}
  }
  data.frame(sigma = num.vec[ , 1], alpha = num.vec[ , 2], sampleSize = num.vec[ , 3], Negative_delta_Value = num.vec[ ,4] )
}

l <- rbind(ty.delta(alpha = 0.7, sigma = 0.5, n = 50), ty.delta(alpha = 0.7, sigma = 0.5, n = 150), ty.delta(alpha = 0.7, sigma = 0.5, n = 300),
           ty.delta(alpha = 0.7, sigma = 1, n = 50), ty.delta(alpha = 0.7, sigma = 1, n = 150), ty.delta(alpha = 0.7, sigma = 1, n = 300),
           ty.delta(alpha = 0.7, sigma = 5, n = 50), ty.delta(alpha = 0.7, sigma = 5, n = 150), ty.delta(alpha = 0.7, sigma = 5, n = 300),
           ty.delta(alpha = 1, sigma = 0.5, n = 50), ty.delta(alpha = 1, sigma = 0.5, n = 150), ty.delta(alpha = 1, sigma = 0.5, n = 300),
           ty.delta(alpha = 1, sigma = 1, n = 50), ty.delta(alpha = 1, sigma = 1, n = 150), ty.delta(alpha = 1, sigma = 1, n = 300),
           ty.delta(alpha = 1, sigma = 5, n = 50), ty.delta(alpha = 1, sigma = 5, n = 150), ty.delta(alpha = 1, sigma = 5, n = 300),
           ty.delta(alpha = 2, sigma = 0.5, n = 50), ty.delta(alpha = 2, sigma = 0.5, n = 150), ty.delta(alpha = 2, sigma = 0.5, n = 300),
           ty.delta(alpha = 2, sigma = 1, n = 50), ty.delta(alpha = 2, sigma = 1, n = 150), ty.delta(alpha = 2, sigma = 1, n = 300),
           ty.delta(alpha = 2, sigma = 5, n = 50), ty.delta(alpha = 2, sigma = 5, n = 150), ty.delta(alpha = 2, sigma = 5, n = 300))

library(ggplot2)

#level_order <- c('50', '150', '300')
#lvl1 <- c('0.7', '1', '2')
#lvl2 <- c('0.5', '1', '5')
gh <- ggplot(l, aes(fill = Negative_delta_Value, x = factor(sampleSize, levels = level_order))) +
  geom_bar( ) + facet_wrap(sigma~ alpha, labeller = label_both ) + xlab("Sample size") + ylab("Number of simulations") +
  guides(fill=guide_legend(title=expression(paste(delta, " < 0")))) + theme_classic()

save(l, file = "MME_case2_delta_countings.RData")

#gh <- ggplot(l, aes(fill = Negative_delta_Value, x = factor(sampleSize, levels = level_order))) +
#  geom_bar( ) + facet_grid(~ factor(alpha, levels = lvl), labeller = label_value) + xlab("Sample size") + ylab("Number of simulations") +
#  guides(fill=guide_legend(title="Positive R value")) + theme_classic()


######## Case 3: MME where sigma and alpha are unknown #################
MME.case3.alpha.sigma <- function(data){
  R <- mean( log(abs(data)) ) - 0.5*(log(2) + digamma(1/2) + digamma(1)) - 2*log(gamma(3/4)*gamma(1/4)/(sqrt(pi)*(2^0.25)*mean(abs(data)^(-0.5))))

  fin <- function(alpha){
    h <- 0.5*(log(1/alpha) - digamma(1/alpha)) + 2*log((alpha^0.25)*gamma(1/alpha + 0.25)/gamma(1/alpha)) - R
    return(h)
  }
  alpha <- uniroot(fin, lower = 0.01, upper = 100, extendInt = "yes")$root
  sigma <- ( (alpha^0.25 * gamma(1/alpha + 0.25)*gamma(3/4)*gamma(1/4)  )/(mean(abs(data)^(-0.5))*gamma(1/2)*(2^0.25)*gamma(1/alpha)))^2
  return(c(alpha, sigma))
}

sim.MME.case3.alpha.sigma <- function(k = 8000, alpha = 0.5, sigma = 1, num = 10, dat.num = 1){
  num.vec <- matrix(NA, nrow = 2*k, ncol = 5)
  num.vec[ , 1] <- rep(num, 2*k)
  num.vec[ , 2] <- c(rep("alpha", k), rep("sigma", k))
  num.vec[ , 3] <- c(rep(alpha, k), rep(sigma, k))

  # Generate independent Y's
  for(i in 1:k){
    cc <- rNPare(n = num, alpha = alpha, sigma = sigma)
    y <- suppressWarnings(try(MME.case3.alpha.sigma(data = cc)))
    if(class(y) == "try-error"){
      num.vec[i, 4] <- NA
    }
    else{
      q <- suppressWarnings(try(x <- MME.case3.alpha.sigma(data = cc)))
      if(class(q) != "try-error"){
      #  parameter1[i] <- q[1]; parameter2[i] <- q[2]
      num.vec[i, 4] <- q[1]; num.vec[i+k, 4] <- q[2]
      }
    }
   num.vec[i, 5] <- (as.numeric(num.vec[i, 4]) - as.numeric(num.vec[i, 3]) )/ as.numeric(num.vec[i, 3])
   num.vec[i+k, 5] <- (as.numeric(num.vec[i+k, 4]) - as.numeric(num.vec[i+k, 3]))/as.numeric(num.vec[i+k, 3])
  }
  num.vec <- num.vec[complete.cases(num.vec), ]
  num.vec <- cbind(num.vec, paste("Simulation set", dat.num))
  colnames(num.vec) <- c("Sample", "Parameter", "True", "Estimate", "RelativeError", "SetNumber")
  num.vec <- num.vec[1:1000, ]
  return(data.frame(num.vec))
}

l <- rbind(sim.MME.case3.alpha.sigma(alpha = 0.7, sigma = 0.5, num = 50, dat.num = 1), sim.MME.case3.alpha.sigma(alpha = 0.7, sigma = 0.5, num = 150,dat.num = 2), sim.MME.case3.alpha.sigma(alpha = 0.7, sigma = 0.5, num = 300, dat.num = 3),
           sim.MME.case3.alpha.sigma(alpha = 1, sigma = 1, num = 50, dat.num = 4), sim.MME.case3.alpha.sigma(alpha = 1, sigma = 1, num = 150, dat.num = 5), sim.MME.case3.alpha.sigma(alpha = 1, sigma = 1, num = 300, dat.num = 6),
           sim.MME.case3.alpha.sigma(alpha = 2, sigma = 5, num = 50, dat.num = 7), sim.MME.case3.alpha.sigma(alpha = 2, sigma = 5, num = 150, dat.num = 8), sim.MME.case3.alpha.sigma(alpha = 2, sigma = 5, num = 300, dat.num = 9),
           sim.MME.case3.alpha.sigma(alpha = 1, sigma = 0.5, num = 50, dat.num = 10), sim.MME.case3.alpha.sigma(alpha = 2, sigma = 0.5, num = 50, dat.num = 11), sim.MME.case3.alpha.sigma(alpha = 1, sigma = 0.5, num = 150, dat.num = 12),
           sim.MME.case3.alpha.sigma(alpha = 2, sigma = 0.5, num = 150, dat.num = 13), sim.MME.case3.alpha.sigma(alpha = 1, sigma = 0.5, num = 300, dat.num = 14), sim.MME.case3.alpha.sigma(alpha = 2, sigma = 0.5, num = 300, dat.num = 15),
           sim.MME.case3.alpha.sigma(alpha = 0.7, sigma = 1, num = 50, dat.num = 16), sim.MME.case3.alpha.sigma(alpha = 2, sigma = 1, num = 50, dat.num = 17), sim.MME.case3.alpha.sigma(alpha = 0.7, sigma = 1, num = 150, dat.num = 18),
           sim.MME.case3.alpha.sigma(alpha = 2, sigma = 1, num = 150, dat.num = 19), sim.MME.case3.alpha.sigma(alpha = 0.7, sigma = 1, num = 300, dat.num = 20), sim.MME.case3.alpha.sigma(alpha = 2, sigma = 1, num = 300, dat.num = 21),
           sim.MME.case3.alpha.sigma(alpha = 0.7, sigma = 5, num = 50, dat.num = 22), sim.MME.case3.alpha.sigma(alpha = 1, sigma = 5, num = 50, dat.num = 23), sim.MME.case3.alpha.sigma(alpha = 0.7, sigma = 5, num = 150, dat.num =24),
           sim.MME.case3.alpha.sigma(alpha = 1, sigma = 5, num = 150, dat.num = 25), sim.MME.case3.alpha.sigma(alpha = 0.7, sigma = 5, num = 300, dat.num = 26), sim.MME.case3.alpha.sigma(alpha = 1, sigma = 5, num = 300, dat.num = 27))

save(l, file = "MME.case3.alpha.sigma.RData")

library(ggplot2)
# Relative error graphs
gh <- ggplot(l, aes(x = RelativeError)) + geom_histogram(aes(y = ..density..)) + facet_grid( alpha ~ sigma + Sample, scales = "free", labeller = label_both) + theme_classic()
gh <- gh + labs(x = "Relative Error")
gh

# Boxplots
th <- ggplot(l, aes(x = as.factor(Sample), y = Estimate)) + geom_boxplot() + geom_hline(aes(yintercept = alpha), color = "red") +
  facet_wrap(alpha ~ sigma, scales = "free", labeller = label_both) + labs(x = "Sample size") + theme_classic()
th

th <- ggplot(l, aes(x = as.factor(Sample), y = Estimate)) + geom_boxplot() + facet_wrap(al)



+ geom_boxplot() + geom_hline(aes(yintercept = alpha), color = "red") +
  facet_wrap(alpha ~ sigma, scales = "free", labeller = label_both) + labs(x = "Sample size") + theme_classic()
th


sim.MME.case3.alpha.sigma.unknown <- function(k = 8000, alpha = 0.5, sigma = 1, num = 10){
  num.vec <- matrix(NA, nrow = k, ncol = 7)
  num.vec[ , 1] <- rep(num, k)
  num.vec[ , 2] <- rep(alpha, k)
  num.vec[ , 3] <- rep(sigma, k)

  # Generate independent Y's
  for(i in 1:k){
    cc <- rNPare(n = num, alpha = alpha, sigma = sigma)
    y <- suppressWarnings(try(MME.case3.alpha.sigma(data = cc)))
    if(class(y) == "try-error"){
      num.vec[i, 4] <- NA
    }
    else{
      q <- suppressWarnings(try(x <- MME.case3.alpha.sigma(data = cc)))
      if(class(q) != "try-error"){
        num.vec[i, 4] <- q[1]; num.vec[i, 5] <- q[2]
      }
    }
    num.vec[i, 6] <- (as.numeric(num.vec[i, 4]) - as.numeric(num.vec[i, 2]) )/ as.numeric(num.vec[i, 2])
    num.vec[i, 7] <- (as.numeric(num.vec[i, 5]) - as.numeric(num.vec[i, 3]))/as.numeric(num.vec[i, 3])
  }
  num.vec <- num.vec[complete.cases(num.vec), ]
  colnames(num.vec) <- c("Sample", "Alpha", "Sigma", "EstimateAlpha", "EstimateSigma", "RelativeErrorAlpha", "RelativeErrorSigma")
  num.vec <- num.vec[1:1000, ]
  return(data.frame(num.vec))
}

l <- rbind(sim.MME.case3.alpha.sigma.unknown(alpha = 0.7, sigma = 0.5, num = 50), sim.MME.case3.alpha.sigma.unknown(alpha = 0.7, sigma = 0.5, num = 150), sim.MME.case3.alpha.sigma.unknown(alpha = 0.7, sigma = 0.5, num = 300),
           sim.MME.case3.alpha.sigma.unknown(alpha = 1, sigma = 1, num = 50), sim.MME.case3.alpha.sigma.unknown(alpha = 1, sigma = 1, num = 150), sim.MME.case3.alpha.sigma.unknown(alpha = 1, sigma = 1, num = 300),
           sim.MME.case3.alpha.sigma.unknown(alpha = 2, sigma = 5, num = 50), sim.MME.case3.alpha.sigma.unknown(alpha = 2, sigma = 5, num = 150), sim.MME.case3.alpha.sigma.unknown(alpha = 2, sigma = 5, num = 300),
           sim.MME.case3.alpha.sigma.unknown(alpha = 1, sigma = 0.5, num = 50), sim.MME.case3.alpha.sigma.unknown(alpha = 2, sigma = 0.5, num = 50), sim.MME.case3.alpha.sigma.unknown(alpha = 1, sigma = 0.5, num = 150),
           sim.MME.case3.alpha.sigma.unknown(alpha = 2, sigma = 0.5, num = 150), sim.MME.case3.alpha.sigma.unknown(alpha = 1, sigma = 0.5, num = 300), sim.MME.case3.alpha.sigma.unknown(alpha = 2, sigma = 0.5, num = 300),
           sim.MME.case3.alpha.sigma.unknown(alpha = 0.7, sigma = 1, num = 50), sim.MME.case3.alpha.sigma.unknown(alpha = 2, sigma = 1, num = 50), sim.MME.case3.alpha.sigma.unknown(alpha = 0.7, sigma = 1, num = 150),
           sim.MME.case3.alpha.sigma.unknown(alpha = 2, sigma = 1, num = 150), sim.MME.case3.alpha.sigma.unknown(alpha = 0.7, sigma = 1, num = 300), sim.MME.case3.alpha.sigma.unknown(alpha = 2, sigma = 1, num = 300),
           sim.MME.case3.alpha.sigma.unknown(alpha = 0.7, sigma = 5, num = 50), sim.MME.case3.alpha.sigma.unknown(alpha = 1, sigma = 5, num = 50), sim.MME.case3.alpha.sigma.unknown(alpha = 0.7, sigma = 5, num = 150),
           sim.MME.case3.alpha.sigma.unknown(alpha = 1, sigma = 5, num = 150), sim.MME.case3.alpha.sigma.unknown(alpha = 0.7, sigma = 5, num = 300), sim.MME.case3.alpha.sigma.unknown(alpha = 1, sigma = 5, num = 300))

save(l, file = "MME.case3.alpha.sigma.1000.RData")

library(ggplot2)
# Plots for alpha
# Relative error graphs
gh <- ggplot(l, aes(x = RelativeErrorAlpha)) + geom_histogram(aes(y = ..density..)) + facet_grid( Alpha ~ Sigma + Sample, scales = "free", labeller = label_both) + theme_classic()
gh <- gh + labs(x = "Relative Error")
gh

# Boxplots
th <- ggplot(l, aes(x = as.factor(Sample), y = EstimateAlpha)) + geom_boxplot() + geom_hline(aes(yintercept = Alpha), color = "red") +
  facet_wrap(Alpha ~ Sigma, scales = "free", labeller = label_both) + labs(x = "Sample size", y = expression(paste("Estimate ", alpha))) + theme_classic()
th

# Plots for sigma
gh <- ggplot(l, aes(x = RelativeErrorSigma)) + geom_histogram(aes(y = ..density..)) + facet_grid(Alpha ~ Sigma + Sample, scales = "free", labeller = label_both) + theme_classic()
gh <- gh + labs(x = "Relative Error")
gh

# Boxplots
th <- ggplot(l, aes(x = as.factor(Sample), y = EstimateSigma)) + geom_boxplot() + geom_hline(aes(yintercept = Sigma), color = "red") +
  facet_wrap(Alpha ~ Sigma, scales = "free", labeller = label_both) + labs(x = "Sample size", y = expression(paste("Estimate ", sigma))) + theme_classic()
th

# Generating table
tab1 <- l[which(l$Alpha == 0.7 & l$Sigma == 0.5 & l$Sample == 50), ]
(esti_alpha <- round(mean(tab1$EstimateAlpha), 4 ))
(esti_sigma <- round(mean(tab1$EstimateSigma), 4))
(MSE <- round(mean( (tab1$Alpha - tab1$EstimateAlpha)^2 ), 4 ))
(MSE <- round(mean( (tab1$Sigma - tab1$EstimateSigma)^2 ), 4 ))

tab1 <- l[which(l$Alpha == 2 & l$Sigma == 5 & l$Sample == 300), ]
(esti_alpha <- round(mean(tab1$EstimateAlpha), 4 ))
(esti_sigma <- round(mean(tab1$EstimateSigma), 4))
(MSE <- round(mean( (tab1$Alpha - tab1$EstimateAlpha)^2 ), 4 ))
(MSE <- round(mean( (tab1$Sigma - tab1$EstimateSigma)^2 ), 4 ))
