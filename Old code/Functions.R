library(mvtnorm)
library(TDA)
library(plyr)
library(caret)
library(truncnorm)
library(LaplacesDemon)
library(coda)
library(abind)
library("cowplot")
library("remotes")
library(compiler)
library(latex2exp)
library("expint")
library("pgnorm")
library("gnorm")
library("dplyr")
library("data.table")
library("egg")
library("latex2exp")
library("sirt")
library(LaplacesDemon)
library("FNN")
library("tram")
library("RColorBrewer")
library(abind)
library(bindata)
library("pracma")
library(bayestestR)


#library("waddR")

#library(waddR)
#stopQuietly <- function(...) {
#  blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "));
#  stop(simpleError(blankMsg));
#} # stopQuietly()
#stopQuietly()

multi_chain <- function(N_sim,burn_in,X, y,initial_theda,true_theta,Lp,M_iter,range,step,times){
  initialtime <- proc.time()
  Beta_chain <- list()
  Lp_chain <- list()
  
  for(i in 1:times){
    result <- Metropolis_Lp_gibbssampler(N_sim,burn_in,
                                         X, y,initial_theda,
                                         true_theta,Lp,M_iter,range,step)
    
    a <-  list()
    a[[paste("chain",i,sep ="")]] <- result$chain
    b <-  list()
    b[[paste("LP_chain",i,sep ="")]] <- result$Lp
    
    Beta_chain <- c(Beta_chain,a)
    Lp_chain <- c(Lp_chain,b)
  }
  Posterior.beta <- data.frame()
  Lp <- c()
  for (i in 1:times) {  
    a <- data.frame(colMeans(Beta_chain[[i]][-(1:burn_in), ]))
    b <- mean(Lp_chain[[i]][-(1:burn_in)])
    Posterior.beta <- data.frame(rbind(Posterior.beta,t(a)))
    row.names(Posterior.beta) <- NULL
    Lp <- c(Lp,b)
  }
  
  temp_array <- abind(Beta_chain,along=3)
  beta.mean.chain <- apply(temp_array, 1:2, mean)
  temp_array <- abind(Lp_chain,along=2)
  Lp.mean.chain <- rowMeans(temp_array)
  Beta.mean <- colMeans(Posterior.beta)
  Lp.mean <- mean(Lp)
  modeling.time <- proc.time() - initialtime 
  
  result <- list("modeling.time"=modeling.time, "Beta_chain"=Beta_chain,"Lp_chain"=Lp_chain,
                 "beta.mean.chain"=beta.mean.chain,"Lp.mean.chain"=Lp.mean.chain,
                 "Posterior.beta"=Posterior.beta,"Lp"=Lp,
                 "Beta.mean"=Beta.mean,"Lp.mean"=Lp.mean)
  return(result)
}
#
#generate.data <- function(N,D){
  #x is generated from a multivariate normal distribution
#  meanvector <- rnorm(D-1,0,5)

 # a <- rmvnorm(n=N,mean=meanvector,sigma=diag(4,D-1,D-1))
  # Create n x D design matrix
  #X <- matrix(c(rep(1, N), a,b), ncol = D)
#  X <- matrix(c(rep(1, N), a), ncol = D)
  
#  return(X)
#}

generate.data <- function(N,D){
  #x is generated from a multivariate normal distribution
  b <- D-1
  sigma <-matrix(rep(NA,b*b),nrow = b)
  for (i in 1:b) {
    for (j in 1:b) {
      sigma[i,j] <- 2*(0.5^abs(i-j))
    }
  }
  
  #a <- rmvnorm(n=N,mean=c(rep(-5,(D-1)/5),rep(5,(D-1)/5),
   #                       rep(0,(D-1)/5)),sigma=sigma)
  
  a <- rmvnorm(n=N,mean=c(rep(-2,(D-1)/5),rep(2,(D-1)/5),
                          rep(-3,(D-1)/5),rep(3,(D-1)/5),
                          rep(0,(D-1)/5)),sigma=sigma)
  # Create n x D design matrix
  #X <- matrix(c(rep(1, N), a,b), ncol = D)
  X <- matrix(c(rep(1, N), a), ncol = D)
  
  return(X)
}

generate.data0 <- function(N,D){
  N1 <- N/2
  #x is generated from a multivariate normal distribution
  b <- D-1
  sigma <-matrix(rep(NA,b*b),nrow = b)
  for (i in 1:b) {
    for (j in 1:b) {
      sigma[i,j] <- 2*(0.5^abs(i-j))
    }
  }
  
  #a <- rmvnorm(n=N,mean=c(rep(-5,(D-1)/5),rep(5,(D-1)/5),
  #                      rep(0,(D-1)/5)),sigma=sigma)
  
  #a <- rmvnorm(n=N1,mean=c(rep(-3,b)),sigma=sigma)
  a <- rmvt(n=N1,mu=c(rep(-4,b)), S=sigma, df=4)
  
  # Create n x D design matrix
  #X <- matrix(c(rep(1, N), a,b), ncol = D)
  X <- matrix(c(rep(1, N1), a), ncol = D)
  
  return(X)
}

generate.data1 <- function(N,D){
  N1 <- N/2
  #x is generated from a multivariate normal distribution
  b <- D-1
  sigma <-matrix(rep(NA,b*b),nrow = b)
  for (i in 1:b) {
    for (j in 1:b) {
      sigma[i,j] <- 4*(0.5^abs(i-j))
    }
  }
  
  #a <- rmvnorm(n=N,mean=c(rep(-5,(D-1)/5),rep(5,(D-1)/5),
  #                       rep(0,(D-1)/5)),sigma=sigma)
  
  #a <- rmvnorm(n=N1,mean=c(rep(3,b)),sigma=sigma)
  a <- rmvt(n=N1,mu=c(rep(8,b)), S=sigma, df=8)
  
  # Create n x D design matrix
  #X <- matrix(c(rep(1, N), a,b), ncol = D)
  X <- matrix(c(rep(1, N1), a), ncol = D)
  
  return(X)
}
# 
# generate.data.t.distribution <-function(N,D){
#   b <- (D-21)/2
#   sigma <-matrix(rep(NA,b*b),nrow = b)
#   for (i in 1:b) {
#     for (j in 1:b) {
#       sigma[i,j] <- 2*(0.5^abs(i-j))
#     }
#   }
#   a1 <- rmvt(n=N,mu=c(rep(-5,b/2),rep(5,b/2)), S=sigma, df=1)
#   a2 <- rmvc(n=N,mu=c(rep(-5,b/2),rep(5,b/2)), S=sigma)
#   a3 <- rmvbin(N, margprob=rep(0.5,D-b*2-1))
#   a <- cbind(a1,a2,a3)
#   #a <- a[,sample(1:(D-1))]
#   X <- matrix(c(rep(1, N),a), ncol = D)
#   return(X)
# }

# generate.data.cauthy <-function(N,D){
#   b <- D-1
#   sigma <-matrix(rep(NA,b*b),nrow = b)
#   for (i in 1:b) {
#     for (j in 1:b) {
#       sigma[i,j] <- 2*(0.5^abs(i-j))
#     }
#   }
#   a <- rmvc(n=N,mu=rep(0,b), S=sigma)
#   X <- matrix(c(rep(1, N), a), ncol = D)
#   return(X)
#   
# }

# Holmes_Gibbs_Sampler<- function(D=D, N=N, N_sim,burn_in, 
#                                 X, y,initial_theta, true_theta=true_theta)
# {
#   theta_0 <- rep(0, D)
#   Q_0 <- diag(10, D)
#   theta <- initial_theta
#   z <- X%*%theta
#   
#   chain <- matrix(0, nrow = N_sim, ncol = D)
#   theta_chain <- matrix(0, nrow = N_sim, ncol = D)
#   theta_chain[1,] <- theta
#   chainlist <- list()
#   prec_0 <- solve(Q_0)
#   V <- solve(prec_0 + crossprod(X, X))
#   y.probit <- y
#   N1 <- sum(y.probit)
#   N0 <- N-N1
#   
#   # Compute the martix S = VX'
#   S <- tcrossprod(V, X)
#   
#   
#   h <- vector(mode = "numeric", length = N)
#   w <- vector(mode = "numeric", length = N)
#   u <- vector(mode = "numeric", length = N)
#   
#   
#   for (j in 1:N){
#     # h stores the diagonal elements of the hat matrix (XS = XVX')
#     h[j] <- X[j, ] %*% S[, j]
#     w[j] <- h[j] / (1 - h[j])
#     u[j] <- w[j] + 1
# 
#   }
#   mu_z <- X %*% theta
#   # Initialize latent variable Z, from truncated normal
#   z[y.probit == 0] <- rtruncnorm(N0, mean = mu_z[y.probit == 0]   , sd = 1, a = -Inf, b = 0)
#   z[y.probit == 1] <- rtruncnorm(N1, mean = mu_z[y.probit == 1], sd = 1, a = 0, b = Inf)
#   
#   # ---------------------------------
#   # Gibbs sampling algorithm
#   # ---------------------------------
#   
#   # Compute the conditional mean of \theta
#   M <- as.vector(S %*% z)
#   initialtime <- proc.time()  # record classification time 
#   
#   
#   for (t in 2:N_sim) {
#     for (j in 1:N){
#       # Store the old value of z
#       z_old <- z[j]
#       
#       # Update mean of latent variable z_i
#       m <- X[j, ] %*% M
#       m <- m - w[j] * (z[j] - m)
#       
#       # Draw latent variable z from full conditional: z_j | z_-j, y, X
#       if (y.probit[j] == 0)
#         z[j] <- rtruncnorm(1, mean = m, sd = u[j], a = -Inf, b = 0)
#       else
#         z[j] <- rtruncnorm(1, mean = m, sd = u[j], a = 0, b = Inf)
#       
#       # Update the posterior mean M
#       M <- as.vector(M + (z[j] - z_old) %*% S[ ,j])
#     }
#     
#     # Posterior of M | Z
#     theta_chain[t, ] <- c(rmvnorm(1, M, V))
#     if(as.integer(100 * (t/N_sim)) == as.numeric(100 * (t/N_sim)))
#       print(paste("Progress %: ", 100 * (t/N_sim)))
#     
#   }
#   post.mean <- colMeans(theta_chain[-(1:burn_in), ])
#   modeling.time <- proc.time() - initialtime 
#   chainlist <- list( "chain"= theta_chain, "modeling.time"=modeling.time,
#                      "true_theta"=true_theta, "post.mean"=post.mean) # record classification time 
#   return(chainlist)
# }



gibbssampler<- function(D=D, N=N, N_sim,burn_in, 
                        X, y,initial_theta, true_theta=true_theta)
  {
  theta_0 <- rep(0, D)
  Q_0 <- diag(10, D)
  theta <- initial_theta
  #z <- rep(0, N)
  z <- X %*% theta
  theta_chain <- matrix(0, nrow = N_sim, ncol = D)
  theta_chain[1,] <- theta
  chainlist <- list()
  prec_0 <- solve(Q_0)
  V <- solve(prec_0 + crossprod(X, X))
  y.probit <- y
  N1 <- sum(y.probit)
  N0 <- N-N1
  initialtime <- proc.time()  # record classification time 
  for (t in 2:N_sim) {
    # Update Mean of z
    mu_z <- X %*% theta
    # Draw latent variable z from its full conditional: z | \theta, y, X
    z[y.probit == 0] <- rtruncnorm(N0, mean = mu_z[y.probit == 0], sd = 1, a = -Inf, b = 0)
    z[y.probit == 1] <- rtruncnorm(N1, mean = mu_z[y.probit == 1], sd = 1, a = 0, b = Inf)
    
    # Compute posterior mean of theta
    M <- V %*% (prec_0 %*% theta_0 + crossprod(X, z))
    # Draw variable \theta from its full conditional: \theta | z, X
    theta <- c(rmvnorm(1, M, V))
    # Store the \theta draws
    theta_chain[t, ] <- theta
    if(as.integer(100 * (t/N_sim)) == as.numeric(100 * (t/N_sim)))
      print(paste("Progress %: ", 100 * (t/N_sim)))
    
  }
  post.mean <- colMeans(theta_chain[-(1:burn_in), ])
  modeling.time <- proc.time() - initialtime 
  chainlist <- list( "chain"= theta_chain, "modeling.time"=modeling.time,
                     "true_theta"=true_theta, "post.mean"=post.mean)                  # record classification time 
  return(chainlist)
}

# p_scale <- function(p){
#   result <- (1/p)^(1/p)
#   return(result)
# }
p_scale <- function(p){
  result <- p^(1/p)
  return(result)
}

truncated_gnorm <- function(n, range, mu,alpha,beta) {
  
  # range is a vector of two values
  
  #F.a <- pgnorm(min(range), alpha=alpha,beta=beta, mu = mu)
  #F.b <- pgnorm(max(range), alpha=alpha,beta=beta, mu = mu)
  F.a <- pgnorm(min(range), alpha=alpha,beta=beta, mu = mu)
  F.b <- pgnorm(max(range), alpha=alpha,beta=beta, mu = mu)
  
  u <- runif(n, min = F.a, max = F.b)
  u[u<0.000001] <- 0+0.000001
  u[u>1-0.000001] <- 1-0.000001
  # 
  #if(is.na(u[1])) stopQuietly("NA")
  #Q
  result <- qgnorm(u, mu = mu,alpha = alpha,beta=beta)
  return(result)
  
}

# range <- range(-Inf,0)
# mu <- 0
# F.a <- pgnorm(min(range), alpha= Phi_p(Lp),beta=Lp, mu = mu)
# F.b <- pgnorm(max(range), alpha=Phi_p(Lp),beta=Lp, mu = mu)
#  
# u <- runif(10, min = F.a, max = F.b)
# u
# result <- qgnorm(u, mu = mu,alpha = alpha,beta=beta)
# 
# 
# 
# rgnorm(u, mu = mu,alpha = alpha,beta=beta)
# F.a <- pgnorm(-Inf, alpha=Phi_p(2),beta=2, mu = 0)
# F.b <- pgnorm(0, alpha=Phi_p(2),beta=2, mu =0)
# F.a
# F.b
# 
# if(F.a==F.b)d=runif(1,1e-20,1e-10)
# d
# u <- runif(10, min = min(F.a,F.b), max = max(F.a,F.b))
# u
# qgnorm(u, mu = 3,alpha = Phi_p(7),beta=7)
# 
# 1e-20
# 
# F.a==F.b
# F.b<-pgnorm(0, alpha=Phi_p(7),beta=7, mu =-2)
# pgnorm(Inf, alpha=Phi_p(2),beta=2, mu =2)
# 
# u <- runif(10, min = 0, max = 0)
# u
# qgnorm(u, mu = 2,alpha = Phi_p(5),beta=5)
# 
# 
# truncated_gnorm(10,mu=2,
#                 range=c(-Inf,0),alpha = Phi_p(7),beta = 7)

Lp_gibbssampler<- function(N_sim,burn_in, 
                        X, y,initial_theda, true_theta=true_theta,Lp)
{ N <- nrow(X)
  D <- ncol(X)
  theta_0 <- rep(0, D)
  Q_0 <- diag(10, D)
  theta <- initial_theda
  z <- rep(0, N)
  theta_chain <- matrix(0, nrow = N_sim, ncol = D)
  chainlist <- list()
  prec_0 <- solve(Q_0)
  # V <- solve(prec_0 + crossprod(X, X))
  V1 <- solve(crossprod(X,X),tol = 1e-25)
  V <- solve(prec_0 + crossprod(X, X))
  tau_p <- (1/Lp)^(2/Lp)*gamma(3/Lp)/gamma(1/Lp)
  basis <- chol(solve(t(X)%*%X))
  V2 <- (basis%*% t(basis))
  N1 <- sum(y)
  N0 <- N-N1
  initialtime <- proc.time()  # record classification time 
  for (t in 2:N_sim) {
    # Update Mean of z
    mu_z <- X %*% theta
    # Draw latent variable z from its full conditional: z | \theta, y, X
    z[y == 0] <- truncated_gnorm(N0,mu=mu_z[y == 0],
                                       range=c(-Inf,0),alpha = p_scale(Lp),beta = Lp)
    z[y == 1] <- truncated_gnorm(N1,mu=mu_z[y == 1],
                                        range=c(0,Inf),alpha = p_scale(Lp),beta = Lp)

    
    #negative.z <- rgnorm(N0, mu = mu_z[y == 0],alpha = p_scale(Lp),beta=Lp)
    #negative.z[sign(negative.z) == 1] <- negative.z[sign(negative.z) == 1]*-1
    #z[y == 0] <- negative.z
    
    #positive.z <- rgnorm(N1, mu = mu_z[y == 1],alpha = p_scale(Lp),beta=Lp)
    #positive.z[sign(positive.z) == -1] <- positive.z[sign(positive.z) ==-1]*-1
    #z[y == 1] <- positive.z

    # Compute posterior mean of theta
    model <- lq_fit(y=z, X=X, pow=Lp,est_pow=F)
    M <- model$coefficients
    # M <- V %*% (prec_0 %*% theta_0 + crossprod(X, z))
    # Draw variable \theta from its full conditional: \theta | z, X
    # theta <- rgnorm(D, mu=M, alpha =p_scale(Lp) ,beta = Lp)
    
    for (i in 1:D) {
      theta[i] <- rgnorm(1, mu=M[i], alpha =V2[i,i] , beta = Lp)
      # theta[i] <- rgnorm(1, mu=M[i], alpha =p_scale(Lp) ,beta = Lp)

    }
    # theta <- c(rmvnorm(1, M, V1))
    # Store the \theta draws
    theta_chain[t, ] <- theta
    if(as.integer(100 * (t/N_sim)) == as.numeric(100 * (t/N_sim)))
      print(paste("Progress %: ", 100 * (t/N_sim)))
  }
  post.mean <- colMeans(theta_chain[-(1:burn_in), ])
  modeling.time <- proc.time() - initialtime 
  chainlist <- list( "chain"= theta_chain, "modeling.time"=modeling.time,
                    "true_theta"=true_theta, "post.mean"=post.mean)                    # record classification time 
  return(chainlist)
}

latent_z <- function(X,theta,y,Lp){
  N <- length(y)
  N1 <- sum(y)
  N0 <- N-N1
  z <- rep(0, N)
  
  mu_z <- X %*% theta
  # Draw latent variable z from its full conditional: z | \theta, y, X
  z[y == 0] <- truncated_gnorm(N0,mu=mu_z[y == 0],
                               range=c(-Inf,0),alpha = p_scale(Lp),beta = Lp)
  z[y == 1] <- truncated_gnorm(N1,mu=mu_z[y == 1],
                               range=c(0,Inf),alpha = p_scale(Lp),beta = Lp)
  return(z)
  
}

calculate_and_plot_posterior <- function(chain) {
  # Calculate the density
  d <- density(chain)
  
  # Plot the posterior distribution
  plot(d, main = "Posterior distribution", xlab = "Values", ylab = "Density")
  
  # create 
  f <- approxfun(d$x, d$y)
  
  # Setting upper limit and lower limit
  lower_limit <- min(d$x)
  upper_limit <- max(d$x)
  
  # Calculate the integral
  result <- tryCatch({
    integrate(f, lower_limit, upper_limit, subdivisions = 1000)
  }, warning = function(w) {
    return("Warning: The integral calculation might not be accurate.")
  }, error = function(e) {
    return("Error: The integral could not be calculated.")
  })
  df <- data.frame(Values = d$x, Density = d$y)
  # Print out integral
return(list(posterior_data=df,result))
  }



proposal <- function(Lp,range,L){
  p <- runif(1,min=max(min(range),(Lp-L)),max=min(max(range),(Lp+L)))
  return(p)
}


Metropolis <- function(X, M_iter, theta, Lp, tol=1e-40,range,L,z){
  Xb <- X%*%theta
  ehat <- z - Xb
  
  u <- runif(M_iter)
  for(i in 2:M_iter){
    proposed <- proposal(Lp,range,L)
    
    R0 <- sum(log(dgnorm(ehat,alpha = p_scale(proposed),beta = proposed)+tol))
    R1 <- sum(log(dgnorm(ehat,alpha = p_scale(Lp),beta = Lp)+tol))
    #R0 <- sum((pgnorm(Xb,alpha = p_scale(proposed),beta = proposed)))
    #R1 <- sum((pgnorm(Xb,alpha = p_scale(Lp),beta = Lp)))
    adj<- max(R0,R1)
    R  <- exp(R0-adj)/exp(R1-adj)
    #R  <- exp(R0)/exp(R1)
    #R  <-R0/(3*R1)
    
    if(u[i] < R) Lp = proposed
    }
  return(Lp)
}

llk <- function(X,b,p,y){
  Xb <- X%*%b
  tol <- 0.001
  pi <- pgnorm(Xb, alpha = p_scale(p),beta = p)
  # a[y==0] <- 1-a[y==0]
  pi[pi==0] <- tol
  pi[pi==1] <- 1-tol
  llk <- sum(y*log(pi)+(1-y)*log(1-pi))
  return(llk)
}

llk_logit <- function(X,b,y){
  Xb <- X%*%b
  tol <- 0.001
  pi <- exp(Xb)/(1+exp(Xb))
  # a[y==0] <- 1-a[y==0]
  pi[pi==0] <- tol
  pi[pi==1] <- 1-tol
  llk <- sum(y*log(pi)+(1-y)*log(1-pi))
  return(llk)
}

llk_clog <- function(X,b,y){
  Xb <- X%*%b
  tol <- 0.001
  pi <- 1-exp(-exp(Xb))  # a[y==0] <- 1-a[y==0]
  pi[pi==0] <- tol
  pi[pi==1] <- 1-tol
  llk <- sum(y*log(pi)+(1-y)*log(1-pi))
  return(llk)
}

Bayes_factor <- function(X,b,p,y){
  XB  <-X*b
  # Z <- -(2*y-1)*X
  #Z <- X
  Xb <- X%*%b
  tol <- 0.0000001
  a <- pgnorm(Xb, alpha = p_scale(p),beta = p)
  a[y==0] <- 1-a[y==0]
  a[a==0] <- tol
   # llk1 <- t(t(log(a))+log(pnorm(b,sd=5)))
  # llk1 <- -(log(a)+sum(log(pgnorm(b, alpha = p_scale(p),beta = p))))
  BF <- -(log(a)+sum(log(pnorm(b,sd=1))))
  return(BF)
}


Metropolis_Lp_gibbssampler <- function(N_sim,burn_in, 
                                       X, y,initial_theda, true_theta=true_theta,Lp,M_iter,range,step){
  D <- ncol(X)
  N <- nrow(X)
  L <- step
  theta_0 <- rep(0, D)
  Q_0 <- diag(10, D)
  theta <- initial_theda
  z <- rep(0, N)
  theta_chain <- matrix(0, nrow = N_sim, ncol = D)
  Lp_chain <- rep(0, N_sim)
  chainlist <- list()
  prec_0 <- solve(Q_0)
  V <- solve(prec_0 + crossprod(X, X))
  V1 <- solve(crossprod(X,X))
  tau_p <- (Lp)^(2/Lp)*gamma(3/Lp)/gamma(1/Lp)
  basis <- chol(solve(t(X)%*%X))
  V2 <- tau_p*(basis%*% t(basis))
  N1 <- sum(y)
  N0 <- N-N1
  initialtime <- proc.time()  # record classification time 
  for (t in 2:N_sim) {
    # Update Mean of z
    z <- latent_z(X,theta,y,Lp)
    #negative.z <- rgnorm(N0, mu = mu_z[y == 0],alpha = p_scale(Lp),beta=Lp)
    #negative.z[sign(negative.z) == 1] <- negative.z[sign(negative.z) == 1]*-1
    #z[y == 0] <- negative.z
    
    #positive.z <- rgnorm(N1, mu = mu_z[y == 1],alpha = p_scale(Lp),beta=Lp)
    #positive.z[sign(positive.z) == -1] <- positive.z[sign(positive.z) ==-1]*-1
    #z[y == 1] <- positive.z
     model <- lq_fit(y=z, X=X, pow=Lp, est_pow=F)
    # model <- lm(z~X)
     M <- model$coefficients

    # Compute posterior mean of theta
    # M <- solve(crossprod(X,X))%*%crossprod(X,z)
    # M <- V %*% (prec_0 %*% theta_0 + crossprod(X, z))
    # Draw variable \theta from its full conditional: \theta | z, X
    #theta <- c(rgnorm(1, mu=M, alpha = ,beta = Lp))
    for (i in 1:D) {
      # theta[i] <- rgnorm(1, mu=M[i], alpha =V2[i,i] , beta = Lp)
      theta[i] <- rgnorm(1, mu=M[i], alpha =V2[i,i], beta = Lp)

    }
    # theta <- c(rmvnorm(1, M, V))
    # Store the \theta draws
    theta_chain[t, ] <- theta
    Lp <- Metropolis(X=X,theta = theta,Lp=Lp,range = range,L=L,M_iter=M_iter,z=z)
    
    Lp_chain[t] <- Lp
    if(as.integer(100 * (t/N_sim)) == as.numeric(100 * (t/N_sim)))
      print(paste("Progress %: ", 100 * (t/N_sim)))
  }
  post.mean <- colMeans(theta_chain[-(1:burn_in), ])
  p.mean <- mean(Lp_chain[-(1:burn_in)])
  modeling.time <- proc.time() - initialtime 
  chainlist <- list( "chain"= theta_chain, "modeling.time"=modeling.time,
                     "true_theta"=true_theta, "post.mean"=post.mean,"Lp"=Lp_chain,"Lp.mean"=p.mean)                    # record classification time 
  return(chainlist)
}

ucomplex <- function(Xb){
  u <- norm(as.matrix(Xb[Xb>0]),type = "1")/norm(as.matrix(Xb[Xb<0]),type="1")
  return(u)
}
qrhat <- function(x, intercept = TRUE)
{
  if(is.qr(x)) n <- nrow(x$qr)
  else {
    if(intercept) x <- cbind(1, x)
    n <- nrow(x)
    x <- qr(x)
  }
  rowSums(qr.qy(x, diag(1, nrow = n, ncol = x$rank))^2)
} 

sampling.size <- function(eps,u,N,D){
  size <- ((u*(D^2))/(eps^2))*log(u*D*N)
  return(size)
}

sample_gibbssampler<- function(D=D, N=N, N_sim,burn_in, data=data,initial_theta,
                         true_theta=true_theta){ 
  X <- as.matrix(data$X_tilda)
  y <- data$y_tilda
  theta_0 <- rep(0, D)
  Q_0 <- diag(10, D)
  theta <- initial_theta
  z <- rep(0, N)
  theta_chain <- matrix(0, nrow = N_sim, ncol = D)
  chainlist <- list()
  prec_0 <- solve(Q_0)
  V <- solve(prec_0 + crossprod(X, X))
  y.probit <- y
  N1 <- sum(y.probit)
  N0 <- N-N1
  initialtime <- proc.time()  # record classification time 
  for (t in 2:N_sim) {
    # Update Mean of z
    mu_z <- X %*% theta
    # Draw latent variable z from its full conditional: z | \theta, y, X
    z[y.probit == 0] <- rtruncnorm(N0, mean = mu_z[y.probit == 0], sd = 1, a = -Inf, b = 0)
    z[y.probit == 1] <- rtruncnorm(N1, mean = mu_z[y.probit == 1], sd = 1, a = 0, b = Inf)
    
    # Compute posterior mean of theta
    M <- V %*% (prec_0 %*% theta_0 + crossprod(X, z))
    # Draw variable \theta from its full conditional: \theta | z, X
    theta <- c(rmvnorm(1, M, V))
    # Store the \theta draws
    theta_chain[t, ] <- theta
    if(as.integer(100 * (t/N_sim)) == as.numeric(100 * (t/N_sim)))
    print(paste("Progress %: ", 100 * (t/N_sim)))
    
  }
  post.mean <- colMeans(theta_chain[-(1:burn_in), ])
  modeling.time <- proc.time() - initialtime 
  chainlist <- list("modeling.time"=modeling.time, "chain"= theta_chain, 
                    "true_theta"=true_theta, "post.mean"=post.mean)                    # record classification time 
  return(chainlist)
}
subsample <- function(N,frac,prob,data){
  sample_data <- as.data.frame(data)
  sample_data <- cbind(sample_data,index=1:N)
  sampled_data <- sample_n(sample_data,size = N*frac,
                    replace = FALSE, weight =prob)
  sampling<- sampled_data$index
  X_tilda <- sampled_data[,1:D]
  pi_tilda <- sampled_data$pi.probit
  y_tilda <- sampled_data$y.probit
  N1_tilda <- sum(y_tilda)
  N0_tilda <- length(pi_tilda)-N1_tilda
  result <- list(X_tilda=X_tilda,pi_tilda=pi_tilda,y_tilda=y_tilda,N1_tilda=N1_tilda,
                 N0_tilda=N0_tilda,sampling=sampling)
  return(result)
}


Compute_sketch <- function(X,Lp,sketch_size){
  N <- nrow(X)
  D <- ncol(X)
  B_i <- sample(1:sketch_size,N,replace = T)
  sigma_i <- sample(c(-1,1),N,replace = T)
  X_b <- matrix(0,nrow = sketch_size,ncol = D)
  lambda_i <- rexp(N,1)
  if(Lp != 2) sigma_i <- sigma_i/(lambda_i^(1/Lp))
  for (i in 1:N) {
    X_b[B_i[i],]=X_b[B_i[i],]+sigma_i[i]*X[i,]
    
  }
  return(X_b)
}

Uniform_coreset <- function(sketch_size,coreset_size,X,Lp){
  N <- nrow(X)
  D <- ncol(X)
  scores <- runif(N)
  sample_index <- sample(1:N,size = coreset_size,
                         replace = FALSE, prob  =scores)
  return(sample_index)
  
}
lp_leverage<- function(x,Lp){
  result<-apply(Matrix::qr.Q(qr(X)),MARGIN = 1,Norm,p = Lp)^(1/2)
  return(result)
}

Lp_coreset <- function(sketch_size,coreset_size,X,Lp){
  N <- nrow(X)
  D <- ncol(X)
  scores <- lp_leverage(x=X,Lp=Lp)
  sample_index <- sample(1:N,size = coreset_size,
                         replace = FALSE, prob  =scores)
  return(sample_index)
  
}

Compute_coreset <- function(sketch_size,coreset_size,X,Lp){
  N <- nrow(X)
  D <- ncol(X)
  X_prime <- matrix(0,nrow = coreset_size,ncol = D)
  X_b <- Compute_sketch(sketch_size=sketch_size,X=X,Lp=Lp)
  R_matrix <- qr.R(qr(X_b))
  if (Lp==2 && log(N) < D) G<-diag(rnorm(D,0,1/log(N))) else G <- diag(D)
  qi <- c()
  XR <-  X%*%(solve(R_matrix)%*%G)
  
  for (i in 1:N) {
    qi[i] <- Norm(XR[i,],p = Lp)^Lp
  }
  scores <- qi +1/N
  sample_index <- sample(1:N,size = coreset_size,
                         replace = FALSE, prob = scores)
  return(sample_index)
}

Sensitivity_score <- function(sketch_size,coreset_size,X,Lp){
  N <- nrow(X)
  D <- ncol(X)
  X_prime <- matrix(0,nrow = coreset_size,ncol = D)
  X_b <- Compute_sketch(sketch_size=sketch_size,X=X,Lp=Lp)
  R_matrix <- qr.R(qr(X_b))
  if (Lp==2 && log(N) < D) G<-diag(rnorm(D,0,1/log(N))) else G <- diag(D)
  qi <- c()
  XR <-  X%*%(solve(R_matrix)%*%G)
  
  for (i in 1:N) {
    qi[i] <- Norm(XR[i,],p = Lp)^Lp
  }
  scores <- qi

  return(scores)
}

one_shot_coreset <- function(sketch_size,coreset_size,X,Lp_max){
  N <- nrow(X)
  delta <- 1/log(N)
  score <- rep(0,N)
  l <- round(log(Lp_max)/log(1+delta))
  for (i in 0:l) {
      p <- (1+delta)^i
      score_i <- Sensitivity_score(sketch_size=sketch_size,coreset_size=coreset_size,X=X,Lp=p)
      score <- cbind(score,score_i)
  }
  score <- rowSums(score)+1/N
  return(score)
}

one_shot_coreset_times <- function(sketch_size,coreset_size,X,Lp_max,times){
  initialtime <- proc.time()
  N <- nrow(X)
  delta <- 1/log(N)
  times <- times
  score <- rep(0,N)
  for (i in 0:times) {
    p <- (1+delta)^i
    score_i <- Sensitivity_score(sketch_size=sketch_size,coreset_size=coreset_size,X=X,Lp=p)
    score <- cbind(score,score_i)
  }
  score <- rowSums(score)+1/N
  modeling.time <- proc.time() - initialtime 
  coreset_list <- list( "score"= score, "modeling.time"=modeling.time)            
  return(coreset_list)
}

gaussian_kernel<- function(x_1,x_2){
  x_1 <- as.matrix(x_1)
  x_2 <- as.matrix(x_2)
  result <- exp(-0.5 * norm(x_1 - x_2) ^ 2)
  return(result)
  
}

polynomial_kernel <- function(x_1, x_2){
  result<- (1 + (x_1%*% x_2))^ 2
  return(result)
}

mmd <- function(sample_x, sample_y){
  m <- length(sample_x)
  n <- length(sample_y)
  
  #sample_x = np.ascontiguousarray(sample_x)
  
  #Computes an estimate of the maximum mean discrepancy between sample_1 and sample_2.
  #See https://arxiv.org/abs/0805.2368 for more info.

  #sample_y = np.ascontiguousarray(sample_y)

sum_x = 0
for (i in 1:m) {
  for (j in 1:m) {
    sum_x = sum_x + gaussian_kernel(sample_x[i], sample_x[j])
  }
}

sum_y = 0
for (i in 1:n) {
  for (j in 1:n) {
    sum_y = sum_y + gaussian_kernel(sample_y[i], sample_y[j])
  }
}
sum_x_y = 0
for (i in 1:m) {
  for (j in 1:n) {
    sum_x_y = sum_x_y + gaussian_kernel(sample_x[i], sample_y[j])
  }
}
result <- sqrt(1 / m ^ 2 * sum_x + 1 / n ^ 2 * sum_y - 2 / (m * n) * sum_x_y)

return(result)
}

generate_p_y <- function(p,Xb,N){
  pi <- pgnorm(Xb,alpha = p_scale(p),beta = p)
  y.p <- rbinom(N,1,pi)
  return(y.p)
}

