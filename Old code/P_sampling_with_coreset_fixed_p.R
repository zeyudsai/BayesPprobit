set.seed(2)
rm(list = ls())
options(scipen = 200)

library(abind)

setwd("/Users/dingzeyu/Documents/文稿 - MacBook Pro/Research Project/Generalized Normal Distribution using MCMC/R_code")
source("Functions.R")

dev.off()
N <- 50000
#D should always be odd number since we use (D-1)/2 to generate the variables
D <- 11
N_sim <- 500
burn_in<- 300
#gengeratedata and true parameterye

#X <- generate.data.cauthy(N,D)
X <- generate.data.cauthy(N,D)
X[X>15] <- sample(-14:14,size = 1)
X[X<(-15)] <- sample(-14:14,size = 1)
#hist(X[,5])
#X[,2:61] <- X[,2:61]/10
X[,2:D] <- X[,2:D]/10
# True values of regression coeffiecients theta
true_theta <- runif(D,-2,2)
#true_theta <- rpois(D,3)
Xb <- X%*%true_theta+rgnorm(N,mu=0,alpha= Phi_p(2),beta=2)
  #rnorm(N,0,sd=1)
  #rt(N,df=10)

  #


hist(Xb)
#Obtain the true probability of pi and the real bernouli distributed Y
pi.logit <- exp(Xb)/(1+exp(Xb))
pi.probit <- pnorm(Xb)
y.logit <- rbinom(n=N, size=1, prob=pi.logit)
#y.logit[is.na(y.logit)] <- 0
y.probit <- rbinom(N, 1, pi.probit)
hist(pi.logit)

hist(pi.probit)

hist(y.probit)
hist(y.logit)
plot(Xb,y.logit)
plot(Xb,pi.logit)
plot(Xb,pi.probit)

plot(Xb,y.probit)


sample_index<- Compute_coreset(sketch_size = 200,coreset_size = 200,X=X,Lp=2)

sample_index

#Using MLE estimator to otain the maximum likelihood result
#Assign a as a symbol states NxxKDxx,which represents the sample size and number of variables.

MLE_model <- glm(y.probit~X[,-1],family = binomial(link = probit))
MLE_logit <- glm(y.logit~X[,-1],family = binomial(link = logit))
MLE_model
MLE_logit
true_theta
MLE_coefficients <- MLE_model$coefficients

#Lp_gibbssampler(D=D, N=N, N_sim=5000,
#                burn_in=4000, X=X, y=y.probit,true_theta = true_theta,Lp=5,initial_theda = MLE_logit$coefficients)
# plot(N50KD11LP_Gibbs2$chain[,1],type="l")
Lp <- c(0.5,1,1.5,2,3,4,5,8)

for (i in Lp) {
  pi_p <- pgnorm(Xb,alpha = Phi_p(i),beta = i)
  assign(paste("pi_",i,sep = ""),pi_p)
  
}
for (i in Lp) {
  y_pi <- generate_p_y(p=i,Xb=Xb,N=N)
  assign(paste("y_p",i,sep = ""),y_pi)
  
}
par(mfrow = c(1, 1))

P_probit <- multi_chain(N_sim=N_sim,burn_in=burn_in,
                        X=X, y=y.probit,initial_theda=MLE_model$coefficients,
                        true_theta=true_theta,Lp=5,M_iter=200,range=c(0.1,5),step=0.05,times = 10)
P_probit$Beta.mean
true_theta
MLE_coefficients
length(P_probit$Lp_chain[[1]])
N_sim
plot(1:N_sim,P_probit$Lp_chain[[1]],type = "l")
P_logit <- multi_chain(N_sim=N_sim,burn_in=burn_in,
                       X=X, y=y.logit,initial_theda=MLE_logit$coefficients,
                       true_theta=true_theta,Lp=5,M_iter=200,range=c(0.1,10),step=0.05,times = 10)
P_logit
P_logit$Beta.mean
MLE_logit$coefficients
plot(1:N_sim,P_logit$Lp_chain[[1]],type = "l")

P_0.5 <- multi_chain(times=10,N_sim=N_sim,burn_in=burn_in,
                     X=X, y=y_p0.5,initial_theda=MLE_model$coefficients,
                     true_theta=true_theta,Lp=5,M_iter=200,range=c(0.1,5),step=0.05)


P_logit
P_logit$chain
P_logit$Lp.mean
plot(P_logit$Lp.mean.chain[1:N_sim],type = "l")
plot(P_logit$Lp.mean.chain[-1],type = "l",main = "Estimated MCMC chain of p with metropolis-hasting", xlab = "iteration",ylab = "p" )
Msubtitle <- "P=logit"
mtext(side = 3, line = 0.4, Msubtitle)
P_logit

P_0.5
P_0.5$Lp.mean

plot(P_0.5$Lp.mean.chain[-1],type = "l",main = "Estimated MCMC chain of p with metropolis-hasting", xlab = "iteration",ylab = "p" )
Msubtitle <- "P=0.5, initial p=5"
mtext(side = 3, line = 0.4, Msubtitle)
P_1
P_1 <- multi_chain(times=1,N_sim=N_sim,burn_in=burn_in,
                   X=X, y=y_p1,initial_theda=MLE_model$coefficients,
                   true_theta=true_theta,Lp=5,M_iter=200,range=c(0.1,5),step=0.05)

P_2 <- multi_chain(times=1,N_sim=N_sim,burn_in=burn_in,
                   X=X, y=y_p2,initial_theda=MLE_model$coefficients,
                   true_theta=true_theta,Lp=5,M_iter=200,range=c(0.1,5),step=0.05)
P_2

P_3 <- multi_chain(times=5,N_sim=N_sim,burn_in=burn_in,
                   X=X, y=y_p3,initial_theda=MLE_model$coefficients,
                   true_theta=true_theta,Lp=5,M_iter=200,range=c(0.1,5),step=0.05)

P_4 <- multi_chain(times=5,N_sim=N_sim,burn_in=burn_in,
                   X=X, y=y_p4,initial_theda=MLE_model$coefficients,
                   true_theta=true_theta,Lp=5,M_iter=200,range=c(0.1,5),step=0.05)

P_5 <- multi_chain(times=5,N_sim=N_sim,burn_in=burn_in,
                   X=X, y=y_p5,initial_theda=MLE_model$coefficients,
                   true_theta=true_theta,Lp=5,M_iter=200,range=c(0.1,5),step=0.05)
P_8 <- multi_chain(times=5,N_sim=N_sim,burn_in=burn_in,
                   X=X, y=y_p8,initial_theda=MLE_model$coefficients,
                   true_theta=true_theta,Lp=5,M_iter=200,range=c(0.1,5),step=0.05)


frac <-seq(0.0001,0.005,0.0001)
frac
sketch_size <- D^2


######probit model
for (i in 1:length(frac)) {
  for (j in 1:10) {
  print(paste("Model_Frac:",i,"chain:",j,sep = ""))
  sample_index<- Compute_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=2)
  Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                       X=X[sample_index,], y=y.probit[sample_index],initial_theda=MLE_model$coefficients,
                       true_theta=true_theta,Lp=2)
  assign(paste("Probit_FP",frac[i],"chain_",j,sep = ""), Model)
  }
}

for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Uniform_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=2)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y.probit[sample_index],initial_theda=MLE_model$coefficients,
                             true_theta=true_theta,Lp=2)
    assign(paste("Probit_FP_uniform",frac[i],"chain_",j,sep = ""), Model)
  }
}

for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Lp_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=2)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y.probit[sample_index],initial_theda=MLE_model$coefficients,
                             true_theta=true_theta,Lp=2)
    assign(paste("Probit_FP_rootl2",frac[i],"chain_",j,sep = ""), Model)
  }
}



probit_fp_pos_llk_median <- c()
probit_fp_pos_llk_max <- c()
probit_fp_pos_llk_min <- c()
probit_fp_pos_llk_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- sum(llk(X,get(paste("Probit_FP",frac[i],"chain_",j,sep = ""))$post.mean,2))/sum(llk(X, P_probit$Beta.mean,2))
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  probit_fp_pos_llk_median <- c(probit_fp_pos_llk_median,median_result)
  probit_fp_pos_llk_max <- c(probit_fp_pos_llk_max,maxresult)
  probit_fp_pos_llk_min <- c(probit_fp_pos_llk_min,minresult)
  probit_fp_pos_llk_average <- c(probit_fp_pos_llk_average,meanresult)
  #assign(paste("probit_fp_pos_mean_diff_norm_",frac[i],"median",sep = ""), result)
  #assign(paste("probit_fp_pos_mean_diff_norm_",frac[i],"max",sep = ""), maxresult)
  #assign(paste("probit_fp_pos_mean_diff_norm_",frac[i],"min",sep = ""), minresult)
}

probit_fp_pos_mean_diff_norm_median <- c()
probit_fp_pos_mean_diff_norm_max <- c()
probit_fp_pos_mean_diff_norm_min <- c()
probit_fp_pos_mean_diff_norm_average <- c()


for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("Probit_FP",frac[i],"chain_",j,sep = ""))$post.mean-P_probit$Beta.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  probit_fp_pos_mean_diff_norm_median <- c(probit_fp_pos_mean_diff_norm_median,median_result)
  probit_fp_pos_mean_diff_norm_max <- c(probit_fp_pos_mean_diff_norm_max,maxresult)
  probit_fp_pos_mean_diff_norm_min <- c(probit_fp_pos_mean_diff_norm_min,minresult)
  probit_fp_pos_mean_diff_norm_average <- c(probit_fp_pos_mean_diff_norm_average,meanresult)
  
  #assign(paste("probit_fp_pos_mean_diff_norm_",frac[i],"median",sep = ""), result)
  #assign(paste("probit_fp_pos_mean_diff_norm_",frac[i],"max",sep = ""), maxresult)
  #assign(paste("probit_fp_pos_mean_diff_norm_",frac[i],"min",sep = ""), minresult)
  
}

probit_fp_uniform_pos_llk_median <- c()
probit_fp_uniform_pos_llk_max <- c()
probit_fp_uniform_pos_llk_min <- c()
probit_fp_uniform_pos_llk_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- sum(llk(X,get(paste("Probit_FP_uniform",frac[i],"chain_",j,sep = ""))$post.mean,2))/sum(llk(X, P_probit$Beta.mean,2))
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  probit_fp_uniform_pos_llk_median <- c(probit_fp_uniform_pos_llk_median,median_result)
  probit_fp_uniform_pos_llk_max <- c(probit_fp_uniform_pos_llk_max,maxresult)
  probit_fp_uniform_pos_llk_min <- c(probit_fp_uniform_pos_llk_min,minresult)
  probit_fp_uniform_pos_llk_average <- c(probit_fp_uniform_pos_llk_average,meanresult)
  #assign(paste("probit_fp_pos_mean_diff_norm_",frac[i],"median",sep = ""), result)
  #assign(paste("probit_fp_pos_mean_diff_norm_",frac[i],"max",sep = ""), maxresult)
  #assign(paste("probit_fp_pos_mean_diff_norm_",frac[i],"min",sep = ""), minresult)
}

probit_fp_uniform_pos_mean_diff_norm_median <- c()
probit_fp_uniform_pos_mean_diff_norm_max <- c()
probit_fp_uniform_pos_mean_diff_norm_min <- c()
probit_fp_uniform_pos_mean_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("Probit_FP_uniform",frac[i],"chain_",j,sep = ""))$post.mean-P_probit$Beta.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  probit_fp_uniform_pos_mean_diff_norm_median <- c(probit_fp_uniform_pos_mean_diff_norm_median,median_result)
  probit_fp_uniform_pos_mean_diff_norm_max <- c(probit_fp_uniform_pos_mean_diff_norm_max,maxresult)
  probit_fp_uniform_pos_mean_diff_norm_min <- c(probit_fp_uniform_pos_mean_diff_norm_min,minresult)
  probit_fp_uniform_pos_mean_diff_norm_average <- c(probit_fp_uniform_pos_mean_diff_norm_average,meanresult)
  
  #assign(paste("probit_fp_uniform_pos_mean_diff_norm_",frac[i],"median",sep = ""), result)
  #assign(paste("probit_fp_uniform_pos_mean_diff_norm_",frac[i],"max",sep = ""), maxresult)
  #assign(paste("probit_fp_uniform_pos_mean_diff_norm_",frac[i],"min",sep = ""), minresult)
  
}

probit_fp_rootl2_pos_llk_median <- c()
probit_fp_rootl2_pos_llk_max <- c()
probit_fp_rootl2_pos_llk_min <- c()
probit_fp_rootl2_pos_llk_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- sum(llk(X,get(paste("Probit_FP_rootl2",frac[i],"chain_",j,sep = ""))$post.mean,2))/sum(llk(X, P_probit$Beta.mean,2))
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  probit_fp_rootl2_pos_llk_median <- c(probit_fp_rootl2_pos_llk_median,median_result)
  probit_fp_rootl2_pos_llk_max <- c(probit_fp_rootl2_pos_llk_max,maxresult)
  probit_fp_rootl2_pos_llk_min <- c(probit_fp_rootl2_pos_llk_min,minresult)
  probit_fp_rootl2_pos_llk_average <- c(probit_fp_rootl2_pos_llk_average,meanresult)
  #assign(paste("probit_fp_pos_mean_diff_norm_",frac[i],"median",sep = ""), result)
  #assign(paste("probit_fp_pos_mean_diff_norm_",frac[i],"max",sep = ""), maxresult)
  #assign(paste("probit_fp_pos_mean_diff_norm_",frac[i],"min",sep = ""), minresult)
}
probit_fp_rootl2_pos_mean_diff_norm_median <- c()
probit_fp_rootl2_pos_mean_diff_norm_max <- c()
probit_fp_rootl2_pos_mean_diff_norm_min <- c()
probit_fp_rootl2_pos_mean_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("Probit_FP_rootl2",frac[i],"chain_",j,sep = ""))$post.mean-P_probit$Beta.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  probit_fp_rootl2_pos_mean_diff_norm_median <- c(probit_fp_rootl2_pos_mean_diff_norm_median,median_result)
  probit_fp_rootl2_pos_mean_diff_norm_max <- c(probit_fp_rootl2_pos_mean_diff_norm_max,maxresult)
  probit_fp_rootl2_pos_mean_diff_norm_min <- c(probit_fp_rootl2_pos_mean_diff_norm_min,minresult)
  probit_fp_rootl2_pos_mean_diff_norm_average <- c(probit_fp_rootl2_pos_mean_diff_norm_average,meanresult)
  
  #assign(paste("probit_fp_uniform_pos_mean_diff_norm_",frac[i],"median",sep = ""), result)
  #assign(paste("probit_fp_uniform_pos_mean_diff_norm_",frac[i],"max",sep = ""), maxresult)
  #assign(paste("probit_fp_uniform_pos_mean_diff_norm_",frac[i],"min",sep = ""), minresult)
  
}


p2_meandiff_probit_data <- data.frame(mean=c(probit_fp_pos_mean_diff_norm_average,
                                                probit_fp_uniform_pos_mean_diff_norm_average,
                                             probit_fp_rootl2_pos_mean_diff_norm_average),
                                         median=c(probit_fp_pos_mean_diff_norm_median,
                                                  probit_fp_uniform_pos_mean_diff_norm_median,
                                                  probit_fp_rootl2_pos_mean_diff_norm_median),
                                         max=c(probit_fp_pos_mean_diff_norm_max,
                                               probit_fp_uniform_pos_mean_diff_norm_max,
                                               probit_fp_rootl2_pos_mean_diff_norm_max),
                                         min=c(probit_fp_pos_mean_diff_norm_min,
                                               probit_fp_uniform_pos_mean_diff_norm_min,
                                               probit_fp_rootl2_pos_mean_diff_norm_min),
                                         label=c(rep("2-probit",length(frac)),rep("uniform",length(frac)),rep("rootl2",length(frac))),
                                         frac=c(rep(frac,3)))



p2_llk_probit_data <- data.frame(mean=c(probit_fp_pos_llk_average,
                                        probit_fp_uniform_pos_llk_average,
                                        probit_fp_rootl2_pos_llk_average),
                                      median=c(probit_fp_pos_llk_median,
                                               probit_fp_uniform_pos_llk_median,
                                               probit_fp_rootl2_pos_llk_median),
                                      max=c(probit_fp_pos_llk_max,
                                            probit_fp_uniform_pos_llk_max,
                                            probit_fp_rootl2_pos_llk_max),
                                      min=c(probit_fp_pos_llk_min,
                                            probit_fp_uniform_pos_llk_min,
                                            probit_fp_rootl2_pos_llk_min),
                                      label=c(rep("2-probit",length(frac)),rep("uniform",length(frac)),rep("rootl2",length(frac))),
                                      frac=c(rep(frac,3)))
p2_llk_probit_data
ggplot(p2_llk_probit_data, aes(frac, mean, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior mean", colour="method", fill="method") +
  theme_bw()




probit_fp_uniform_pos_cov_diff_norm_median <- c()
probit_fp_uniform_pos_cov_diff_norm_max <- c()
probit_fp_uniform_pos_cov_diff_norm_min <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    b <- norm(cov(get(paste("Probit_FP_uniform",frac[i],"chain_",j,sep = ""))$chain)-
                cov(P_probit$beta.mean.chain),type = "2")
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  probit_fp_uniform_pos_cov_diff_norm_median <- c(probit_fp_uniform_pos_cov_diff_norm_median,median_result)
  probit_fp_uniform_pos_cov_diff_norm_max <- c(probit_fp_uniform_pos_cov_diff_norm_max,maxresult)
  probit_fp_uniform_pos_cov_diff_norm_min <- c(probit_fp_uniform_pos_cov_diff_norm_min,minresult)
  
}

probit_fp_pos_cov_diff_norm_median <- c()
probit_fp_pos_cov_diff_norm_max <- c()
probit_fp_pos_cov_diff_norm_min <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    b <- norm(cov(get(paste("Probit_FP",frac[i],"chain_",j,sep = ""))$chain)-
                cov(P_probit$beta.mean.chain),type = "2")
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  probit_fp_pos_cov_diff_norm_median <- c(probit_fp_pos_cov_diff_norm_median,median_result)
  probit_fp_pos_cov_diff_norm_max <- c(probit_fp_pos_cov_diff_norm_max,maxresult)
  probit_fp_pos_cov_diff_norm_min <- c(probit_fp_pos_cov_diff_norm_min,minresult)
  
}

probit_fp_rootl2_cov_diff_norm_median <- c()
probit_fp_rootl2_cov_diff_norm_max <- c()
probit_fp_rootl2_cov_diff_norm_min <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    b <- norm(cov(get(paste("Probit_FP_rootl2",frac[i],"chain_",j,sep = ""))$chain)-
                cov(P_probit$beta.mean.chain),type = "2")
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  probit_fp_rootl2_cov_diff_norm_median <- c(probit_fp_rootl2_cov_diff_norm_median,median_result)
  probit_fp_rootl2_cov_diff_norm_max <- c(probit_fp_rootl2_cov_diff_norm_max,maxresult)
  probit_fp_rootl2_cov_diff_norm_min <- c(probit_fp_rootl2_cov_diff_norm_min,minresult)
  
}

p2_covdiff_covertype_data <- data.frame(
  median=c(probit_fp_pos_cov_diff_norm_median,
           probit_fp_uniform_pos_cov_diff_norm_median,
           probit_fp_rootl2_cov_diff_norm_median),
  max=c(probit_fp_pos_cov_diff_norm_max,
        probit_fp_uniform_pos_cov_diff_norm_max,
        probit_fp_rootl2_cov_diff_norm_max),
  min=c(probit_fp_pos_cov_diff_norm_min,
        probit_fp_uniform_pos_cov_diff_norm_min,
        probit_fp_rootl2_cov_diff_norm_min),
  label=c(rep("2-probit",length(frac)),rep("uniform",length(frac)),rep("rootl2",length(frac))),
  frac=c(rep(frac,3)))
p2_covdiff_covertype_data
ggplot(p2_covdiff_covertype_data, aes(frac, median, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior covariance", colour="method", fill="method") +
  theme_bw()

probit_fp_uniform_mmd_median <- c()
probit_fp_uniform_mmd_max <- c()
probit_fp_uniform_mmd_min <- c()
for (i in 1:length(frac)){
  print(paste("MMD_Frac:",i,sep = ""))
  medianchain <- c()
  for (j in 1:5) {
    b <- mmd(sample_x=tail(get(paste("Probit_FP_uniform",frac[i],"chain_",j,sep = ""))$chain,100),
             sample_y = tail(P_probit$beta.mean.chain,100))
    medianchain<- rbind(medianchain,b)
  }
  mean_result <- apply(medianchain, MARGIN = 2, mean)
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  probit_fp_uniform_mmd_median <- c(probit_fp_uniform_mmd_median,median_result)
  probit_fp_uniform_mmd_max <- c(probit_fp_uniform_mmd_max,maxresult)
  probit_fp_uniform_mmd_min <- c(probit_fp_uniform_mmd_min,maxresult)
  
}



p2_meandiff_covertype_data <- data.frame(mean=c(covertype_p_2_frac_mean_diff_norm_average,
                                                covertype_p_2_lpcoreset_frac_mean_diff_norm_average,
                                                covertype_p_2_uniform_frac_mean_diff_norm_average),
                                         median=c(covertype_p_2_frac_mean_diff_norm_median,
                                                  covertype_p_2_lpcoreset_frac_mean_diff_norm_median,
                                                  covertype_p_2_uniform_frac_mean_diff_norm_median),
                                         max=c(covertype_p_2_frac_mean_diff_norm_max,
                                               covertype_p_2_lpcoreset_frac_mean_diff_norm_max,
                                               covertype_p_2_uniform_frac_mean_diff_norm_max),
                                         min=c(covertype_p_2_frac_mean_diff_norm_min,
                                               covertype_p_2_lpcoreset_frac_mean_diff_norm_min,
                                               covertype_p_2_uniform_frac_mean_diff_norm_min),
                                         label=c(rep("2-probit",length(frac)),rep("root l2",length(frac)),rep("uniform",length(frac))),
                                         frac=c(rep(frac,3)))
p2_meandiff_covertype_data
ggplot(p2_meandiff_covertype_data, aes(frac, mean, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior mean", colour="method", fill="method") +
  theme_bw()
probit_fp_pos_cov_diff_norm_median <- c()
probit_fp_pos_cov_diff_norm_max <- c()
probit_fp_pos_cov_diff_norm_min <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("Probit_FP",frac[i],"chain_",j,sep = ""))$chain)-
                cov(P_probit$beta.mean.chain),type = "2")
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  probit_fp_pos_cov_diff_norm_median <- c(probit_fp_pos_cov_diff_norm_median,median_result)
  probit_fp_pos_cov_diff_norm_max <- c(probit_fp_pos_cov_diff_norm_max,maxresult)
  probit_fp_pos_cov_diff_norm_min <- c(probit_fp_pos_cov_diff_norm_min,minresult)
  
}

plot(frac,probit_fp_pos_cov_diff_norm_median,type = "l",xlab = "percentage",ylab = "Norm", main="L2 norm difference of posterior covariance")
Msubtitle <- "P=Probit"
mtext(side = 3, line = 0.4, Msubtitle)


probit_fp_mmd_median <- c()
probit_fp_mmd_max <- c()
probit_fp_mmd_min <- c()
for (i in 1:length(frac)){
  print(paste("MMD_Frac:",i,sep = ""))
  medianchain <- c()
  for (j in 1:5) {
    b <- mmd(sample_x=tail(get(paste("Probit_FP",frac[i],"chain_",j,sep = ""))$chain,100),
             sample_y = tail(P_probit$beta.mean.chain,100))
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  probit_fp_mmd_median <- c(probit_fp_mmd_median,median_result)
  probit_fp_mmd_max <- c(probit_fp_mmd_max,maxresult)
  probit_fp_mmd_min <- c(probit_fp_mmd_min,maxresult)
  
}

plot(frac,probit_mmd_FP,type = "l",xlab = "percentage",ylab = "MMD", main="Maximum Mean Discrepancy")
Msubtitle <- "P=Probit"
mtext(side = 3, line = 0.4, Msubtitle)
probit_time_fp <- c()

for (i in 1:length(frac)){
  print(paste("MMD_Frac:",i,sep = ""))
  
  b <- get(paste("Probit_FP",frac[i],"chain_1",sep = ""))$modeling.time[3]
  probit_time_fp<- rbind(probit_time_fp,b)
}
probit_time_fp
plot(frac,probit_time_fp,type = "l",xlab = "percentage",ylab = "Time(seconds)", main="Modeling time for different sample")
Msubtitle <- "P=Probit"
mtext(side = 3, line = 0.4, Msubtitle)

#logit model


for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Compute_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=1)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y.logit[sample_index],initial_theda=MLE_model$coefficients,
                             true_theta=true_theta,Lp=1)
    assign(paste("Logit_FP",frac[i],"chain_",j,sep = ""), Model)
  }
}

for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Uniform_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=1)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y.logit[sample_index],initial_theda=MLE_model$coefficients,
                             true_theta=true_theta,Lp=1)
    assign(paste("Logit_FP_uniform",frac[i],"chain_",j,sep = ""), Model)
  }
}

for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Lp_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=1)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y.logit[sample_index],initial_theda=MLE_model$coefficients,
                             true_theta=true_theta,Lp=1)
    assign(paste("Logit_FP_rootl2",frac[i],"chain_",j,sep = ""), Model)
  }
}


Logit_fp_pos_mean_diff_norm_median <- c()
Logit_fp_pos_mean_diff_norm_max <- c()
Logit_fp_pos_mean_diff_norm_min <- c()
Logit_fp_pos_mean_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("Logit_FP",frac[i],"chain_",j,sep = ""))$post.mean-P_logit$Beta.mean,type = "2")
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  Logit_fp_pos_mean_diff_norm_median <- c(Logit_fp_pos_mean_diff_norm_median,median_result)
  Logit_fp_pos_mean_diff_norm_max <- c(Logit_fp_pos_mean_diff_norm_max,maxresult)
  Logit_fp_pos_mean_diff_norm_min <- c(Logit_fp_pos_mean_diff_norm_min,minresult)
  Logit_fp_pos_mean_diff_norm_average <- c(Logit_fp_pos_mean_diff_norm_average,meanresult)
  
  #assign(paste("Logit_fp_pos_mean_diff_norm_",frac[i],"median",sep = ""), result)
  #assign(paste("Logit_fp_pos_mean_diff_norm_",frac[i],"max",sep = ""), maxresult)
  #assign(paste("Logit_fp_pos_mean_diff_norm_",frac[i],"min",sep = ""), minresult)
  
}

Logit_fp_uniform_pos_mean_diff_norm_median <- c()
Logit_fp_uniform_pos_mean_diff_norm_max <- c()
Logit_fp_uniform_pos_mean_diff_norm_min <- c()
Logit_fp_uniform_pos_mean_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("Logit_FP_uniform",frac[i],"chain_",j,sep = ""))$post.mean-P_logit$Beta.mean,type = "2")
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  Logit_fp_uniform_pos_mean_diff_norm_median <- c(Logit_fp_uniform_pos_mean_diff_norm_median,median_result)
  Logit_fp_uniform_pos_mean_diff_norm_max <- c(Logit_fp_uniform_pos_mean_diff_norm_max,maxresult)
  Logit_fp_uniform_pos_mean_diff_norm_min <- c(Logit_fp_uniform_pos_mean_diff_norm_min,minresult)
  Logit_fp_uniform_pos_mean_diff_norm_average <- c(Logit_fp_uniform_pos_mean_diff_norm_average,meanresult)
  
  #assign(paste("Logit_fp_uniform_pos_mean_diff_norm_",frac[i],"median",sep = ""), result)
  #assign(paste("Logit_fp_uniform_pos_mean_diff_norm_",frac[i],"max",sep = ""), maxresult)
  #assign(paste("Logit_fp_uniform_pos_mean_diff_norm_",frac[i],"min",sep = ""), minresult)
  
}

Logit_fp_rootl2_pos_mean_diff_norm_median <- c()
Logit_fp_rootl2_pos_mean_diff_norm_max <- c()
Logit_fp_rootl2_pos_mean_diff_norm_min <- c()
Logit_fp_rootl2_pos_mean_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("Logit_FP_rootl2",frac[i],"chain_",j,sep = ""))$post.mean-P_logit$Beta.mean,type = "2")
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  Logit_fp_rootl2_pos_mean_diff_norm_median <- c(Logit_fp_rootl2_pos_mean_diff_norm_median,median_result)
  Logit_fp_rootl2_pos_mean_diff_norm_max <- c(Logit_fp_rootl2_pos_mean_diff_norm_max,maxresult)
  Logit_fp_rootl2_pos_mean_diff_norm_min <- c(Logit_fp_rootl2_pos_mean_diff_norm_min,minresult)
  Logit_fp_rootl2_pos_mean_diff_norm_average <- c(Logit_fp_rootl2_pos_mean_diff_norm_average,meanresult)
  
  #assign(paste("Logit_fp_uniform_pos_mean_diff_norm_",frac[i],"median",sep = ""), result)
  #assign(paste("Logit_fp_uniform_pos_mean_diff_norm_",frac[i],"max",sep = ""), maxresult)
  #assign(paste("Logit_fp_uniform_pos_mean_diff_norm_",frac[i],"min",sep = ""), minresult)
  
}


p2_meandiff_Logit_data <- data.frame(mean=c(Logit_fp_pos_mean_diff_norm_average,
                                             Logit_fp_uniform_pos_mean_diff_norm_average,
                                             Logit_fp_rootl2_pos_mean_diff_norm_average),
                                      median=c(Logit_fp_pos_mean_diff_norm_median,
                                               Logit_fp_uniform_pos_mean_diff_norm_median,
                                               Logit_fp_rootl2_pos_mean_diff_norm_median),
                                      max=c(Logit_fp_pos_mean_diff_norm_max,
                                            Logit_fp_uniform_pos_mean_diff_norm_max,
                                            Logit_fp_rootl2_pos_mean_diff_norm_max),
                                      min=c(Logit_fp_pos_mean_diff_norm_min,
                                            Logit_fp_pos_mean_diff_norm_min,
                                            Logit_fp_rootl2_pos_mean_diff_norm_min),
                                      label=c(rep("2-probit",length(frac)),rep("uniform",length(frac)),rep("rootl2",length(frac))),
                                      frac=c(rep(frac,3)))
p2_meandiff_Logit_data
ggplot(p2_meandiff_Logit_data, aes(frac, mean, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior mean", colour="method", fill="method") +
  theme_bw()

Logit_fp_uniform_pos_cov_diff_norm_median <- c()
Logit_fp_uniform_pos_cov_diff_norm_max <- c()
Logit_fp_uniform_pos_cov_diff_norm_min <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("Logit_FP_uniform",frac[i],"chain_",j,sep = ""))$chain)-
                cov(P_logit$beta.mean.chain),type = "2")
    
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  Logit_fp_uniform_pos_cov_diff_norm_median <- c(Logit_fp_uniform_pos_cov_diff_norm_median,median_result)
  Logit_fp_uniform_pos_cov_diff_norm_max <- c(Logit_fp_uniform_pos_cov_diff_norm_max,maxresult)
  Logit_fp_uniform_pos_cov_diff_norm_min <- c(Logit_fp_uniform_pos_cov_diff_norm_min,minresult)
  
}

Logit_fp_pos_cov_diff_norm_median <- c()
Logit_fp_pos_cov_diff_norm_max <- c()
Logit_fp_pos_cov_diff_norm_min <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("Logit_FP",frac[i],"chain_",j,sep = ""))$chain)-
                cov(P_logit$beta.mean.chain),type = "2")
    
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  Logit_fp_pos_cov_diff_norm_median <- c(Logit_fp_pos_cov_diff_norm_median,median_result)
  Logit_fp_pos_cov_diff_norm_max <- c(Logit_fp_pos_cov_diff_norm_max,maxresult)
  Logit_fp_pos_cov_diff_norm_min <- c(Logit_fp_pos_cov_diff_norm_min,minresult)
  
}

Logit_fp_rootl2_cov_diff_norm_median <- c()
Logit_fp_rootl2_cov_diff_norm_max <- c()
Logit_fp_rootl2_cov_diff_norm_min <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("Logit_FP_rootl2",frac[i],"chain_",j,sep = ""))$chain)-
                cov(P_logit$beta.mean.chain),type = "2")
    
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  Logit_fp_rootl2_cov_diff_norm_median <- c(Logit_fp_rootl2_cov_diff_norm_median,median_result)
  Logit_fp_rootl2_cov_diff_norm_max <- c(Logit_fp_rootl2_cov_diff_norm_max,maxresult)
  Logit_fp_rootl2_cov_diff_norm_min <- c(Logit_fp_rootl2_cov_diff_norm_min,minresult)
  
}

p2_covdiff_covertype_data <- data.frame(
  median=c(Logit_fp_pos_cov_diff_norm_median,
           Logit_fp_uniform_pos_cov_diff_norm_median,
           Logit_fp_rootl2_cov_diff_norm_median),
  max=c(Logit_fp_pos_cov_diff_norm_max,
        Logit_fp_uniform_pos_cov_diff_norm_max,
        Logit_fp_rootl2_cov_diff_norm_max),
  min=c(Logit_fp_pos_cov_diff_norm_min,
        Logit_fp_uniform_pos_cov_diff_norm_min,
        Logit_fp_rootl2_cov_diff_norm_min),
  label=c(rep("2-probit",length(frac)),rep("uniform",length(frac)),rep("rootl2",length(frac))),
  frac=c(rep(frac,3)))
p2_covdiff_covertype_data
ggplot(p2_covdiff_covertype_data, aes(frac, median, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior covariance", colour="method", fill="method") +
  theme_bw()

Logit_fp_mmd_median <- c()
Logit_fp_mmd_max <- c()
Logit_fp_mmd_min <- c()
for (i in 1:length(frac)){
  print(paste("MMD_Frac:",i,sep = ""))
  medianchain <- c()
  for (j in 1:5) {
    b <- mmd(sample_x=tail(get(paste("Logit_FP",frac[i],"chain_",j,sep = ""))$chain,100),
             sample_y = tail(P_logit$beta.mean.chain,100))
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  Logit_fp_mmd_median <- c(Logit_fp_mmd_median,median_result)
  Logit_fp_mmd_max <- c(Logit_fp_mmd_max,maxresult)
  Logit_fp_mmd_min <- c(Logit_fp_mmd_min,minresult)
  
}

plot(frac,Logit_fp_mmd_median,type = "l",xlab = "percentage",ylab = "MMD", main="Maximum Mean Discrepancy")
Msubtitle <- "P=Logit"
mtext(side = 3, line = 0.4, Msubtitle)
Logit_time_fp <- c()

for (i in 1:length(frac)){
  print(paste("MMD_Frac:",i,sep = ""))
  
  b <- get(paste("Logit_FP",frac[i],"chain_1",sep = ""))$modeling.time[3]
  Logit_time_fp<- rbind(Logit_time_fp,b)
}
Logit_time_fp
plot(frac,Logit_time_fp,type = "l",xlab = "percentage",ylab = "Time(seconds)", main="Modeling time for different sample")
Msubtitle <- "P=Logit"
mtext(side = 3, line = 0.4, Msubtitle)


logit_mmd_FP
logit_pos_mean_diff_norm_FP
logit_pos_cov_diff_norm_FP
logit_time_fp
plot(frac,logit_pos_mean_diff_norm_FP,type = "l",xlab = "percentage",ylab = "Norm", main="L2 norm difference of posterior mean")
Msubtitle <- "P=logit"
mtext(side = 3, line = 0.4, Msubtitle)
plot(frac,logit_pos_cov_diff_norm_FP,type = "l",xlab = "percentage",ylab = "Norm", main="L2 norm difference of posterior covariance")
Msubtitle <- "P=logit"
mtext(side = 3, line = 0.4, Msubtitle)
plot(frac,Logit_mmd_FP_median,type = "l",xlab = "percentage",ylab = "MMD", main="Maximum Mean Discrepancy")
Msubtitle <- "P=logit"
mtext(side = 3, line = 0.4, Msubtitle)
plot(frac,logit_time_fp,type = "l",xlab = "percentage",ylab = "Time(seconds)", main="Modeling time for different sample")
Msubtitle <- "P=logit"
mtext(side = 3, line = 0.4, Msubtitle)

######################P=0.5#####################

for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Compute_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=0.5)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y_p0.5[sample_index],initial_theda=MLE_model$coefficients,
                             true_theta=true_theta,Lp=0.5)
    assign(paste("P_0.5_FP",frac[i],"chain_",j,sep = ""), Model)
  }
}

for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Uniform_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=0.5)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y_p0.5[sample_index],initial_theda=MLE_model$coefficients,
                             true_theta=true_theta,Lp=0.5)
    assign(paste("P_0.5_FP_uniform",frac[i],"chain_",j,sep = ""), Model)
  }
}

for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Lp_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=0.5)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y_p0.5[sample_index],initial_theda=MLE_model$coefficients,
                             true_theta=true_theta,Lp=0.5)
    assign(paste("P_0.5_FP_rootl2",frac[i],"chain_",j,sep = ""), Model)
  }
}

P_0.5_fp_pos_mean_diff_norm_median <- c()
P_0.5_fp_pos_mean_diff_norm_max <- c()
P_0.5_fp_pos_mean_diff_norm_min <- c()
P_0.5_fp_pos_mean_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("P_0.5_FP",frac[i],"chain_",j,sep = ""))$post.mean-P_0.5$Beta.mean,type = "2")
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  P_0.5_fp_pos_mean_diff_norm_median <- c(P_0.5_fp_pos_mean_diff_norm_median,median_result)
  P_0.5_fp_pos_mean_diff_norm_max <- c(P_0.5_fp_pos_mean_diff_norm_max,maxresult)
  P_0.5_fp_pos_mean_diff_norm_min <- c(P_0.5_fp_pos_mean_diff_norm_min,minresult)
  P_0.5_fp_pos_mean_diff_norm_average <- c(P_0.5_fp_pos_mean_diff_norm_average,meanresult)
  
  #assign(paste("P_0.5_fp_pos_mean_diff_norm_",frac[i],"median",sep = ""), result)
  #assign(paste("P_0.5_fp_pos_mean_diff_norm_",frac[i],"max",sep = ""), maxresult)
  #assign(paste("P_0.5_fp_pos_mean_diff_norm_",frac[i],"min",sep = ""), minresult)
  
}

P_0.5_fp_uniform_pos_mean_diff_norm_median <- c()
P_0.5_fp_uniform_pos_mean_diff_norm_max <- c()
P_0.5_fp_uniform_pos_mean_diff_norm_min <- c()
P_0.5_fp_uniform_pos_mean_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("P_0.5_FP_uniform",frac[i],"chain_",j,sep = ""))$post.mean-P_0.5$Beta.mean,type = "2")
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  P_0.5_fp_uniform_pos_mean_diff_norm_median <- c(P_0.5_fp_uniform_pos_mean_diff_norm_median,median_result)
  P_0.5_fp_uniform_pos_mean_diff_norm_max <- c(P_0.5_fp_uniform_pos_mean_diff_norm_max,maxresult)
  P_0.5_fp_uniform_pos_mean_diff_norm_min <- c(P_0.5_fp_uniform_pos_mean_diff_norm_min,minresult)
  P_0.5_fp_uniform_pos_mean_diff_norm_average <- c(P_0.5_fp_uniform_pos_mean_diff_norm_average,meanresult)
  
  #assign(paste("P_0.5_fp_uniform_pos_mean_diff_norm_",frac[i],"median",sep = ""), result)
  #assign(paste("P_0.5_fp_uniform_pos_mean_diff_norm_",frac[i],"max",sep = ""), maxresult)
  #assign(paste("P_0.5_fp_uniform_pos_mean_diff_norm_",frac[i],"min",sep = ""), minresult)
  
}

P_0.5_fp_rootl2_pos_mean_diff_norm_median <- c()
P_0.5_fp_rootl2_pos_mean_diff_norm_max <- c()
P_0.5_fp_rootl2_pos_mean_diff_norm_min <- c()
P_0.5_fp_rootl2_pos_mean_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("P_0.5_FP_rootl2",frac[i],"chain_",j,sep = ""))$post.mean-P_0.5$Beta.mean,type = "2")
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  P_0.5_fp_rootl2_pos_mean_diff_norm_median <- c(P_0.5_fp_rootl2_pos_mean_diff_norm_median,median_result)
  P_0.5_fp_rootl2_pos_mean_diff_norm_max <- c(P_0.5_fp_rootl2_pos_mean_diff_norm_max,maxresult)
  P_0.5_fp_rootl2_pos_mean_diff_norm_min <- c(P_0.5_fp_rootl2_pos_mean_diff_norm_min,minresult)
  P_0.5_fp_rootl2_pos_mean_diff_norm_average <- c(P_0.5_fp_rootl2_pos_mean_diff_norm_average,meanresult)
  
  #assign(paste("P_0.5_fp_uniform_pos_mean_diff_norm_",frac[i],"median",sep = ""), result)
  #assign(paste("P_0.5_fp_uniform_pos_mean_diff_norm_",frac[i],"max",sep = ""), maxresult)
  #assign(paste("P_0.5_fp_uniform_pos_mean_diff_norm_",frac[i],"min",sep = ""), minresult)
  
}


P_0.5_meandiff <- data.frame(mean=c(P_0.5_fp_pos_mean_diff_norm_average,
                                            P_0.5_fp_uniform_pos_mean_diff_norm_average,
                                            P_0.5_fp_rootl2_pos_mean_diff_norm_average),
                                     median=c(P_0.5_fp_pos_mean_diff_norm_median,
                                              P_0.5_fp_uniform_pos_mean_diff_norm_median,
                                              P_0.5_fp_rootl2_pos_mean_diff_norm_median),
                                     max=c(P_0.5_fp_pos_mean_diff_norm_max,
                                           P_0.5_fp_uniform_pos_mean_diff_norm_max,
                                           P_0.5_fp_rootl2_pos_mean_diff_norm_max),
                                     min=c(P_0.5_fp_pos_mean_diff_norm_min,
                                           P_0.5_fp_pos_mean_diff_norm_min,
                                           P_0.5_fp_rootl2_pos_mean_diff_norm_min),
                                     label=c(rep("2-probit",length(frac)),rep("uniform",length(frac)),rep("rootl2",length(frac))),
                                     frac=c(rep(frac,3)))
P_0.5_meandiff
ggplot(P_0.5_meandiff, aes(frac, mean, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior mean (P=0.5)", colour="method", fill="method") +
  theme_bw()

P_0.5_fp_uniform_pos_cov_diff_norm_median <- c()
P_0.5_fp_uniform_pos_cov_diff_norm_max <- c()
P_0.5_fp_uniform_pos_cov_diff_norm_min <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("P_0.5_FP_uniform",frac[i],"chain_",j,sep = ""))$chain)-
                cov(P_0.5$beta.mean.chain),type = "2")
    
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  P_0.5_fp_uniform_pos_cov_diff_norm_median <- c(P_0.5_fp_uniform_pos_cov_diff_norm_median,median_result)
  P_0.5_fp_uniform_pos_cov_diff_norm_max <- c(P_0.5_fp_uniform_pos_cov_diff_norm_max,maxresult)
  P_0.5_fp_uniform_pos_cov_diff_norm_min <- c(P_0.5_fp_uniform_pos_cov_diff_norm_min,minresult)
  
}

P_0.5_fp_pos_cov_diff_norm_median <- c()
P_0.5_fp_pos_cov_diff_norm_max <- c()
P_0.5_fp_pos_cov_diff_norm_min <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("P_0.5_FP",frac[i],"chain_",j,sep = ""))$chain)-
                cov(P_0.5$beta.mean.chain),type = "2")
    
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  P_0.5_fp_pos_cov_diff_norm_median <- c(P_0.5_fp_pos_cov_diff_norm_median,median_result)
  P_0.5_fp_pos_cov_diff_norm_max <- c(P_0.5_fp_pos_cov_diff_norm_max,maxresult)
  P_0.5_fp_pos_cov_diff_norm_min <- c(P_0.5_fp_pos_cov_diff_norm_min,minresult)
  
}

P_0.5_fp_rootl2_cov_diff_norm_median <- c()
P_0.5_fp_rootl2_cov_diff_norm_max <- c()
P_0.5_fp_rootl2_cov_diff_norm_min <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("P_0.5_FP_rootl2",frac[i],"chain_",j,sep = ""))$chain)-
                cov(P_0.5$beta.mean.chain),type = "2")
    
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  P_0.5_fp_rootl2_cov_diff_norm_median <- c(P_0.5_fp_rootl2_cov_diff_norm_median,median_result)
  P_0.5_fp_rootl2_cov_diff_norm_max <- c(P_0.5_fp_rootl2_cov_diff_norm_max,maxresult)
  P_0.5_fp_rootl2_cov_diff_norm_min <- c(P_0.5_fp_rootl2_cov_diff_norm_min,minresult)
  
}

P_0.5_covdiff <- data.frame(
  median=c(P_0.5_fp_pos_cov_diff_norm_median,
           P_0.5_fp_uniform_pos_cov_diff_norm_median,
           P_0.5_fp_rootl2_cov_diff_norm_median),
  max=c(P_0.5_fp_pos_cov_diff_norm_max,
        P_0.5_fp_uniform_pos_cov_diff_norm_max,
        P_0.5_fp_rootl2_cov_diff_norm_max),
  min=c(P_0.5_fp_pos_cov_diff_norm_min,
        P_0.5_fp_uniform_pos_cov_diff_norm_min,
        P_0.5_fp_rootl2_cov_diff_norm_min),
  label=c(rep("2-probit",length(frac)),rep("uniform",length(frac)),rep("rootl2",length(frac))),
  frac=c(rep(frac,3)))
ggplot(P_0.5_covdiff, aes(frac, median, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior covariance (P=0.5)", colour="method", fill="method") +
  theme_bw()



######################P=1#####################

for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Compute_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=1)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y_p1[sample_index],initial_theda=MLE_model$coefficients,
                             true_theta=true_theta,Lp=1)
    assign(paste("P_1_FP",frac[i],"chain_",j,sep = ""), Model)
  }
}

for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Uniform_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=1)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y_p1[sample_index],initial_theda=MLE_model$coefficients,
                             true_theta=true_theta,Lp=1)
    assign(paste("P_1_FP_uniform",frac[i],"chain_",j,sep = ""), Model)
  }
}

for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Lp_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=1)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y_p1[sample_index],initial_theda=MLE_model$coefficients,
                             true_theta=true_theta,Lp=1)
    assign(paste("P_1_FP_rootl2",frac[i],"chain_",j,sep = ""), Model)
  }
}

P_1_fp_pos_mean_diff_norm_median <- c()
P_1_fp_pos_mean_diff_norm_max <- c()
P_1_fp_pos_mean_diff_norm_min <- c()
P_1_fp_pos_mean_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("P_1_FP",frac[i],"chain_",j,sep = ""))$post.mean-p_1$Beta.mean,type = "2")
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  P_1_fp_pos_mean_diff_norm_median <- c(P_1_fp_pos_mean_diff_norm_median,median_result)
  P_1_fp_pos_mean_diff_norm_max <- c(P_1_fp_pos_mean_diff_norm_max,maxresult)
  P_1_fp_pos_mean_diff_norm_min <- c(P_1_fp_pos_mean_diff_norm_min,minresult)
  P_1_fp_pos_mean_diff_norm_average <- c(P_1_fp_pos_mean_diff_norm_average,meanresult)
  
  #assign(paste("P_1_fp_pos_mean_diff_norm_",frac[i],"median",sep = ""), result)
  #assign(paste("P_1_fp_pos_mean_diff_norm_",frac[i],"max",sep = ""), maxresult)
  #assign(paste("P_1_fp_pos_mean_diff_norm_",frac[i],"min",sep = ""), minresult)
  
}

P_1_fp_uniform_pos_mean_diff_norm_median <- c()
P_1_fp_uniform_pos_mean_diff_norm_max <- c()
P_1_fp_uniform_pos_mean_diff_norm_min <- c()
P_1_fp_uniform_pos_mean_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("P_1_FP_uniform",frac[i],"chain_",j,sep = ""))$post.mean-P_1$Beta.mean,type = "2")
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  P_1_fp_uniform_pos_mean_diff_norm_median <- c(P_1_fp_uniform_pos_mean_diff_norm_median,median_result)
  P_1_fp_uniform_pos_mean_diff_norm_max <- c(P_1_fp_uniform_pos_mean_diff_norm_max,maxresult)
  P_1_fp_uniform_pos_mean_diff_norm_min <- c(P_1_fp_uniform_pos_mean_diff_norm_min,minresult)
  P_1_fp_uniform_pos_mean_diff_norm_average <- c(P_1_fp_uniform_pos_mean_diff_norm_average,meanresult)
  
  #assign(paste("P_1_fp_uniform_pos_mean_diff_norm_",frac[i],"median",sep = ""), result)
  #assign(paste("P_1_fp_uniform_pos_mean_diff_norm_",frac[i],"max",sep = ""), maxresult)
  #assign(paste("P_1_fp_uniform_pos_mean_diff_norm_",frac[i],"min",sep = ""), minresult)
  
}

P_1_fp_rootl2_pos_mean_diff_norm_median <- c()
P_1_fp_rootl2_pos_mean_diff_norm_max <- c()
P_1_fp_rootl2_pos_mean_diff_norm_min <- c()
P_1_fp_rootl2_pos_mean_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("P_1_FP_rootl2",frac[i],"chain_",j,sep = ""))$post.mean-P_1$Beta.mean,type = "2")
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  P_1_fp_rootl2_pos_mean_diff_norm_median <- c(P_1_fp_rootl2_pos_mean_diff_norm_median,median_result)
  P_1_fp_rootl2_pos_mean_diff_norm_max <- c(P_1_fp_rootl2_pos_mean_diff_norm_max,maxresult)
  P_1_fp_rootl2_pos_mean_diff_norm_min <- c(P_1_fp_rootl2_pos_mean_diff_norm_min,minresult)
  P_1_fp_rootl2_pos_mean_diff_norm_average <- c(P_1_fp_rootl2_pos_mean_diff_norm_average,meanresult)
  
  #assign(paste("P_1_fp_uniform_pos_mean_diff_norm_",frac[i],"median",sep = ""), result)
  #assign(paste("P_1_fp_uniform_pos_mean_diff_norm_",frac[i],"max",sep = ""), maxresult)
  #assign(paste("P_1_fp_uniform_pos_mean_diff_norm_",frac[i],"min",sep = ""), minresult)
  
}


P_1_meandiff <- data.frame(mean=c(P_1_fp_pos_mean_diff_norm_average,
                                    P_1_fp_uniform_pos_mean_diff_norm_average,
                                    P_1_fp_rootl2_pos_mean_diff_norm_average),
                             median=c(P_1_fp_pos_mean_diff_norm_median,
                                      P_1_fp_uniform_pos_mean_diff_norm_median,
                                      P_1_fp_rootl2_pos_mean_diff_norm_median),
                             max=c(P_1_fp_pos_mean_diff_norm_max,
                                   P_1_fp_uniform_pos_mean_diff_norm_max,
                                   P_1_fp_rootl2_pos_mean_diff_norm_max),
                             min=c(P_1_fp_pos_mean_diff_norm_min,
                                   P_1_fp_pos_mean_diff_norm_min,
                                   P_1_fp_rootl2_pos_mean_diff_norm_min),
                             label=c(rep("2-probit",length(frac)),rep("uniform",length(frac)),rep("rootl2",length(frac))),
                             frac=c(rep(frac,3)))
P_1_meandiff
ggplot(P_1_meandiff, aes(frac, mean, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior mean (P=1)", colour="method", fill="method") +
  theme_bw()

P_1_fp_uniform_pos_cov_diff_norm_median <- c()
P_1_fp_uniform_pos_cov_diff_norm_max <- c()
P_1_fp_uniform_pos_cov_diff_norm_min <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("P_1_FP_uniform",frac[i],"chain_",j,sep = ""))$chain)-
                cov(P_1$beta.mean.chain),type = "2")
    
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  P_1_fp_uniform_pos_cov_diff_norm_median <- c(P_1_fp_uniform_pos_cov_diff_norm_median,median_result)
  P_1_fp_uniform_pos_cov_diff_norm_max <- c(P_1_fp_uniform_pos_cov_diff_norm_max,maxresult)
  P_1_fp_uniform_pos_cov_diff_norm_min <- c(P_1_fp_uniform_pos_cov_diff_norm_min,minresult)
  
}

P_1_fp_pos_cov_diff_norm_median <- c()
P_1_fp_pos_cov_diff_norm_max <- c()
P_1_fp_pos_cov_diff_norm_min <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("P_1_FP",frac[i],"chain_",j,sep = ""))$chain)-
                cov(P_1$beta.mean.chain),type = "2")
    
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  P_1_fp_pos_cov_diff_norm_median <- c(P_1_fp_pos_cov_diff_norm_median,median_result)
  P_1_fp_pos_cov_diff_norm_max <- c(P_1_fp_pos_cov_diff_norm_max,maxresult)
  P_1_fp_pos_cov_diff_norm_min <- c(P_1_fp_pos_cov_diff_norm_min,minresult)
  
}

P_1_fp_rootl2_cov_diff_norm_median <- c()
P_1_fp_rootl2_cov_diff_norm_max <- c()
P_1_fp_rootl2_cov_diff_norm_min <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("P_1_FP_rootl2",frac[i],"chain_",j,sep = ""))$chain)-
                cov(P_1$beta.mean.chain),type = "2")
    
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  P_1_fp_rootl2_cov_diff_norm_median <- c(P_1_fp_rootl2_cov_diff_norm_median,median_result)
  P_1_fp_rootl2_cov_diff_norm_max <- c(P_1_fp_rootl2_cov_diff_norm_max,maxresult)
  P_1_fp_rootl2_cov_diff_norm_min <- c(P_1_fp_rootl2_cov_diff_norm_min,minresult)
  
}

P_1_covdiff <- data.frame(
  median=c(P_1_fp_pos_cov_diff_norm_median,
           P_1_fp_uniform_pos_cov_diff_norm_median,
           P_1_fp_rootl2_cov_diff_norm_median),
  max=c(P_1_fp_pos_cov_diff_norm_max,
        P_1_fp_uniform_pos_cov_diff_norm_max,
        P_1_fp_rootl2_cov_diff_norm_max),
  min=c(P_1_fp_pos_cov_diff_norm_min,
        P_1_fp_uniform_pos_cov_diff_norm_min,
        P_1_fp_rootl2_cov_diff_norm_min),
  label=c(rep("2-probit",length(frac)),rep("uniform",length(frac)),rep("rootl2",length(frac))),
  frac=c(rep(frac,3)))
ggplot(P_1_covdiff, aes(frac, median, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior covariance (P=1)", colour="method", fill="method") +
  theme_bw()






plot(P_probit$Lp.mean.chain[-1],type = "l",main = "Estimated MCMC chain of p with metropolis-hasting", xlab = "iteration",ylab = "p" )
Msubtitle <- "P=probit, initial p = 5"
mtext(side = 3, line = 0.4, Msubtitle)
lines(1:2000,rep(2,2000),col=2)
legend("topright",legend=c("p=2"), col=c(2), 
       bty = 'n', lwd = 2, inset = c(-0.15, 0), lty = 1, cex = 0.73)

par(mfrow = c(2, 2))

plot(P_1.5$Lp.mean.chain[-1],type = "l",main = "Estimated MCMC chain of p with metropolis-hasting", xlab = "iteration",ylab = "p" )
Msubtitle <- "P=1.5, initial p = 5"
mtext(side = 3, line = 0.4, Msubtitle)
lines(1:2000,rep(1.5,2000),col=2)
legend("topright",legend=c("p=1.5"), col=c(2), 
       bty = 'n', lwd = 2, inset = c(-0.4, 0), lty = 1, cex = 0.73)


plot(P_1$Lp.mean.chain[-1],type = "l",main = "Estimated MCMC chain of p with metropolis-hasting", xlab = "iteration",ylab = "p" )
Msubtitle <- "P=1, initial p = 5"
mtext(side = 3, line = 0.4, Msubtitle)
lines(1:2000,rep(1,2000),col=2)
legend("topright",legend=c("p=1"), col=c(2), 
       bty = 'n', lwd = 2, inset = c(-0.3, 0), lty = 1, cex = 0.73)


plot(P_1.5$Lp.mean.chain[-1],type = "l",main = "Estimated MCMC chain of p with metropolis-hasting", xlab = "iteration",ylab = "p" )
Msubtitle <- "P=1.5, initial p = 5"
mtext(side = 3, line = 0.4, Msubtitle)
lines(1:2000,rep(1.5,2000),col=2)
legend("topright",legend=c("p=1.5"), col=c(2), 
       bty = 'n', lwd = 2, inset = c(-0.4, 0), lty = 1, cex = 0.73)


plot(P_2$Lp.mean.chain[-1],type = "l",main = "Estimated MCMC chain of p with metropolis-hasting", xlab = "iteration",ylab = "p" )
Msubtitle <- "P=2, initial p = 5"
mtext(side = 3, line = 0.4, Msubtitle)
lines(1:2000,rep(2,2000),col=2)
legend("topright",legend=c("p=2"), col=c(2), 
       bty = 'n', lwd = 2, inset = c(-0.3, 0), lty = 1, cex = 0.73)

plot(P_3$Lp.mean.chain[-1],type = "l",main = "Estimated MCMC chain of p with metropolis-hasting", xlab = "iteration",ylab = "p" )
Msubtitle <- "P=3, initial p = 5"
mtext(side = 3, line = 0.4, Msubtitle)
lines(1:2000,rep(3,2000),col=2)
legend("topright",legend=c("p=3"), col=c(2), 
       bty = 'n', lwd = 2, inset = c(-0.33, 0), lty = 1, cex = 0.73)



plot(P_4$Lp.mean.chain[-1],type = "l",main = "Estimated MCMC chain of p with metropolis-hasting", xlab = "iteration",ylab = "p" )
Msubtitle <- "P=4, initial p = 5"
mtext(side = 3, line = 0.4, Msubtitle)
lines(1:2000,rep(4,2000),col=2)
legend("topright",legend=c("p=4"), col=c(2), 
       bty = 'n', lwd = 2, inset = c(-0.33, 0), lty = 1, cex = 0.73)


plot(P_5$Lp.mean.chain[-1],type = "l",main = "Estimated MCMC chain of p with metropolis-hasting",ylim = c(3.5,7), xlab = "iteration",ylab = "p" )
Msubtitle <- "P=5, initial p = 5"
mtext(side = 3, line = 0.4, Msubtitle)
lines(1:2000,rep(5,2000),col=2)
legend("topright",legend=c("p=5"), col=c(2), 
       bty = 'n', lwd = 2, inset = c(-0.33, 0), lty = 1, cex = 0.73)


plot(P_8$Lp.mean.chain[-1],type = "l",main = "Estimated MCMC chain of p with metropolis-hasting",ylim = c(3.5,8.5), xlab = "iteration",ylab = "p" )
Msubtitle <- "P=8, initial p = 5"
mtext(side = 3, line = 0.4, Msubtitle)
lines(1:2000,rep(8,2000),col=2)
legend("topright",legend=c("p=8"), col=c(2), 
       bty = 'n', lwd = 2, inset = c(-0.33, 0), lty = 1, cex = 0.73)

P_logit$Lp.mean
var(P_logit$Lp.mean.chain)
norm(as.matrix(P_logit$Beta.mean-true_theta,2))
pi_logit_predict <- exp(X%*%P_logit$Beta.mean)/(1+exp(X%*%P_logit$Beta.mean))
sum((y.logit-pi_logit_predict)^2)

P_probit$Lp.mean
var(P_probit$Lp.mean.chain)
norm(as.matrix(P_probit$Beta.mean-true_theta,2))
pi_probit_predict <- pgnorm(X%*%P_probit$Beta.mean,alpha = Phi_p(P_probit$Lp.mean),beta = P_probit$Lp.mean)
sum((y.probit-pi_probit_predict)^2)

P_1.5$Lp.mean
var(P_1.5$Lp.mean.chain)
norm(as.matrix(P_0.5$Beta.mean-true_theta,2))
pi_probit_predict <- pgnorm(X%*%P_0.5$Beta.mean,alpha = Phi_p(P_0.5$Lp.mean),beta = P_0.5$Lp.mean)
sum((y.probit-pi_probit_predict)^2)

P_1$Lp.mean
var(P_1$Lp.mean.chain)
norm(as.matrix(P_1$Beta.mean-true_theta,2))
pi_probit_predict <- pgnorm(X%*%P_1$Beta.mean,alpha = Phi_p(P_1$Lp.mean),beta = P_1$Lp.mean)
sum((y.probit-pi_probit_predict)^2)

P_1.5$Lp.mean
var(P_1.5$Lp.mean.chain)
norm(as.matrix(P_1.5$Beta.mean-true_theta,2))
pi_probit_predict <- pgnorm(X%*%P_1.5$Beta.mean,alpha = Phi_p(P_1.5$Lp.mean),beta = P_1.5$Lp.mean)
sum((y.probit-pi_probit_predict)^2)

P_2$Lp.mean
var(P_2$Lp.mean.chain)
norm(as.matrix(P_2$Beta.mean-true_theta,2))
pi_probit_predict <- pgnorm(X%*%P_2$Beta.mean,alpha = Phi_p(P_2$Lp.mean),beta = P_2$Lp.mean)
sum((y.probit-pi_probit_predict)^2)

P_3$Lp.mean
var(P_3$Lp.mean.chain)
norm(as.matrix(P_3$Beta.mean-true_theta,2))
pi_probit_predict <- pgnorm(X%*%P_3$Beta.mean,alpha = Phi_p(P_3$Lp.mean),beta = P_3$Lp.mean)
sum((y.probit-pi_probit_predict)^2)

P_4$Lp.mean
var(P_4$Lp.mean.chain)
norm(as.matrix(P_4$Beta.mean-true_theta,2))
pi_probit_predict <- pgnorm(X%*%P_4$Beta.mean,alpha = Phi_p(P_4$Lp.mean),beta = P_4$Lp.mean)
sum((y.probit-pi_probit_predict)^2)

P_5$Lp.mean
var(P_5$Lp.mean.chain)
norm(as.matrix(P_5$Beta.mean-true_theta,2))
pi_probit_predict <- pgnorm(X%*%P_5$Beta.mean,alpha = Phi_p(P_5$Lp.mean),beta = P_5$Lp.mean)
sum((y.probit-pi_probit_predict)^2)

P_8$Lp.mean
var(P_8$Lp.mean.chain)
norm(as.matrix(P_8$Beta.mean-true_theta,2))
pi_probit_predict <- pgnorm(X%*%P_8$Beta.mean,alpha = Phi_p(P_8$Lp.mean),beta = P_8$Lp.mean)
sum((y.probit-pi_probit_predict)^2)

plot(density(tail(P_8$Lp.mean.chain,1000)))

plot(Xb, dgnorm(Xb, mu = 0, alpha = sqrt(8), beta = 8),
     xlab = "x", ylab = expression(p(x)))
Xb <- subset(Xb,Xb<10)
Xb <- subset(Xb,-10<Xb)

density0.5 <- data.frame(Xb=Xb,density=dgnorm(Xb,alpha = Phi_p(0.5),beta = 0.5))
density1 <- data.frame(Xb=Xb,density=dgnorm(Xb,alpha = Phi_p(1),beta = 1))
density1.5 <- data.frame(Xb=Xb,density=dgnorm(Xb,alpha = Phi_p(1.5),beta = 1.5))
density2 <- data.frame(Xb=Xb,density=dgnorm(Xb,alpha = Phi_p(2),beta = 2))
density2.5 <- data.frame(Xb=Xb,density=dgnorm(Xb,alpha = Phi_p(2.5),beta = 2.5))
density3 <- data.frame(Xb=Xb,density=dgnorm(Xb,alpha = Phi_p(3),beta = 3))
density3.5 <- data.frame(Xb=Xb,density=dgnorm(Xb,alpha = Phi_p(3.5),beta = 3.5))
density4 <- data.frame(Xb=Xb,density=dgnorm(Xb,alpha = Phi_p(4),beta = 4))
density5 <- data.frame(Xb=Xb,density=dgnorm(Xb,alpha = Phi_p(5),beta = 5))
density8<- data.frame(Xb=Xb,density=dgnorm(Xb,alpha = Phi_p(8),beta = 8))
density0.5$p <- "0.5"
density1$p <- "1"
density1.5$p <- "1.5"
density2$p <- "2"
density2.5$p <- "2.5"
density3$p <- "3"
density3.5$p <- "3.5"
density4$p <- "4"
density5$p <- "5"
density8$p <- "8"

probit_density_plot <- rbind(density0.5,density1,density1.5,density2,density2.5,density3,density3.5,
                             density4,density5,density8)
ggplot(probit_density_plot, aes(x=Xb,y=density, colour = p)) + geom_line()+
  ggtitle("Plot of XB versus density for different P") +
  xlab("XB") + ylab("density")

prob0.5 <- data.frame(Xb=Xb,pi=pgnorm(Xb,alpha = Phi_p(0.5),beta = 0.5))
prob1 <- data.frame(Xb=Xb,pi=pgnorm(Xb,alpha = Phi_p(1),beta = 1))
prob1.5 <- data.frame(Xb=Xb,pi=pgnorm(Xb,alpha = Phi_p(1.5),beta = 1.5))
prob2 <- data.frame(Xb=Xb,pi=pgnorm(Xb,alpha = Phi_p(2),beta = 2))
prob2.5 <- data.frame(Xb=Xb,pi=pgnorm(Xb,alpha = Phi_p(2.5),beta = 2.5))
prob3 <- data.frame(Xb=Xb,pi=pgnorm(Xb,alpha = Phi_p(3),beta = 3))
prob3.5 <- data.frame(Xb=Xb,pi=pgnorm(Xb,alpha = Phi_p(3.5),beta = 3.5))
prob4 <- data.frame(Xb=Xb,pi=pgnorm(Xb,alpha = Phi_p(4),beta = 4))
prob5 <- data.frame(Xb=Xb,pi=pgnorm(Xb,alpha = Phi_p(5),beta = 5))
prob8<- data.frame(Xb=Xb,pi=pgnorm(Xb,alpha = Phi_p(8),beta = 8))
prob0.5$p <- "0.5"
prob1$p <- "1"
prob1.5$p <- "1.5"
prob2$p <- "2"
prob2.5$p <- "2.5"
prob3$p <- "3"
prob3.5$p <- "3.5"
prob4$p <- "4"
prob5$p <- "5"
prob8$p <- "8"

probit_pi_plot <- rbind(prob0.5,prob1,prob1.5,prob2,prob2.5,prob3,prob3.5,
                        prob4,prob5,prob8)
ggplot(probit_pi_plot, aes(x=Xb,y=pi, colour = p)) + geom_line()+ggtitle("Plot of XB versus binomial probability for different P") +
  xlab("XB") + ylab("Probability")
#geom_text(aes(x=true_theta[2], label="true theta", y=20),colour="blue", angle=90, vjust = 1.2)


df <- data.frame(Xb=Xb,aaa=dgnorm(Xb,alpha = Phi_p(2),beta = 2))
df2 <- data.frame(Xb=Xb,aaa=pgnorm(Xb,alpha = Phi_p(2),beta = 2))


ggplot(data=df, aes(x=Xb,y=aaa)) +  geom_line()

plot(Xb,pi.probit,type = "l")

probit_probit <- data.frame(data=tail(P_probit$Lp.mean.chain,1000))
probit_logit <- data.frame(data=tail(P_probit$Lp.mean.chain,1000))
probit0.5 <- data.frame(data=tail(P_0.5$Lp.mean.chain,1000))
probit1 <- data.frame(data=tail(P_1$Lp.mean.chain,1000))
probit1.5 <- data.frame(data=tail(P_1.5$Lp.mean.chain,1000))
probit2 <- data.frame(data=tail(P_2$Lp.mean.chain,1000))
probit2.5 <- data.frame(data=tail(P_2.5$Lp.mean.chain,1000))
probit3 <- data.frame(data=tail(P_3$Lp.mean.chain,1000))
probit3.5 <- data.frame(data=tail(P_3.5$Lp.mean.chain,1000))
probit4 <- data.frame(data=tail(P_4$Lp.mean.chain,1000))
probit5 <- data.frame(data=tail(P_5$Lp.mean.chain,1000))
probit8 <- data.frame(data=tail(P_8$Lp.mean.chain,1000))
probit_probit$p <- "probit"
probit_probit$p <- "logit"
probit0.5$p <- "0.5"
probit1$p <- "1"
probit1.5$p <- "1.5"
probit2$p <- "2"
probit2.5$p <- "2.5"
probit3$p <- "3"
probit3.5$p <- "3.5"
probit4$p <- "4"
probit5$p <- "5"
probit8$p <- "8"
probit_probit
probit5

posterior_plot <- rbind(probit0.5, 
                        probit1,probit1.5,probit2,probit2.5,
                        probit3,probit3.5,probit4,probit5,probit8)
posterior_plot$x <- rep(1:1000,10)
Lp_use <- data.frame(x=c(0.5,1,1.5,2,2.5,3,3.5,4,5,8),p=c(0.5,1,1.5,2,2.5,3,3.5,4,5,8))
Lp_use
ggplot(posterior_plot, aes(data, fill = p)) + geom_density(alpha = 1)+
  geom_segment(aes(x=2,xend=2,y=0,yend=25),colour="blue",size=0.35,linetype=2)+
  geom_segment(aes(x=0.5,xend=0.5,y=0,yend=25),colour="blue",size=0.35,linetype=2)+
  geom_segment(aes(x=1,xend=1,y=0,yend=25),colour="blue",size=0.35,linetype=2)+
  geom_segment(aes(x=1.5,xend=1.5,y=0,yend=25),colour="blue",size=0.35,linetype=2)+
  geom_segment(aes(x=2,xend=2,y=0,yend=25),colour="blue",size=0.35,linetype=2)+
  geom_segment(aes(x=2.5,xend=2.5,y=0,yend=25),colour="blue",size=0.35,linetype=2)+
  geom_segment(aes(x=3,xend=3,y=0,yend=25),colour="blue",size=0.35,linetype=2)+
  geom_segment(aes(x=3.5,xend=3.5,y=0,yend=25),colour="blue",size=0.35,linetype=2)+
  geom_segment(aes(x=4,xend=4,y=0,yend=25),colour="blue",size=0.35,linetype=2)+
  geom_segment(aes(x=5,xend=5,y=0,yend=25),colour="blue",size=0.35,linetype=2)+
  geom_segment(aes(x=8,xend=8,y=0,yend=25),colour="blue",size=0.35,linetype=2)+
  ggtitle("Posterior Distribution of the MCMC chain of different P") +
  xlab("XB") + ylab("density")+
  geom_text(size=2.5,aes(x=0.5, label="p=0.5", y=-0.5), 
            colour="blue", angle=0, vjust = 1.2)+
  geom_text(size=2.5,aes(x=1, label="p=1", y=-0.5), 
            colour="blue", angle=0, vjust = 1.2)+
  geom_text(size=2.5,aes(x=1.5, label="p=1,5", y=-0.5), 
            colour="blue", angle=0, vjust = 1.2)+
  geom_text(size=2.5,aes(x=2, label="p=2", y=-0.5), 
            colour="blue", angle=0, vjust = 1.2)+
  geom_text(size=2.5,aes(x=2.5, label="p=2.5", y=-0.5), 
            colour="blue", angle=0, vjust = 1.2)+
  geom_text(size=2.5,aes(x=3, label="p=3", y=-0.5), 
            colour="blue", angle=0, vjust = 1.2)+
  geom_text(size=2.5,aes(x=3.5, label="p=3.5", y=-0.5), 
            colour="blue", angle=0, vjust = 1.2)+
  geom_text(size=2.5,aes(x=4, label="p=4", y=-0.5), 
            colour="blue", angle=0, vjust = 1.2)+
  geom_text(size=2.5,aes(x=5, label="p=5", y=-0.5), 
            colour="blue", angle=0, vjust = 1.2)+
  geom_text(size=2.5,aes(x=8, label="p=8", y=-0.5), 
            colour="blue", angle=0, vjust = 1.2)

ggplot(posterior_plot, aes(x=x,y=data, colour = p)) + geom_line()+
  scale_y_continuous(breaks = seq(0,8,1)) +
  ggtitle("MCMC estimated chain of different P") +
  xlab("Iteration") + ylab("P value")


geom_line(data = data.frame(x = c(2, Inf), y = c(2, 2)), aes(x = x , y = y))
#geom_line(data=Lp_use, aes(x=Lp,y=p),
#           linetype="dashed")
#geom_vline(xintercept=true_theta[2],colour="red")+
geom_text(aes(x=true_theta[2], label="true theta", y=20), 
          colour="blue", angle=90, vjust = 1.2)
Lp
probit_p <- data.frame(data=tail(P_probit$Lp.mean.chain,1000))

ggplot(probit_p,aes(data)) + geom_density(alpha = 0.5)
#+geom_vline(xintercept=true_theta[2],colour="red")+geom_text(aes(x=true_theta[2], label="true theta", y=20), colour="blue", angle=90, vjust = 1.2)





######################P=2#####################

for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Compute_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=2)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y_p2[sample_index],initial_theda=MLE_model$coefficients,
                             true_theta=true_theta,Lp=2)
    assign(paste("P_2_FP",frac[i],"chain_",j,sep = ""), Model)
  }
}

for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Uniform_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=2)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y_p2[sample_index],initial_theda=MLE_model$coefficients,
                             true_theta=true_theta,Lp=2)
    assign(paste("P_2_FP_uniform",frac[i],"chain_",j,sep = ""), Model)
  }
}

for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Lp_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=2)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y_p2[sample_index],initial_theda=MLE_model$coefficients,
                             true_theta=true_theta,Lp=2)
    assign(paste("P_2_FP_rootl2",frac[i],"chain_",j,sep = ""), Model)
  }
}

P_2_fp_pos_mean_diff_norm_median <- c()
P_2_fp_pos_mean_diff_norm_max <- c()
P_2_fp_pos_mean_diff_norm_min <- c()
P_2_fp_pos_mean_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("P_2_FP",frac[i],"chain_",j,sep = ""))$post.mean-P_2$Beta.mean,type = "2")
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  P_2_fp_pos_mean_diff_norm_median <- c(P_2_fp_pos_mean_diff_norm_median,median_result)
  P_2_fp_pos_mean_diff_norm_max <- c(P_2_fp_pos_mean_diff_norm_max,maxresult)
  P_2_fp_pos_mean_diff_norm_min <- c(P_2_fp_pos_mean_diff_norm_min,minresult)
  P_2_fp_pos_mean_diff_norm_average <- c(P_2_fp_pos_mean_diff_norm_average,meanresult)
  
  #assign(paste("P_2_fp_pos_mean_diff_norm_",frac[i],"median",sep = ""), result)
  #assign(paste("P_2_fp_pos_mean_diff_norm_",frac[i],"max",sep = ""), maxresult)
  #assign(paste("P_2_fp_pos_mean_diff_norm_",frac[i],"min",sep = ""), minresult)
  
}

P_2_fp_uniform_pos_mean_diff_norm_median <- c()
P_2_fp_uniform_pos_mean_diff_norm_max <- c()
P_2_fp_uniform_pos_mean_diff_norm_min <- c()
P_2_fp_uniform_pos_mean_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("P_2_FP_uniform",frac[i],"chain_",j,sep = ""))$post.mean-P_2$Beta.mean,type = "2")
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  P_2_fp_uniform_pos_mean_diff_norm_median <- c(P_2_fp_uniform_pos_mean_diff_norm_median,median_result)
  P_2_fp_uniform_pos_mean_diff_norm_max <- c(P_2_fp_uniform_pos_mean_diff_norm_max,maxresult)
  P_2_fp_uniform_pos_mean_diff_norm_min <- c(P_2_fp_uniform_pos_mean_diff_norm_min,minresult)
  P_2_fp_uniform_pos_mean_diff_norm_average <- c(P_2_fp_uniform_pos_mean_diff_norm_average,meanresult)
  
  #assign(paste("P_2_fp_uniform_pos_mean_diff_norm_",frac[i],"median",sep = ""), result)
  #assign(paste("P_2_fp_uniform_pos_mean_diff_norm_",frac[i],"max",sep = ""), maxresult)
  #assign(paste("P_2_fp_uniform_pos_mean_diff_norm_",frac[i],"min",sep = ""), minresult)
  
}

P_2_fp_rootl2_pos_mean_diff_norm_median <- c()
P_2_fp_rootl2_pos_mean_diff_norm_max <- c()
P_2_fp_rootl2_pos_mean_diff_norm_min <- c()
P_2_fp_rootl2_pos_mean_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("P_2_FP_rootl2",frac[i],"chain_",j,sep = ""))$post.mean-P_2$Beta.mean,type = "2")
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  P_2_fp_rootl2_pos_mean_diff_norm_median <- c(P_2_fp_rootl2_pos_mean_diff_norm_median,median_result)
  P_2_fp_rootl2_pos_mean_diff_norm_max <- c(P_2_fp_rootl2_pos_mean_diff_norm_max,maxresult)
  P_2_fp_rootl2_pos_mean_diff_norm_min <- c(P_2_fp_rootl2_pos_mean_diff_norm_min,minresult)
  P_2_fp_rootl2_pos_mean_diff_norm_average <- c(P_2_fp_rootl2_pos_mean_diff_norm_average,meanresult)
  
  #assign(paste("P_2_fp_uniform_pos_mean_diff_norm_",frac[i],"median",sep = ""), result)
  #assign(paste("P_2_fp_uniform_pos_mean_diff_norm_",frac[i],"max",sep = ""), maxresult)
  #assign(paste("P_2_fp_uniform_pos_mean_diff_norm_",frac[i],"min",sep = ""), minresult)
  
}


P_2_meandiff <- data.frame(mean=c(P_2_fp_pos_mean_diff_norm_average,
                                    P_2_fp_uniform_pos_mean_diff_norm_average,
                                    P_2_fp_rootl2_pos_mean_diff_norm_average),
                             median=c(P_2_fp_pos_mean_diff_norm_median,
                                      P_2_fp_uniform_pos_mean_diff_norm_median,
                                      P_2_fp_rootl2_pos_mean_diff_norm_median),
                             max=c(P_2_fp_pos_mean_diff_norm_max,
                                   P_2_fp_uniform_pos_mean_diff_norm_max,
                                   P_2_fp_rootl2_pos_mean_diff_norm_max),
                             min=c(P_2_fp_pos_mean_diff_norm_min,
                                   P_2_fp_pos_mean_diff_norm_min,
                                   P_2_fp_rootl2_pos_mean_diff_norm_min),
                             label=c(rep("2-probit",length(frac)),rep("uniform",length(frac)),rep("rootl2",length(frac))),
                             frac=c(rep(frac,3)))
P_2_meandiff
ggplot(P_2_meandiff, aes(frac, mean, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior mean (P=2)", colour="method", fill="method") +
  theme_bw()

P_2_fp_uniform_pos_cov_diff_norm_median <- c()
P_2_fp_uniform_pos_cov_diff_norm_max <- c()
P_2_fp_uniform_pos_cov_diff_norm_min <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("P_2_FP_uniform",frac[i],"chain_",j,sep = ""))$chain)-
                cov(P_2$beta.mean.chain),type = "2")
    
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  P_2_fp_uniform_pos_cov_diff_norm_median <- c(P_2_fp_uniform_pos_cov_diff_norm_median,median_result)
  P_2_fp_uniform_pos_cov_diff_norm_max <- c(P_2_fp_uniform_pos_cov_diff_norm_max,maxresult)
  P_2_fp_uniform_pos_cov_diff_norm_min <- c(P_2_fp_uniform_pos_cov_diff_norm_min,minresult)
  
}

P_2_fp_pos_cov_diff_norm_median <- c()
P_2_fp_pos_cov_diff_norm_max <- c()
P_2_fp_pos_cov_diff_norm_min <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("P_2_FP",frac[i],"chain_",j,sep = ""))$chain)-
                cov(P_2$beta.mean.chain),type = "2")
    
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  P_2_fp_pos_cov_diff_norm_median <- c(P_2_fp_pos_cov_diff_norm_median,median_result)
  P_2_fp_pos_cov_diff_norm_max <- c(P_2_fp_pos_cov_diff_norm_max,maxresult)
  P_2_fp_pos_cov_diff_norm_min <- c(P_2_fp_pos_cov_diff_norm_min,minresult)
  
}

P_2_fp_rootl2_cov_diff_norm_median <- c()
P_2_fp_rootl2_cov_diff_norm_max <- c()
P_2_fp_rootl2_cov_diff_norm_min <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("P_2_FP_rootl2",frac[i],"chain_",j,sep = ""))$chain)-
                cov(P_2$beta.mean.chain),type = "2")
    
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  P_2_fp_rootl2_cov_diff_norm_median <- c(P_2_fp_rootl2_cov_diff_norm_median,median_result)
  P_2_fp_rootl2_cov_diff_norm_max <- c(P_2_fp_rootl2_cov_diff_norm_max,maxresult)
  P_2_fp_rootl2_cov_diff_norm_min <- c(P_2_fp_rootl2_cov_diff_norm_min,minresult)
  
}

P_2_covdiff <- data.frame(
  median=c(P_2_fp_pos_cov_diff_norm_median,
           P_2_fp_uniform_pos_cov_diff_norm_median,
           P_2_fp_rootl2_cov_diff_norm_median),
  max=c(P_2_fp_pos_cov_diff_norm_max,
        P_2_fp_uniform_pos_cov_diff_norm_max,
        P_2_fp_rootl2_cov_diff_norm_max),
  min=c(P_2_fp_pos_cov_diff_norm_min,
        P_2_fp_uniform_pos_cov_diff_norm_min,
        P_2_fp_rootl2_cov_diff_norm_min),
  label=c(rep("2-probit",length(frac)),rep("uniform",length(frac)),rep("rootl2",length(frac))),
  frac=c(rep(frac,3)))
ggplot(P_2_covdiff, aes(frac, median, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior covariance (P=2)", colour="method", fill="method") +
  theme_bw()

######################P=3#####################

for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Compute_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=3)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y_p3[sample_index],initial_theda=MLE_model$coefficients,
                             true_theta=true_theta,Lp=3)
    assign(paste("P_3_FP",frac[i],"chain_",j,sep = ""), Model)
  }
}

for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Uniform_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=3)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y_p3[sample_index],initial_theda=MLE_model$coefficients,
                             true_theta=true_theta,Lp=3)
    assign(paste("P_3_FP_uniform",frac[i],"chain_",j,sep = ""), Model)
  }
}

for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Lp_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=3)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y_p3[sample_index],initial_theda=MLE_model$coefficients,
                             true_theta=true_theta,Lp=3)
    assign(paste("P_3_FP_rootl2",frac[i],"chain_",j,sep = ""), Model)
  }
}

P_3_fp_pos_mean_diff_norm_median <- c()
P_3_fp_pos_mean_diff_norm_max <- c()
P_3_fp_pos_mean_diff_norm_min <- c()
P_3_fp_pos_mean_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("P_3_FP",frac[i],"chain_",j,sep = ""))$post.mean-P_3$Beta.mean,type = "2")
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  P_3_fp_pos_mean_diff_norm_median <- c(P_3_fp_pos_mean_diff_norm_median,median_result)
  P_3_fp_pos_mean_diff_norm_max <- c(P_3_fp_pos_mean_diff_norm_max,maxresult)
  P_3_fp_pos_mean_diff_norm_min <- c(P_3_fp_pos_mean_diff_norm_min,minresult)
  P_3_fp_pos_mean_diff_norm_average <- c(P_3_fp_pos_mean_diff_norm_average,meanresult)
  
  #assign(paste("P_3_fp_pos_mean_diff_norm_",frac[i],"median",sep = ""), result)
  #assign(paste("P_3_fp_pos_mean_diff_norm_",frac[i],"max",sep = ""), maxresult)
  #assign(paste("P_3_fp_pos_mean_diff_norm_",frac[i],"min",sep = ""), minresult)
  
}

P_3_fp_uniform_pos_mean_diff_norm_median <- c()
P_3_fp_uniform_pos_mean_diff_norm_max <- c()
P_3_fp_uniform_pos_mean_diff_norm_min <- c()
P_3_fp_uniform_pos_mean_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("P_3_FP_uniform",frac[i],"chain_",j,sep = ""))$post.mean-P_3$Beta.mean,type = "2")
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  P_3_fp_uniform_pos_mean_diff_norm_median <- c(P_3_fp_uniform_pos_mean_diff_norm_median,median_result)
  P_3_fp_uniform_pos_mean_diff_norm_max <- c(P_3_fp_uniform_pos_mean_diff_norm_max,maxresult)
  P_3_fp_uniform_pos_mean_diff_norm_min <- c(P_3_fp_uniform_pos_mean_diff_norm_min,minresult)
  P_3_fp_uniform_pos_mean_diff_norm_average <- c(P_3_fp_uniform_pos_mean_diff_norm_average,meanresult)
  
  #assign(paste("P_3_fp_uniform_pos_mean_diff_norm_",frac[i],"median",sep = ""), result)
  #assign(paste("P_3_fp_uniform_pos_mean_diff_norm_",frac[i],"max",sep = ""), maxresult)
  #assign(paste("P_3_fp_uniform_pos_mean_diff_norm_",frac[i],"min",sep = ""), minresult)
  
}

P_3_fp_rootl2_pos_mean_diff_norm_median <- c()
P_3_fp_rootl2_pos_mean_diff_norm_max <- c()
P_3_fp_rootl2_pos_mean_diff_norm_min <- c()
P_3_fp_rootl2_pos_mean_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("P_3_FP_rootl2",frac[i],"chain_",j,sep = ""))$post.mean-P_3$Beta.mean,type = "2")
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  P_3_fp_rootl2_pos_mean_diff_norm_median <- c(P_3_fp_rootl2_pos_mean_diff_norm_median,median_result)
  P_3_fp_rootl2_pos_mean_diff_norm_max <- c(P_3_fp_rootl2_pos_mean_diff_norm_max,maxresult)
  P_3_fp_rootl2_pos_mean_diff_norm_min <- c(P_3_fp_rootl2_pos_mean_diff_norm_min,minresult)
  P_3_fp_rootl2_pos_mean_diff_norm_average <- c(P_3_fp_rootl2_pos_mean_diff_norm_average,meanresult)
  
  #assign(paste("P_3_fp_uniform_pos_mean_diff_norm_",frac[i],"median",sep = ""), result)
  #assign(paste("P_3_fp_uniform_pos_mean_diff_norm_",frac[i],"max",sep = ""), maxresult)
  #assign(paste("P_3_fp_uniform_pos_mean_diff_norm_",frac[i],"min",sep = ""), minresult)
  
}


P_3_meandiff <- data.frame(mean=c(P_3_fp_pos_mean_diff_norm_average,
                                    P_3_fp_uniform_pos_mean_diff_norm_average,
                                    P_3_fp_rootl2_pos_mean_diff_norm_average),
                             median=c(P_3_fp_pos_mean_diff_norm_median,
                                      P_3_fp_uniform_pos_mean_diff_norm_median,
                                      P_3_fp_rootl2_pos_mean_diff_norm_median),
                             max=c(P_3_fp_pos_mean_diff_norm_max,
                                   P_3_fp_uniform_pos_mean_diff_norm_max,
                                   P_3_fp_rootl2_pos_mean_diff_norm_max),
                             min=c(P_3_fp_pos_mean_diff_norm_min,
                                   P_3_fp_pos_mean_diff_norm_min,
                                   P_3_fp_rootl2_pos_mean_diff_norm_min),
                             label=c(rep("2-probit",length(frac)),rep("uniform",length(frac)),rep("rootl2",length(frac))),
                             frac=c(rep(frac,3)))
P_3_meandiff
ggplot(P_3_meandiff, aes(frac, mean, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior mean (P=3)", colour="method", fill="method") +
  theme_bw()

P_3_fp_uniform_pos_cov_diff_norm_median <- c()
P_3_fp_uniform_pos_cov_diff_norm_max <- c()
P_3_fp_uniform_pos_cov_diff_norm_min <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("P_3_FP_uniform",frac[i],"chain_",j,sep = ""))$chain)-
                cov(P_3$beta.mean.chain),type = "2")
    
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  P_3_fp_uniform_pos_cov_diff_norm_median <- c(P_3_fp_uniform_pos_cov_diff_norm_median,median_result)
  P_3_fp_uniform_pos_cov_diff_norm_max <- c(P_3_fp_uniform_pos_cov_diff_norm_max,maxresult)
  P_3_fp_uniform_pos_cov_diff_norm_min <- c(P_3_fp_uniform_pos_cov_diff_norm_min,minresult)
  
}

P_3_fp_pos_cov_diff_norm_median <- c()
P_3_fp_pos_cov_diff_norm_max <- c()
P_3_fp_pos_cov_diff_norm_min <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("P_3_FP",frac[i],"chain_",j,sep = ""))$chain)-
                cov(P_3$beta.mean.chain),type = "2")
    
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  P_3_fp_pos_cov_diff_norm_median <- c(P_3_fp_pos_cov_diff_norm_median,median_result)
  P_3_fp_pos_cov_diff_norm_max <- c(P_3_fp_pos_cov_diff_norm_max,maxresult)
  P_3_fp_pos_cov_diff_norm_min <- c(P_3_fp_pos_cov_diff_norm_min,minresult)
  
}

P_3_fp_rootl2_cov_diff_norm_median <- c()
P_3_fp_rootl2_cov_diff_norm_max <- c()
P_3_fp_rootl2_cov_diff_norm_min <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("P_3_FP_rootl2",frac[i],"chain_",j,sep = ""))$chain)-
                cov(P_3$beta.mean.chain),type = "2")
    
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  P_3_fp_rootl2_cov_diff_norm_median <- c(P_3_fp_rootl2_cov_diff_norm_median,median_result)
  P_3_fp_rootl2_cov_diff_norm_max <- c(P_3_fp_rootl2_cov_diff_norm_max,maxresult)
  P_3_fp_rootl2_cov_diff_norm_min <- c(P_3_fp_rootl2_cov_diff_norm_min,minresult)
  
}

P_3_covdiff <- data.frame(
  median=c(P_3_fp_pos_cov_diff_norm_median,
           P_3_fp_uniform_pos_cov_diff_norm_median,
           P_3_fp_rootl2_cov_diff_norm_median),
  max=c(P_3_fp_pos_cov_diff_norm_max,
        P_3_fp_uniform_pos_cov_diff_norm_max,
        P_3_fp_rootl2_cov_diff_norm_max),
  min=c(P_3_fp_pos_cov_diff_norm_min,
        P_3_fp_uniform_pos_cov_diff_norm_min,
        P_3_fp_rootl2_cov_diff_norm_min),
  label=c(rep("2-probit",length(frac)),rep("uniform",length(frac)),rep("rootl2",length(frac))),
  frac=c(rep(frac,3)))
ggplot(P_3_covdiff, aes(frac, median, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior covariance (P=3)", colour="method", fill="method") +
  theme_bw()


######################P=5#####################

for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Compute_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=5)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y_p5[sample_index],initial_theda=MLE_model$coefficients,
                             true_theta=true_theta,Lp=5)
    assign(paste("P_5_FP",frac[i],"chain_",j,sep = ""), Model)
  }
}

for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Uniform_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=5)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y_p5[sample_index],initial_theda=MLE_model$coefficients,
                             true_theta=true_theta,Lp=5)
    assign(paste("P_5_FP_uniform",frac[i],"chain_",j,sep = ""), Model)
  }
}

for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Lp_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=5)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y_p5[sample_index],initial_theda=MLE_model$coefficients,
                             true_theta=true_theta,Lp=5)
    assign(paste("P_5_FP_rootl2",frac[i],"chain_",j,sep = ""), Model)
  }
}

P_5_fp_pos_mean_diff_norm_median <- c()
P_5_fp_pos_mean_diff_norm_max <- c()
P_5_fp_pos_mean_diff_norm_min <- c()
P_5_fp_pos_mean_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("P_5_FP",frac[i],"chain_",j,sep = ""))$post.mean-P_5$Beta.mean,type = "2")
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  P_5_fp_pos_mean_diff_norm_median <- c(P_5_fp_pos_mean_diff_norm_median,median_result)
  P_5_fp_pos_mean_diff_norm_max <- c(P_5_fp_pos_mean_diff_norm_max,maxresult)
  P_5_fp_pos_mean_diff_norm_min <- c(P_5_fp_pos_mean_diff_norm_min,minresult)
  P_5_fp_pos_mean_diff_norm_average <- c(P_5_fp_pos_mean_diff_norm_average,meanresult)
  
  #assign(paste("P_5_fp_pos_mean_diff_norm_",frac[i],"median",sep = ""), result)
  #assign(paste("P_5_fp_pos_mean_diff_norm_",frac[i],"max",sep = ""), maxresult)
  #assign(paste("P_5_fp_pos_mean_diff_norm_",frac[i],"min",sep = ""), minresult)
  
}

P_5_fp_uniform_pos_mean_diff_norm_median <- c()
P_5_fp_uniform_pos_mean_diff_norm_max <- c()
P_5_fp_uniform_pos_mean_diff_norm_min <- c()
P_5_fp_uniform_pos_mean_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("P_5_FP_uniform",frac[i],"chain_",j,sep = ""))$post.mean-P_5$Beta.mean,type = "2")
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  P_5_fp_uniform_pos_mean_diff_norm_median <- c(P_5_fp_uniform_pos_mean_diff_norm_median,median_result)
  P_5_fp_uniform_pos_mean_diff_norm_max <- c(P_5_fp_uniform_pos_mean_diff_norm_max,maxresult)
  P_5_fp_uniform_pos_mean_diff_norm_min <- c(P_5_fp_uniform_pos_mean_diff_norm_min,minresult)
  P_5_fp_uniform_pos_mean_diff_norm_average <- c(P_5_fp_uniform_pos_mean_diff_norm_average,meanresult)
  
  #assign(paste("P_5_fp_uniform_pos_mean_diff_norm_",frac[i],"median",sep = ""), result)
  #assign(paste("P_5_fp_uniform_pos_mean_diff_norm_",frac[i],"max",sep = ""), maxresult)
  #assign(paste("P_5_fp_uniform_pos_mean_diff_norm_",frac[i],"min",sep = ""), minresult)
  
}

P_5_fp_rootl2_pos_mean_diff_norm_median <- c()
P_5_fp_rootl2_pos_mean_diff_norm_max <- c()
P_5_fp_rootl2_pos_mean_diff_norm_min <- c()
P_5_fp_rootl2_pos_mean_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("P_5_FP_rootl2",frac[i],"chain_",j,sep = ""))$post.mean-P_5$Beta.mean,type = "2")
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  P_5_fp_rootl2_pos_mean_diff_norm_median <- c(P_5_fp_rootl2_pos_mean_diff_norm_median,median_result)
  P_5_fp_rootl2_pos_mean_diff_norm_max <- c(P_5_fp_rootl2_pos_mean_diff_norm_max,maxresult)
  P_5_fp_rootl2_pos_mean_diff_norm_min <- c(P_5_fp_rootl2_pos_mean_diff_norm_min,minresult)
  P_5_fp_rootl2_pos_mean_diff_norm_average <- c(P_5_fp_rootl2_pos_mean_diff_norm_average,meanresult)
  
  #assign(paste("P_5_fp_uniform_pos_mean_diff_norm_",frac[i],"median",sep = ""), result)
  #assign(paste("P_5_fp_uniform_pos_mean_diff_norm_",frac[i],"max",sep = ""), maxresult)
  #assign(paste("P_5_fp_uniform_pos_mean_diff_norm_",frac[i],"min",sep = ""), minresult)
  
}


P_5_meandiff <- data.frame(mean=c(P_5_fp_pos_mean_diff_norm_average,
                                    P_5_fp_uniform_pos_mean_diff_norm_average,
                                    P_5_fp_rootl2_pos_mean_diff_norm_average),
                             median=c(P_5_fp_pos_mean_diff_norm_median,
                                      P_5_fp_uniform_pos_mean_diff_norm_median,
                                      P_5_fp_rootl2_pos_mean_diff_norm_median),
                             max=c(P_5_fp_pos_mean_diff_norm_max,
                                   P_5_fp_uniform_pos_mean_diff_norm_max,
                                   P_5_fp_rootl2_pos_mean_diff_norm_max),
                             min=c(P_5_fp_pos_mean_diff_norm_min,
                                   P_5_fp_pos_mean_diff_norm_min,
                                   P_5_fp_rootl2_pos_mean_diff_norm_min),
                             label=c(rep("2-probit",length(frac)),rep("uniform",length(frac)),rep("rootl2",length(frac))),
                             frac=c(rep(frac,3)))
P_5_meandiff
ggplot(P_5_meandiff, aes(frac, mean, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior mean (P=5)", colour="method", fill="method") +
  theme_bw()

P_5_fp_uniform_pos_cov_diff_norm_median <- c()
P_5_fp_uniform_pos_cov_diff_norm_max <- c()
P_5_fp_uniform_pos_cov_diff_norm_min <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("P_5_FP_uniform",frac[i],"chain_",j,sep = ""))$chain)-
                cov(P_5$beta.mean.chain),type = "2")
    
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  P_5_fp_uniform_pos_cov_diff_norm_median <- c(P_5_fp_uniform_pos_cov_diff_norm_median,median_result)
  P_5_fp_uniform_pos_cov_diff_norm_max <- c(P_5_fp_uniform_pos_cov_diff_norm_max,maxresult)
  P_5_fp_uniform_pos_cov_diff_norm_min <- c(P_5_fp_uniform_pos_cov_diff_norm_min,minresult)
  
}

P_5_fp_pos_cov_diff_norm_median <- c()
P_5_fp_pos_cov_diff_norm_max <- c()
P_5_fp_pos_cov_diff_norm_min <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("P_5_FP",frac[i],"chain_",j,sep = ""))$chain)-
                cov(P_5$beta.mean.chain),type = "2")
    
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  P_5_fp_pos_cov_diff_norm_median <- c(P_5_fp_pos_cov_diff_norm_median,median_result)
  P_5_fp_pos_cov_diff_norm_max <- c(P_5_fp_pos_cov_diff_norm_max,maxresult)
  P_5_fp_pos_cov_diff_norm_min <- c(P_5_fp_pos_cov_diff_norm_min,minresult)
  
}

P_5_fp_rootl2_cov_diff_norm_median <- c()
P_5_fp_rootl2_cov_diff_norm_max <- c()
P_5_fp_rootl2_cov_diff_norm_min <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("P_5_FP_rootl2",frac[i],"chain_",j,sep = ""))$chain)-
                cov(P_5$beta.mean.chain),type = "2")
    
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  P_5_fp_rootl2_cov_diff_norm_median <- c(P_5_fp_rootl2_cov_diff_norm_median,median_result)
  P_5_fp_rootl2_cov_diff_norm_max <- c(P_5_fp_rootl2_cov_diff_norm_max,maxresult)
  P_5_fp_rootl2_cov_diff_norm_min <- c(P_5_fp_rootl2_cov_diff_norm_min,minresult)
  
}

P_5_covdiff <- data.frame(
  median=c(P_5_fp_pos_cov_diff_norm_median,
           P_5_fp_uniform_pos_cov_diff_norm_median,
           P_5_fp_rootl2_cov_diff_norm_median),
  max=c(P_5_fp_pos_cov_diff_norm_max,
        P_5_fp_uniform_pos_cov_diff_norm_max,
        P_5_fp_rootl2_cov_diff_norm_max),
  min=c(P_5_fp_pos_cov_diff_norm_min,
        P_5_fp_uniform_pos_cov_diff_norm_min,
        P_5_fp_rootl2_cov_diff_norm_min),
  label=c(rep("2-probit",length(frac)),rep("uniform",length(frac)),rep("rootl2",length(frac))),
  frac=c(rep(frac,3)))
ggplot(P_5_covdiff, aes(frac, median, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior covariance (P=5)", colour="method", fill="method") +
  theme_bw()


######################P=8#####################

for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Compute_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=8)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y_p8[sample_index],initial_theda=MLE_model$coefficients,
                             true_theta=true_theta,Lp=8)
    assign(paste("P_8_FP",frac[i],"chain_",j,sep = ""), Model)
  }
}

for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Uniform_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=8)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y_p8[sample_index],initial_theda=MLE_model$coefficients,
                             true_theta=true_theta,Lp=8)
    assign(paste("P_8_FP_uniform",frac[i],"chain_",j,sep = ""), Model)
  }
}

for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Lp_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=8)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y_p8[sample_index],initial_theda=MLE_model$coefficients,
                             true_theta=true_theta,Lp=8)
    assign(paste("P_8_FP_rootl2",frac[i],"chain_",j,sep = ""), Model)
  }
}

P_8_fp_pos_mean_diff_norm_median <- c()
P_8_fp_pos_mean_diff_norm_max <- c()
P_8_fp_pos_mean_diff_norm_min <- c()
P_8_fp_pos_mean_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("P_8_FP",frac[i],"chain_",j,sep = ""))$post.mean-P_8$Beta.mean,type = "2")
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  P_8_fp_pos_mean_diff_norm_median <- c(P_8_fp_pos_mean_diff_norm_median,median_result)
  P_8_fp_pos_mean_diff_norm_max <- c(P_8_fp_pos_mean_diff_norm_max,maxresult)
  P_8_fp_pos_mean_diff_norm_min <- c(P_8_fp_pos_mean_diff_norm_min,minresult)
  P_8_fp_pos_mean_diff_norm_average <- c(P_8_fp_pos_mean_diff_norm_average,meanresult)
  
  #assign(paste("P_8_fp_pos_mean_diff_norm_",frac[i],"median",sep = ""), result)
  #assign(paste("P_8_fp_pos_mean_diff_norm_",frac[i],"max",sep = ""), maxresult)
  #assign(paste("P_8_fp_pos_mean_diff_norm_",frac[i],"min",sep = ""), minresult)
  
}

P_8_fp_uniform_pos_mean_diff_norm_median <- c()
P_8_fp_uniform_pos_mean_diff_norm_max <- c()
P_8_fp_uniform_pos_mean_diff_norm_min <- c()
P_8_fp_uniform_pos_mean_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("P_8_FP_uniform",frac[i],"chain_",j,sep = ""))$post.mean-P_8$Beta.mean,type = "2")
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  P_8_fp_uniform_pos_mean_diff_norm_median <- c(P_8_fp_uniform_pos_mean_diff_norm_median,median_result)
  P_8_fp_uniform_pos_mean_diff_norm_max <- c(P_8_fp_uniform_pos_mean_diff_norm_max,maxresult)
  P_8_fp_uniform_pos_mean_diff_norm_min <- c(P_8_fp_uniform_pos_mean_diff_norm_min,minresult)
  P_8_fp_uniform_pos_mean_diff_norm_average <- c(P_8_fp_uniform_pos_mean_diff_norm_average,meanresult)
  
  #assign(paste("P_8_fp_uniform_pos_mean_diff_norm_",frac[i],"median",sep = ""), result)
  #assign(paste("P_8_fp_uniform_pos_mean_diff_norm_",frac[i],"max",sep = ""), maxresult)
  #assign(paste("P_8_fp_uniform_pos_mean_diff_norm_",frac[i],"min",sep = ""), minresult)
  
}

P_8_fp_rootl2_pos_mean_diff_norm_median <- c()
P_8_fp_rootl2_pos_mean_diff_norm_max <- c()
P_8_fp_rootl2_pos_mean_diff_norm_min <- c()
P_8_fp_rootl2_pos_mean_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("P_8_FP_rootl2",frac[i],"chain_",j,sep = ""))$post.mean-P_8$Beta.mean,type = "2")
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  P_8_fp_rootl2_pos_mean_diff_norm_median <- c(P_8_fp_rootl2_pos_mean_diff_norm_median,median_result)
  P_8_fp_rootl2_pos_mean_diff_norm_max <- c(P_8_fp_rootl2_pos_mean_diff_norm_max,maxresult)
  P_8_fp_rootl2_pos_mean_diff_norm_min <- c(P_8_fp_rootl2_pos_mean_diff_norm_min,minresult)
  P_8_fp_rootl2_pos_mean_diff_norm_average <- c(P_8_fp_rootl2_pos_mean_diff_norm_average,meanresult)
  
  #assign(paste("P_8_fp_uniform_pos_mean_diff_norm_",frac[i],"median",sep = ""), result)
  #assign(paste("P_8_fp_uniform_pos_mean_diff_norm_",frac[i],"max",sep = ""), maxresult)
  #assign(paste("P_8_fp_uniform_pos_mean_diff_norm_",frac[i],"min",sep = ""), minresult)
  
}


P_8_meandiff <- data.frame(mean=c(P_8_fp_pos_mean_diff_norm_average,
                                    P_8_fp_uniform_pos_mean_diff_norm_average,
                                    P_8_fp_rootl2_pos_mean_diff_norm_average),
                             median=c(P_8_fp_pos_mean_diff_norm_median,
                                      P_8_fp_uniform_pos_mean_diff_norm_median,
                                      P_8_fp_rootl2_pos_mean_diff_norm_median),
                             max=c(P_8_fp_pos_mean_diff_norm_max,
                                   P_8_fp_uniform_pos_mean_diff_norm_max,
                                   P_8_fp_rootl2_pos_mean_diff_norm_max),
                             min=c(P_8_fp_pos_mean_diff_norm_min,
                                   P_8_fp_pos_mean_diff_norm_min,
                                   P_8_fp_rootl2_pos_mean_diff_norm_min),
                             label=c(rep("2-probit",length(frac)),rep("uniform",length(frac)),rep("rootl2",length(frac))),
                             frac=c(rep(frac,3)))
P_8_meandiff
ggplot(P_8_meandiff, aes(frac, mean, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior mean (P=8)", colour="method", fill="method") +
  theme_bw()

P_8_fp_uniform_pos_cov_diff_norm_median <- c()
P_8_fp_uniform_pos_cov_diff_norm_max <- c()
P_8_fp_uniform_pos_cov_diff_norm_min <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("P_8_FP_uniform",frac[i],"chain_",j,sep = ""))$chain)-
                cov(P_8$beta.mean.chain),type = "2")
    
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  P_8_fp_uniform_pos_cov_diff_norm_median <- c(P_8_fp_uniform_pos_cov_diff_norm_median,median_result)
  P_8_fp_uniform_pos_cov_diff_norm_max <- c(P_8_fp_uniform_pos_cov_diff_norm_max,maxresult)
  P_8_fp_uniform_pos_cov_diff_norm_min <- c(P_8_fp_uniform_pos_cov_diff_norm_min,minresult)
  
}

P_8_fp_pos_cov_diff_norm_median <- c()
P_8_fp_pos_cov_diff_norm_max <- c()
P_8_fp_pos_cov_diff_norm_min <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("P_8_FP",frac[i],"chain_",j,sep = ""))$chain)-
                cov(P_8$beta.mean.chain),type = "2")
    
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  P_8_fp_pos_cov_diff_norm_median <- c(P_8_fp_pos_cov_diff_norm_median,median_result)
  P_8_fp_pos_cov_diff_norm_max <- c(P_8_fp_pos_cov_diff_norm_max,maxresult)
  P_8_fp_pos_cov_diff_norm_min <- c(P_8_fp_pos_cov_diff_norm_min,minresult)
  
}

P_8_fp_rootl2_cov_diff_norm_median <- c()
P_8_fp_rootl2_cov_diff_norm_max <- c()
P_8_fp_rootl2_cov_diff_norm_min <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("P_8_FP_rootl2",frac[i],"chain_",j,sep = ""))$chain)-
                cov(P_8$beta.mean.chain),type = "2")
    
    medianchain<- rbind(medianchain,b^2)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  P_8_fp_rootl2_cov_diff_norm_median <- c(P_8_fp_rootl2_cov_diff_norm_median,median_result)
  P_8_fp_rootl2_cov_diff_norm_max <- c(P_8_fp_rootl2_cov_diff_norm_max,maxresult)
  P_8_fp_rootl2_cov_diff_norm_min <- c(P_8_fp_rootl2_cov_diff_norm_min,minresult)
  
}

P_8_covdiff <- data.frame(
  median=c(P_8_fp_pos_cov_diff_norm_median,
           P_8_fp_uniform_pos_cov_diff_norm_median,
           P_8_fp_rootl2_cov_diff_norm_median),
  max=c(P_8_fp_pos_cov_diff_norm_max,
        P_8_fp_uniform_pos_cov_diff_norm_max,
        P_8_fp_rootl2_cov_diff_norm_max),
  min=c(P_8_fp_pos_cov_diff_norm_min,
        P_8_fp_uniform_pos_cov_diff_norm_min,
        P_8_fp_rootl2_cov_diff_norm_min),
  label=c(rep("2-probit",length(frac)),rep("uniform",length(frac)),rep("rootl2",length(frac))),
  frac=c(rep(frac,3)))
ggplot(P_8_covdiff, aes(frac, median, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior covariance (P=8)", colour="method", fill="method") +
  theme_bw()


