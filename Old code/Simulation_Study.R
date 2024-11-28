set.seed(2)
rm(list = ls())
options(scipen = 8)

library(abind)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("Functions.R")

dev.off()
N <- 2000
#D should always be odd number since we use (D-1)/2 to generate the variables
D <- 11
N_sim <- 2000
burn_in<- 1500
#gengerate data and true parameters
X <- generate.data(N,D)/10
# True values of regression coeffiecients theta
true_theta <- runif(D,-3,3)
Xb <- X%*%true_theta
#####This step should be repeated several times to avoid the apperance of the
#perfectly separation.
hist(Xb)
#Obtain the true probability of pi and the real bernouli distributed Y
pi.logit <- exp(Xb)/(1+exp(Xb))
pi.probit <- pnorm(Xb)
pi.cloglog <- 1-exp(-exp(Xb))
#Here we check if the data are perfectly linear separation. If it's then the
#data need to be regenerated
hist(pi.logit)
hist(pi.probit)
hist(pi.cloglog)
#Assign ys
y.logit <- rbinom(n=N, size=1, prob=pi.logit)
y.probit <- rbinom(N, 1, pi.probit)
y.cloglog <- rbinom(N, 1, pi.cloglog)

#sample_index<- Compute_coreset(sketch_size = 100,coreset_size = 100,X=X,Lp=2)


#Using MLE estimator to otain the maximum likelihood result
#Assign a as a symbol states NxxKDxx,which represents the sample size and number of variables.

MLE_probit <- glm(y.probit~X[,-1],family = binomial(link = probit))
MLE_logit <- glm(y.logit~X[,-1],family = binomial(link = logit))

MLE_probit_coefficients <- MLE_probit$coefficients

summary(MLE_logit)
summary(MLE_probit)
#Setting LP
Lp <- c(0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8)

for (i in Lp) {
  pi_p <- pgnorm(Xb,alpha = p_scale(i),beta = i)
  assign(paste("pi_",i,sep = ""),pi_p)

}
#Setting Yp


for (i in Lp) {
  y_pi <- generate_p_y(p=i,Xb=Xb,N=N)
  assign(paste("y_p",i,sep = ""),y_pi)

}

#frac <- c(seq(0.001,0.01,0.001))
#sample_index<- Compute_coreset(sketch_size = 100,coreset_size = 100,X=X,Lp=2)
#sample logit model
P_logit <- multi_chain(N_sim=200,burn_in=100,
                                    X=X, y=y.logit,initial_theda=MLE_logit$coefficients,
                                    true_theta=true_theta,
                       Lp=1,M_iter=100,
                       range=c(0.1,5),step=0.01,
                       times = 1)


#Convergence check
mcmc_logit <- lapply(P_logit$Lp_chain, mcmc)
gelman_res <- gelman.diag(mcmc_logit,autoburnin = F)
gelman.plot(mcmc_logit, autoburnin = F)
gelman_res

plot(P_logit$Lp.mean.chain[-1],type = "l",main = "Estimated MCMC chain of p with metropolis-hasting", xlab = "iteration",ylab = "p" )
Msubtitle <- "P=logit"
mtext(side = 3, line = 0.4, Msubtitle)

#sample probit model
P_probit <- multi_chain(N_sim=2000,burn_in=1000,
                                    X=X, y=y.probit,initial_theda=MLE_probit_coefficients,
                                    true_theta=true_theta,Lp=2,M_iter=100,range=c(0.1,5),step=0.01,
                        times = 5)
P_probit$Beta.mean
P_probit$Lp.mean
#Convergence check
mcmc_probit <- lapply(P_probit$Lp_chain, mcmc)
gelman_res <- gelman.diag(mcmc_probit,autoburnin = FALSE)
gelman.plot(mcmc_probit,autoburnin = FALSE)
gelman_res

plot(P_probit$Lp.mean.chain[-1],type = "l",main = "Estimated MCMC chain of p with metropolis-hasting",
     xlab = "iteration",ylab = "p" ,ylim=c(1.5,2.5))
Msubtitle <- "P=probit, initial p=5"
mtext(side = 3, line = 0.4, Msubtitle)
#Sample p0.5 model
# plot(1:10000,test2$Lp,type = "l")
P_0.5 <- multi_chain(times=5,N_sim=2000,burn_in=1000,
                                    X=X, y=y_p0.5,initial_theda=MLE_probit_coefficients,
                                    true_theta=true_theta,Lp=2,M_iter=50,range=c(0.1,5),step=0.01)
P_0.5$Lp.mean
plot(P_0.5$Lp.mean.chain[-1],type = "l",main = "Estimated MCMC chain of p with metropolis-hasting", xlab = "iteration",ylab = "p" )
Msubtitle <- "P=0.5, initial p=5"
mtext(side = 3, line = 0.4, Msubtitle)

mcmc_0.5 <- lapply(P_0.5$Lp_chain, mcmc)
gelman_res <- gelman.diag(mcmc_0.5, autoburnin = FALSE)
gelman.plot(mcmc_0.5,autoburnin = FALSE)
gelman_res

P_1 <- multi_chain(times=5,N_sim=2000,burn_in=1000,
                                    X=X,
                   y=y_p1,initial_theda=MLE_probit_coefficients,
                                    true_theta=true_theta,Lp=5,
                   M_iter=100,range=c(0.1,5),step=0.01)

plot(P_1$Lp.mean.chain[-1],type = "l",main = "Estimated MCMC chain of p with metropolis-hasting", xlab = "iteration",ylab = "p" )
Msubtitle <- "P=1, initial p=5"
mtext(side = 3, line = 0.4, Msubtitle)

mcmc_1 <- lapply(P_1$Lp_chain, mcmc)
gelman_res <- gelman.diag(mcmc_1, autoburnin = FALSE)
gelman.plot(mcmc_1, autoburnin = F)
gelman_res

P_1.5 <- multi_chain(times=5,N_sim=2000,burn_in=1000,
                                    X=X, y=y_p1.5,initial_theda=MLE_probit_coefficients,
                                    true_theta=true_theta,Lp=5,M_iter=200,range=c(0.1,5),step=0.05)
P_1.5
P_1.5$Lp.mean
plot(P_1.5$Lp.mean.chain[-1],type = "l",main = "Estimated MCMC chain of p with metropolis-hasting", xlab = "iteration",ylab = "p" )
Msubtitle <- "P=1.5, initial p=5"
mtext(side = 3, line = 0.4, Msubtitle)
mcmc_1.5 <- lapply(P_1.5$Lp_chain, mcmc)
gelman_res <- gelman.diag(mcmc_1.5, autoburnin = FALSE)
gelman.plot(mcmc_1, autoburnin = F)
gelman_res
#
# P_2 <- multi_chain(times=1,N_sim=N_sim,burn_in=burn_in,
#                                     X=X, y=y_p2,initial_theda=MLE_probit_coefficients,
#                                     true_theta=true_theta,Lp=5,M_iter=200,range=c(0.1,5),step=0.05)
# P_2
# P_2$Lp.mean
#
# plot(P_2$Lp.mean.chain[-1],type = "l",main = "Estimated MCMC chain of p with metropolis-hasting", xlab = "iteration",ylab = "p" )
# Msubtitle <- "P=2, initial p=5"
# mtext(side = 3, line = 0.4, Msubtitle)

#
# P_2.5 <- multi_chain(times=1,N_sim=N_sim,burn_in=burn_in,
#                                     X=X, y=y_p2.5,initial_theda=MLE_probit_coefficients,
#                                     true_theta=true_theta,Lp=5,M_iter=200,range=c(0.1,5),step=0.05)
# P_2.5
# P_2.5$Lp.mean
# plot(P_2.5$Lp.mean.chain[-1],type = "l",main = "Estimated MCMC chain of p with metropolis-hasting", xlab = "iteration",ylab = "p" )
# Msubtitle <- "P=2.5, initial p=5"
# mtext(side = 3, line = 0.4, Msubtitle)

P_3 <- multi_chain(times=5,N_sim=N_sim,burn_in=burn_in,
                                    X=X, y=y_p3,
                   initial_theda=MLE_probit_coefficients,
                   true_theta=true_theta,Lp=5,M_iter=200,
                   range=c(0.1,5),step=0.01)
P_3
P_3$Lp.mean
plot(P_3$Lp.mean.chain[-1],type = "l",main = "Estimated MCMC chain of p with metropolis-hasting", xlab = "iteration",ylab = "p" )
Msubtitle <- "P=3, initial p=5"
mtext(side = 3, line = 0.4, Msubtitle)

mcmc_3 <- lapply(P_3$Lp_chain, mcmc)
gelman_res <- gelman.diag(mcmc_3, autoburnin = FALSE)
gelman.plot(mcmc_3)
gelman_res

# P_3.5 <- multi_chain(times=1,N_sim=N_sim,burn_in=burn_in,
#                                     X=X, y=y_p3.5,initial_theda=MLE_probit_coefficients,
#                                     true_theta=true_theta,Lp=5,M_iter=200,range=c(0.1,5),step=0.05)
# P_3.5
# P_3.5$Lp.mean
# plot(P_3.5$Lp.mean.chain[-1],type = "l",main = "Estimated MCMC chain of p with metropolis-hasting", xlab = "iteration",ylab = "p" )
# Msubtitle <- "P=3.5, initial p=5"
# mtext(side = 3, line = 0.4, Msubtitle)



P_4 <- multi_chain(times=5,N_sim=N_sim,burn_in=burn_in,
                                     X=X, y=y_p4,initial_theda=MLE_probit_coefficients,
                                     true_theta=true_theta,Lp=5,M_iter=200,range=c(0.1,5),step=0.05)
P_4
P_4$Lp.mean
plot(P_4$Lp.mean.chain[-1],type = "l",main = "Estimated MCMC chain of p with metropolis-hasting", xlab = "iteration",ylab = "p" )
Msubtitle <- "P=4, initial p=5"
mtext(side = 3, line = 0.4, Msubtitle)
mcmc_4 <- lapply(P_4$Lp_chain, mcmc)
gelman_res <- gelman.diag(mcmc_4, autoburnin = FALSE)
gelman.plot(mcmc_4)
gelman_res

# P_4.5 <- multi_chain(times=1, N_sim=N_sim,burn_in=burn_in,
#                                      X=X, y=y_p4.5,initial_theda=MLE_probit_coefficients,
#                                      true_theta=true_theta,Lp=5,M_iter=200,range=c(0.1,10),step=0.01)
#
# P_4.5
# P_4.5$Lp.mean
# plot(P_4.5$Lp.mean.chain[-1],type = "l",main = "Estimated MCMC chain of p with metropolis-hasting", xlab = "iteration",ylab = "p" )
# Msubtitle <- "P=4.5, initial p=5"
# mtext(side = 3, line = 0.4, Msubtitle)


#5 m200
P_5 <- multi_chain(times=5, N_sim=N_sim,burn_in=burn_in,
                                     X=X, y=y_p5,initial_theda=MLE_probit_coefficients,
                                     true_theta=true_theta,Lp=5,M_iter=100,range=c(0.1,8),step=0.05)


P_5
P_5$Lp.mean
plot(P_5$Lp.mean.chain[-1],type = "l",main = "Estimated MCMC chain of p with metropolis-hasting", xlab = "iteration",ylab = "p" )
Msubtitle <- "P=5, initial p=5"
mtext(side = 3, line = 0.4, Msubtitle)
mcmc_5 <- lapply(P_5$Lp_chain, mcmc)
gelman_res <- gelman.diag(mcmc_5, autoburnin = FALSE)
gelman.plot(mcmc_5,autoburnin = FALSE)
gelman_res


P_8 <- multi_chain(times=5,N_sim=2000,burn_in=1000,
                                     X=X, y=y_p8,initial_theda=MLE_probit_coefficients,
                                     true_theta=true_theta,Lp=5,
                   M_iter=200,range=c(0.1,10),
                   step=0.05)
P_8
P_8$Lp.mean
plot(P_8$Lp.mean.chain[-1],type = "l",main = "Estimated MCMC chain of p with metropolis-hasting", xlab = "iteration",ylab = "p" )
Msubtitle <- "P=8, initial p=5"
mtext(side = 3, line = 0.4, Msubtitle)

mcmc_8 <- lapply(P_8$Lp_chain, mcmc)
gelman_res <- gelman.diag(mcmc_8, autoburnin = FALSE)
gelman.plot(mcmc_8, autoburnin = F)
gelman_res


#
# plot_a <- ggplot(data=data.frame(Lp.mean.chain=P_logit$Lp.mean.chain[1:N_sim],
#                                  x=1:length(P_logit$Lp.mean.chain)),
#                  aes(x=x, y=Lp.mean.chain)) +
#   geom_line()+
#   theme_bw(base_size = 12) +
#   labs(title="The estimated p for the logit data",
#        # subtitle="Estimated p=0.81",
#        x="MCMC iteration",
#        y="p")+ylim(c(0.5,2))
#
# plot_b <- ggplot(data=data.frame(Lp.mean.chain=P_probit$Lp.mean.chain[1:N_sim],
#                                  x=1:length(P_probit$Lp.mean.chain)),
#                  aes(x=x, y=Lp.mean.chain)) +
#   geom_line()+
#   theme_bw(base_size = 12) +
#   labs(title="The estimated p for the probit data",
#        # subtitle="Estimated p=1.94",
#        x="MCMC iteration",
#        y="p")+ylim(c(1,3))
#
# plot_c <- ggplot(data=data.frame(Lp.mean.chain=P_1$Lp.mean.chain[1:N_sim],
#                                  x=1:length(P_1$Lp.mean.chain)),
#                  aes(x=x, y=Lp.mean.chain)) +
#   geom_line()+
#   theme_bw(base_size = 12) +
#   labs(title="The estimated p for the p=1 data",
#        # subtitle="Estimated p=1.00",
#        x="MCMC iteration",
#        y="p")+ylim(c(0.5,4))
#
# plot_d <- ggplot(data=data.frame(Lp.mean.chain=P_2$Lp.mean.chain[1:N_sim],
#                                  x=1:length(P_2$Lp.mean.chain)),
#                  aes(x=x, y=Lp.mean.chain)) +
#   geom_line()+
#   theme_bw(base_size = 12) +
#   labs(title="The estimated p for the p=2 data",
#        # subtitle="Estimated p=2.06",
#        x="MCMC iteration",
#        y="p")+ylim(c(1.5,4))
#
# ggarrange(plot_a, plot_b, plot_c ,plot_d,
#           labels = c("A", "B", "C","D"),
#           ncol = 2, nrow = 2)
#
#
# jpeg("Simulation_p_logit_probit_p=2(version2).jpeg", units="in", width=8, height=5, res=300)
# ggarrange(plot_a, plot_b, plot_c ,plot_d,
#           labels = c("A", "B", "C","D"),
#           ncol = 2, nrow = 2)
# jpeg("p_chain_probit_N50K.jpeg", units="in", width=8, height=5, res=300)
#Plot the chain
start <- 300
chain_list <- P_probit$Lp_chain
chain_df <- data.frame(do.call(cbind, chain_list))
chain_df <- chain_df[1:2000,]
plotdata <- data.frame(mean=apply(chain_df[-(1:start),], MARGIN = 1, mean),
                       max=apply(chain_df[-(1:start),], MARGIN = 1, max),
                       min=apply(chain_df[-(1:start),], MARGIN = 1, min),
                       iteration = seq(start+1,nrow(chain_df)))
ggplot(plotdata, aes(iteration, mean)) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.4, colour=NA) +
  geom_line() +
  labs(x="MCMC Iteration",
       y="p") +
  annotate("text", x=700, y=2.7,size=6 ,
           label= TeX("$\\hat{p} = 2.033",
                      output='character'),parse=TRUE) +
  annotate("text", x=700, y=2.9,size=6 ,
           label= "Estimation of p for probit data") +
  theme_bw(base_size = 12) +
  ylim(c(1.5,3))
# dev.off()

#jpeg("PSRF_logit_N50K.jpeg", units="in", width=8, height=5, res=300)
gelman.plot(mcmc_logit, autoburnin = F,main="PSRF of MCMC chain for logit data \n PSRF=1.06",
            xlab = "Iteration",
            ylab = "PSRF")
dev.off()

# jpeg("p_chain_logit_n50k.jpeg", units="in", width=8, height=5, res=300)
start <- 300
chain_list <- P_logit$Lp_chain
chain_df <- data.frame(do.call(cbind, chain_list))
plotdata <- data.frame(mean=apply(chain_df[-(1:start),], MARGIN = 1, mean),
                       max=apply(chain_df[-(1:start),], MARGIN = 1, max),
                       min=apply(chain_df[-(1:start),], MARGIN = 1, min),
                       iteration = seq(start+1,nrow(chain_df)))
ggplot(plotdata, aes(iteration, mean)) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.4, colour=NA) +
  geom_line() +
  labs(x="MCMC Iteration",
       y="p") +
  annotate("text", x=7500, y=1.0,size=6 ,
           label= TeX("$\\hat{p} = 1.472",
                      output='character'),parse=TRUE) +
  annotate("text", x=7500, y=1.2,size=6 ,
           label= "Estimation of p for logit data") +
  theme_bw(base_size = 12) +
  ylim(c(0.5,2))
dev.off()


# jpeg("PSRF_probit_N50K.jpeg", units="in", width=8, height=5, res=300)
mcmc_probit <- lapply(P_probit$Lp_chain, mcmc, end=2000)
gelman.diag(mcmc_probit,autoburnin = F)
gelman.plot(mcmc_probit, autoburnin = F,main="PSRF of MCMC chain for probit data \n PSRF=1.11",
            xlab = "Iteration",
            ylab = "PSRF")
dev.off()
# ggplot(data=data.frame(Lp.mean.chain=P_probit$Lp.mean.chain[seq(50,1000,5)],
#                        x=seq(50,1000,5)),
#        aes(x=x, y=Lp.mean.chain)) +
#   geom_line()+
#   annotate("text", x=700, y=2.7,size=6 ,label= TeX("$\\hat{p} = 2.033", output='character'),parse=TRUE) +
#   annotate("text", x=700, y=2.9,,size=6 ,label= "Estimation of p for probit data") +
#   theme_bw(base_size = 12) +
#   labs(
#     # title="Estimation of p for p=3 data",
#     # subtitle="Estimated p=3.03",
#     x="MCMC iteration",
#     y="p")+ylim(c(1.5,3))


# ggplot(data=data.frame(Lp.mean.chain=P_logit$Lp.mean.chain[seq(50,5000,20)],
#                        x=seq(50,5000,20)),
#        aes(x=x, y=Lp.mean.chain)) +
#   geom_line()+
#   annotate("text", x=1300, y=2.89,,size=6, label= TeX("$\\hat{p} = 1.647", output='character'),parse=TRUE) +
#   annotate("text", x=1300, y=2.7,,size=6 ,label= "Estimation of p for logit data") +
#
#   theme_bw(base_size = 12) +
#   labs(
#     # title="Estimation of p for p=3 data",
#     # subtitle="Estimated p=3.03",
#     x="MCMC iteration",
#     y="p")+ylim(c(0.5,3))
# dev.off()
#checking if the posterior of p can be integral to 1.
calculate_and_plot_posterior(P_logit$Lp.mean.chain[500:2000])
calculate_and_plot_posterior(P_probit$Lp.mean.chain[500:2000])
calculate_and_plot_posterior(P_0.5$Lp.mean.chain[500:2000])
calculate_and_plot_posterior(P_1$Lp.mean.chain[500:2000])
calculate_and_plot_posterior(P_1.5$Lp.mean.chain[500:2000])
calculate_and_plot_posterior(P_3$Lp.mean.chain[500:2000])
calculate_and_plot_posterior(P_4$Lp.mean.chain[500:2000])
calculate_and_plot_posterior(P_5$Lp.mean.chain[500:2000])
calculate_and_plot_posterior(P_8$Lp.mean.chain[500:2000])
dflogit <- calculate_and_plot_posterior(P_logit$Lp.mean.chain[500:10000])$posterior_data
dfprobit <- calculate_and_plot_posterior(P_probit$Lp.mean.chain[500:2000])$posterior_data
df0.5 <- calculate_and_plot_posterior(P_0.5$Lp.mean.chain[500:2000])$posterior_data
df1 <- calculate_and_plot_posterior(P_1$Lp.mean.chain[500:2000])$posterior_data
df1.5 <- calculate_and_plot_posterior(P_1.5$Lp.mean.chain[500:2000])$posterior_data
df3 <- calculate_and_plot_posterior(P_3$Lp.mean.chain[500:2000])$posterior_data
df4 <- calculate_and_plot_posterior(P_4$Lp.mean.chain[500:2000])$posterior_data
df5 <- calculate_and_plot_posterior(P_5$Lp.mean.chain[500:2000])$posterior_data
df8 <- calculate_and_plot_posterior(P_8$Lp.mean.chain[500:10000])$posterior_data
ggplot(df0.5, aes(x = Values, y = Density)) +
  geom_line() +xlim(0.5,0.55)+
  labs(title = "Posterior distribution of p=0.5", x = "Values", y = "Density")
# jpeg("posterior_logit.jpeg", units="in", width=8, height=5, res=300)
ggplot(dflogit, aes(x = Values, y = Density)) +
  geom_line()+
  labs(title = "Posterior distribution of logit data", x = "Values", y = "Density")+
  xlim(1.3,1.55)
dev.off()
# jpeg("posterior_probit.jpeg", units="in", width=8, height=5, res=300)
ggplot(dfprobit, aes(x = Values, y = Density)) +
  geom_line()+
  labs(title = "Posterior distribution of probit data", x = "Values", y = "Density")+
  xlim(1.9,2.1)
dev.off()
# jpeg("posterior_p1.jpeg", units="in", width=8, height=5, res=300)
ggplot(df1, aes(x = Values, y = Density)) +
  geom_line()+
  labs(title = "Posterior distribution of p=1 data", x = "Values", y = "Density")+
  xlim(1,1.1)
dev.off()
ggplot(df1.5, aes(x = Values, y = Density)) +
  geom_line()+
  labs(title = "Posterior distribution of p=1.5 data", x = "Values", y = "Density")+
  xlim(1.2,1.6)
# jpeg("posterior_p3.jpeg", units="in", width=8, height=5, res=300)
ggplot(df3, aes(x = Values, y = Density)) +
  geom_line()+
  labs(title = "Posterior distribution of p=3 data", x = "Values", y = "Density")+
  xlim(2.7,3.2)
dev.off()
# jpeg("posterior_p5.jpeg", units="in", width=8, height=5, res=300)
ggplot(df5, aes(x = Values, y = Density)) +
  geom_line()+
  labs(title = "Posterior distribution of p=5 data", x = "Values", y = "Density")+
  xlim(4.8,6.2)
dev.off()
# jpeg("posterior_p8.jpeg", units="in", width=8, height=5, res=300)
ggplot(df8, aes(x = Values, y = Density)) +
  geom_line()+
  labs(title = "Posterior distribution of p=8 data", x = "Values", y = "Density")+
  xlim(7,10)
dev.off()

jpeg("PSRF_logit_N50K.jpeg", units="in", width=8, height=5, res=300)
mcmc_logit <- lapply(P_logit$Lp_chain, mcmc)
gelman.diag(mcmc_logit,autoburnin = F)
gelman.plot(mcmc_logit, autoburnin =F,main="PSRF of MCMC chain for logit data \n PSRF=1.00",
            xlab = "Iteration",
            ylab = "PSRF")
dev.off()

jpeg("PSRF_probit_N50K.jpeg", units="in", width=8, height=5, res=300)
mcmc_probit <- lapply(P_probit$Lp_chain, mcmc, end=2000)
gelman.diag(mcmc_probit,autoburnin = F)
gelman.plot(mcmc_probit, autoburnin =F,main="PSRF of MCMC chain for probit data \n PSRF=1.12",
            xlab = "Iteration",
            ylab = "PSRF", xlim =c(0,2000))
dev.off()
jpeg("PSRF_p1_N50K.jpeg", units="in", width=8, height=5, res=300)
gelman.diag(mcmc_1,autoburnin = F)
gelman.plot(mcmc_1, autoburnin = F,main="PSRF of MCMC chain for p=1 data \n PSRF=1.00",
            xlab = "Iteration",
            ylab = "PSRF")
dev.off()

jpeg("PSRF_p3_N50K.jpeg", units="in", width=8, height=5, res=300)
gelman.diag(mcmc_3,autoburnin = F)
gelman.plot(mcmc_3, autoburnin = F,main="PSRF of MCMC chain for p=3 data \n PSRF=1.16",
            xlab = "Iteration",
            ylab = "PSRF")
dev.off()

jpeg("PSRF_p5_N50K.jpeg", units="in", width=8, height=5, res=300)
gelman.diag(mcmc_5,autoburnin = F)
gelman.plot(mcmc_5, autoburnin = F,main="PSRF of MCMC chain for p=5 data \n PSRF=1.07",
            xlab = "Iteration",
            ylab = "PSRF")
dev.off()

jpeg("PSRF_p8_N50K.jpeg", units="in", width=8, height=5, res=300)
mcmc_8 <- lapply(P_8$Lp_chain, mcmc)
gelman.diag(mcmc_8,autoburnin = F)
gelman.plot(mcmc_8, autoburnin = F,main="PSRF of MCMC chain for p=8 data \n PSRF=1.27",
            xlab = "Iteration",
            ylab = "PSRF")
dev.off()


# jpeg("p_chain_1_n50k.jpeg", units="in", width=8, height=5, res=300)
start <- 300
chain_list <- P_1$Lp_chain
chain_df <- data.frame(do.call(cbind, chain_list))
plotdata <- data.frame(mean=apply(chain_df[-(1:start),], MARGIN = 1, mean),
                       max=apply(chain_df[-(1:start),], MARGIN = 1, max),
                       min=apply(chain_df[-(1:start),], MARGIN = 1, min),
                       iteration = seq(start+1,nrow(chain_df)))
ggplot(plotdata, aes(iteration, mean)) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.4, colour=NA) +
  geom_line() +
  labs(x="MCMC Iteration",
       y="p") +
  annotate("text", x=1500, y=1.3,size=6 ,
           label= TeX("$\\hat{p} = 1.069",
                      output='character'),parse=TRUE) +
  annotate("text", x=1500, y=1.4,size=6 ,
           label= "Estimation of p for p=1 data") +
  theme_bw(base_size = 12) +
  ylim(c(0.75,1.5))
dev.off()


# jpeg("p_chain_1_n10k.jpeg", units="in", width=8, height=5, res=300)
# ggplot(data=data.frame(Lp.mean.chain=P_1$Lp.mean.chain[seq(50,1000,5)],
#                        x=seq(50,1000,5)),
#        aes(x=x, y=Lp.mean.chain)) +
#   geom_line()+
#   annotate("text", x=700, y=1.7,size=6, label= TeX("$\\hat{p} = 0.930", output='character'),parse=TRUE) +
#   annotate("text", x=700, y=1.6, ,size=6,label= "Estimation of p for p=1 data") +
#
#   theme_bw(base_size = 12) +
#   labs(
#     # title="Estimation of p for p=3 data",
#     # subtitle="Estimated p=3.03",
#     x="MCMC iteration",
#     y="p")+ylim(c(0.5,2))
# dev.off()

# jpeg("p=2.jpeg", units="in", width=8, height=5, res=300)
#
# ggplot(data=data.frame(Lp.mean.chain=P_2$Lp.mean.chain[seq(50,1000,5)],
#                        x=seq(50,1000,5)),
#        aes(x=x, y=Lp.mean.chain)) +
#   geom_line()+
#   annotate("text", x=800, y=2.7,size=6, label= TeX("$\\hat{p} = 1.984", output='character'),parse=TRUE) +
#   annotate("text", x=800, y=2.9, size=6, label= "Estimation of p for p=2 data") +
#
#   theme_bw(base_size = 12) +
#   labs(
#     # title="Estimation of p for p=3 data",
#     # subtitle="Estimated p=3.03",
#     x="MCMC iteration",
#     y="p")+ylim(c(1.5,3))
# dev.off()


# jpeg("p_chain_3_n50k.jpeg", units="in", width=8, height=5, res=300)
start <- 300
chain_list <- P_3$Lp_chain
chain_df <- data.frame(do.call(cbind, chain_list))
plotdata <- data.frame(mean=apply(chain_df[-(1:start),], MARGIN = 1, mean),
                       max=apply(chain_df[-(1:start),], MARGIN = 1, max),
                       min=apply(chain_df[-(1:start),], MARGIN = 1, min),
                       iteration = seq(start+1,nrow(chain_df)))
ggplot(plotdata, aes(iteration, mean)) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.4, colour=NA) +
  geom_line() +
  labs(x="MCMC Iteration",
       y="p") +
  annotate("text", x=1500, y=2.1,size=6 ,
           label= TeX("$\\hat{p} = 3.021",
                      output='character'),parse=TRUE) +
  annotate("text", x=1500, y=2.3,size=6 ,
           label= "Estimation of p for p=3 data") +
  theme_bw(base_size = 12) +
  ylim(c(2,4))
dev.off()

#
# jpeg("p_chain_3_n10k.jpeg", units="in", width=8, height=5, res=300)
#
#   ggplot(data=data.frame(Lp.mean.chain=P_3$Lp.mean.chain[seq(50,1000,5)],
#                                  x=seq(50,1000,5)),
#                  aes(x=x, y=Lp.mean.chain)) +
#   geom_line()+
#   annotate("text", x=800, y=2.9,size=6, label= TeX("$\\hat{p} = 4.039", output='character'),parse=TRUE) +
#   annotate("text", x=800, y=3.15,size=6, label= "Estimation of p for p=3 data") +
#
#   theme_bw(base_size = 12) +
#   labs(
#     # title="Estimation of p for p=3 data",
#     # subtitle="Estimated p=3.03",
#        x="MCMC iteration",
#        y="p")+ylim(c(2.5,5.5))
#   dev.off()

  # jpeg("p_chain_4_n10k.jpeg", units="in", width=8, height=5, res=300)
  #
  #
  # ggplot(data=data.frame(Lp.mean.chain=P_4$Lp.mean.chain[seq(50,1000,5)],
  #                        x=seq(50,1000,5)),
  #        aes(x=x, y=Lp.mean.chain)) +
  #   geom_line()+
  #   annotate("text", x=700, y=4.9,size=6, label= TeX("$\\hat{p} = 4.070", output='character'),parse=TRUE) +
  #   annotate("text", x=700, y=4.7,size=6, label= "Estimation of p for p=4 data") +
  #
  #   theme_bw(base_size = 12) +
  #   labs(
  #     # title="Estimation of p for p=3 data",
  #     # subtitle="Estimated p=3.03",
  #     x="MCMC iteration",
  #     y="p")+ylim(c(3.5,5))
  #
  # dev.off()

# jpeg("p_chain_5_n50k.jpeg", units="in", width=8, height=5, res=300)
start <- 300
chain_list <- P_5$Lp_chain
chain_df <- data.frame(do.call(cbind, chain_list))
plotdata <- data.frame(mean=apply(chain_df[-(1:start),], MARGIN = 1, mean),
                       max=apply(chain_df[-(1:start),], MARGIN = 1, max),
                       min=apply(chain_df[-(1:start),], MARGIN = 1, min),
                       iteration = seq(start+1,nrow(chain_df)))
ggplot(plotdata, aes(iteration, mean)) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.4, colour=NA) +
  geom_line() +
  labs(x="MCMC Iteration",
       y="p") +
  annotate("text", x=1500, y=4.2,size=6 ,
           label= TeX("$\\hat{p} = 5.507",
                      output='character'),parse=TRUE) +
  annotate("text", x=1500, y=4.5,size=6 ,
           label= "Estimation of p for p=5 data") +
  theme_bw(base_size = 12) +
  ylim(c(4,7))
dev.off()

# jpeg("p_chain_8_n50k.jpeg", units="in", width=8, height=5, res=300)
start <- 50
chain_list <- P_8$Lp_chain
chain_df <- data.frame(do.call(cbind, chain_list))
plotdata <- data.frame(mean=apply(chain_df[-(1:start),], MARGIN = 1, mean),
                       max=apply(chain_df[-(1:start),], MARGIN = 1, max),
                       min=apply(chain_df[-(1:start),], MARGIN = 1, min),
                       iteration = seq(start+1,nrow(chain_df)))
ggplot(plotdata, aes(iteration, mean)) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.4, colour=NA) +
  geom_line() +
  labs(x="MCMC Iteration",
       y="p") +
  annotate("text", x=3500, y=10,size=6 ,
           label= TeX("$\\hat{p} = 8.525",
                      output='character'),parse=TRUE) +
  annotate("text", x=3500, y=10.4,size=6 ,
           label= "Estimation of p for p=8 data") +
  theme_bw(base_size = 12) +
  ylim(c(7,11))
dev.off()


  # jpeg("p_chain_5_n10k.jpeg", units="in", width=8, height=5, res=300)
  #
  # ggplot(data=data.frame(Lp.mean.chain=P_5$Lp.mean.chain[seq(50,1000,5)],
  #                        x=seq(50,1000,5)),
  #        aes(x=x, y=Lp.mean.chain)) +
  #   geom_line()+
  #   annotate("text", x=700, y=5.4,size=6, label= TeX("$\\hat{p} = 4.234 ", output='character'),parse=TRUE) +
  #   annotate("text", x=700, y=5.6,size=6, label= "Estimation of p for p=5 data") +
  #
  #   theme_bw(base_size = 12) +
  #   labs(
  #     # title="Estimation of p for p=3 data",
  #     # subtitle="Estimated p=3.03",
  #     x="MCMC iteration",
  #     y="p")+ylim(c(4,6.5))
  # dev.off()

  # jpeg("p_chain_8_n10k.jpeg", units="in", width=8, height=5, res=300)

  ggplot(data=data.frame(Lp.mean.chain=P_8$Lp.mean.chain[seq(50,10000,20)],
                                                    x=seq(50,10000,20)),
                                    aes(x=x, y=Lp.mean.chain)) +
    geom_line()+
    annotate("text", x=7500, y=6.7,size=6, label= TeX("$\\hat{p} = 9.351", output='character'),parse=TRUE) +
    annotate("text", x=7500, y=7.2,size=6, label= "Estimation of p for p=8 data") +

    theme_bw(base_size = 12) +
    labs(
      # title="Estimation of p for p=3 data",
      # subtitle="Estimated p=3.03",
      x="MCMC iteration",
      y="p")+ylim(c(6.5,10.5))
  dev.off()

#
# plot_b <- ggplot(data=data.frame(Lp.mean.chain=P_4$Lp.mean.chain[1:N_sim],
#                                  x=1:length(P_4$Lp.mean.chain)),
#                  aes(x=x, y=Lp.mean.chain)) +
#   geom_line()+
#   theme_bw(base_size = 12) +
#   labs(title="Estimation of p for p=4 data", subtitle="Estimated p=3.68", x="MCMC iteration",
#        y="p")+ylim(c(3,5))

# plot_c <- ggplot(data=data.frame(Lp.mean.chain=P_5$Lp.mean.chain[1:N_sim],
#                                  x=1:length(P_5$Lp.mean.chain)),
#                  aes(x=x, y=Lp.mean.chain)) +
#   geom_line()+
#   theme_bw(base_size = 12) +
#   labs(title="Estimation p for p=5 data", subtitle="Estimated p=4.13", x="MCMC iteration",
#        y="p")+ylim(c(4,6))
#
# plot_d <- ggplot(data=data.frame(Lp.mean.chain=P_8$Lp.mean.chain[1:N_sim],
#                                  x=1:length(P_8$Lp.mean.chain)),
#                  aes(x=x, y=Lp.mean.chain)) +
#   geom_line()+
#   theme_bw(base_size = 12) +
#   labs(title="The estimated p for the p=8 data", subtitle="Estimated p=7.25", x="MCMC iteration",
#        y="p")+ylim(c(6,10))

# ggarrange(plot_a, plot_b, plot_c ,plot_d,
#           labels = c("A", "B", "C","D"),
#           ncol = 2, nrow = 2)
#
# jpeg("Simulation_p_5_8.jpeg(version2)", units="in", width=8, height=5, res=300)
# ggarrange(plot_a, plot_b, plot_c ,plot_d,
#           labels = c("A", "B", "C","D"),
#           ncol = 2, nrow = 2)
# dev.off()
#
# probit_llk_ratio<- c()
# for (i in 1:10) {
#
#   b<- llk(X, credit_p2$post.mean,2,)-
#                    llk(X,get(paste("credit_p_2_2probit_frac",frac[i],"chain_",j,sep = ""))$post.mean,2)
#
# }






# plot_a <- ggplot() +
#   geom_line(data=data.frame(Lp.mean.chain=probit_llk_ratio,
#                             x=llk_plot_grid),
#             aes(x=x, y=Lp.mean.chain))+
#   geom_hline(yintercept=1, linetype="dashed", color = "red")+
#   theme_bw(base_size = 12) +
#   labs(title="likelihood ratio using different p for probit data",
#         x="p",
#        y="lilelihood ratio")+scale_x_continuous(breaks=seq(0,8,1))
#

# plot_b <- ggplot() +
#   geom_line(data=data.frame(Lp.mean.chain=logit_llk_ratio,
#                             x=llk_plot_grid),
#             aes(x=x, y=Lp.mean.chain))+
#   geom_hline(yintercept=1, linetype="dashed", color = "red")+
#   theme_bw(base_size = 12) +
#   labs(title="likelihood ratio using different p for logit data",
#        x="p",
#        y="lilelihood ratio")+scale_x_continuous(breaks=seq(0,8,1))
#


# plot_c<-
# ggplot() +
#   geom_line(data=data.frame(Lp.mean.chain=p_1_llk_ratio,
#                             x=llk_plot_grid),
#             aes(x=x, y=Lp.mean.chain))+
#   geom_hline(yintercept=1, linetype="dashed", color = "red")+
#   theme_bw(base_size = 12) +
#   labs(title="likelihood ratio using different p for p=1 data",
#        x="p",
#        y="lilelihood ratio") +scale_x_continuous(breaks=seq(0,8,1))
# plot_d<-
#   ggplot() +
#   geom_line(data=data.frame(Lp.mean.chain=p_2_llk_ratio,
#                             x=llk_plot_grid),
#             aes(x=x, y=Lp.mean.chain))+
#   geom_hline(yintercept=1, linetype="dashed", color = "red")+
#   theme_bw(base_size = 12) +
#   labs(title="likelihood ratio using different p for p=2 data",
#        x="p",
#        y="lilelihood ratio") +scale_x_continuous(breaks=seq(0,8,1))
#

#
#
P_logit$Lp.mean
var(P_logit$Lp.mean.chain)
norm(as.matrix(P_logit$Beta.mean-true_theta,2))
pi_logit_predict <- exp(X%*%P_logit$Beta.mean)/(1+exp(X%*%P_logit$Beta.mean))
sum((y.logit-pi_logit_predict)^2)

P_probit$Lp.mean
var(P_probit$Lp.mean.chain)
norm(as.matrix(P_probit$Beta.mean-true_theta,2))
pi_probit_predict <- pgnorm(X%*%P_probit$Beta.mean,alpha = p_scale(P_probit$Lp.mean),beta = P_probit$Lp.mean)
sum((y.probit-pi_probit_predict)^2)


P_0.5$Lp.mean
var(P_0.5$Lp.mean.chain)
norm(as.matrix(P_0.5$Beta.mean-true_theta,2))
pi_probit_predict <- pgnorm(X%*%P_0.5$Beta.mean,alpha = p_scale(P_0.5$Lp.mean),beta = P_0.5$Lp.mean)
sum((y.probit-pi_probit_predict)^2)

P_1$Lp.mean
var(P_1$Lp.mean.chain)
norm(as.matrix(P_1$Beta.mean-true_theta,2))
pi_probit_predict <- pgnorm(X%*%P_1$Beta.mean,alpha = p_scale(P_1$Lp.mean),beta = P_1$Lp.mean)
sum((y.probit-pi_probit_predict)^2)

P_1.5$Lp.mean
var(P_1.5$Lp.mean.chain)
norm(as.matrix(P_1.5$Beta.mean-true_theta,2))
pi_probit_predict <- pgnorm(X%*%P_1.5$Beta.mean,alpha = p_scale(P_1.5$Lp.mean),beta = P_1.5$Lp.mean)
sum((y.probit-pi_probit_predict)^2)

P_3$Lp.mean
var(P_3$Lp.mean.chain)
norm(as.matrix(P_3$Beta.mean-true_theta,2))
pi_probit_predict <- pgnorm(X%*%P_3$Beta.mean,alpha = p_scale(P_3$Lp.mean),beta = P_3$Lp.mean)
sum((y.probit-pi_probit_predict)^2)

P_4$Lp.mean
var(P_4$Lp.mean.chain)
norm(as.matrix(P_4$Beta.mean-true_theta,2))
pi_probit_predict <- pgnorm(X%*%P_4$Beta.mean,alpha = p_scale(P_4$Lp.mean),beta = P_4$Lp.mean)
sum((y.probit-pi_probit_predict)^2)

P_5$Lp.mean
var(P_5$Lp.mean.chain)
norm(as.matrix(P_5$Beta.mean-true_theta,2))
pi_probit_predict <- pgnorm(X%*%P_5$Beta.mean,alpha = p_scale(P_5$Lp.mean),beta = P_5$Lp.mean)
sum((y.probit-pi_probit_predict)^2)


P_8$Lp.mean
var(P_8$Lp.mean.chain)
norm(as.matrix(P_8$Beta.mean-true_theta,2))
pi_probit_predict <- pgnorm(X%*%P_8$Beta.mean,alpha = p_scale(P_8$Lp.mean),beta = P_8$Lp.mean)
sum((y.probit-pi_probit_predict)^2)

plot(density(tail(P_8$Lp.mean.chain,1000)))

#####Fix p using same data

#Logit model
logit_probit <- glm(y.probit~X[,-1],family = binomial(link = "logit"))
logit_logit <- glm(y.logit~X[,-1],family = binomial(link = "logit"))
logit0.5 <- glm(y_p0.5~X[,-1],family = binomial(link = "logit"))
logit1 <- glm(y_p1~X[,-1],family = binomial(link = "logit"))
logit1.5 <- glm(y_p1.5~X[,-1],family = binomial(link = "logit"))
logit3 <- glm(y_p3~X[,-1],family = binomial(link = "logit"))
logit4 <- glm(y_p4~X[,-1],family = binomial(link = "logit"))
logit5 <- glm(y_p5~X[,-1],family = binomial(link = "logit"))
logit8 <- glm(y_p8~X[,-1],family = binomial(link = "logit"))
#Estimation probit model for all scenarios
probit_probit <- glm(y.probit~X[,-1],family = binomial(link = "probit"))
probit_logit <- glm(y.logit~X[,-1],family = binomial(link = "probit"))
probit0.5 <- glm(y_p0.5~X[,-1],family = binomial(link = "probit"))
probit1 <- glm(y_p1~X[,-1],family = binomial(link = "probit"))
probit1.5 <- glm(y_p1.5~X[,-1],family = binomial(link = "probit"))
probit3 <- glm(y_p3~X[,-1],family = binomial(link = "probit"))
probit4 <- glm(y_p4~X[,-1],family = binomial(link = "probit"))
probit5 <- glm(y_p5~X[,-1],family = binomial(link = "probit"))
probit8 <- glm(y_p8~X[,-1],family = binomial(link = "probit"))
#Estimation cloglog model for all scenarios
cloglog_probit <- glm(y.probit~X[,-1],family = binomial(link = "cloglog"))
cloglog_logit <- glm(y.logit~X[,-1],family = binomial(link = "cloglog"))
cloglog_0.5 <- glm(y_p0.5~X[,-1],family = binomial(link = "cloglog"))
cloglog_1 <- glm(y_p1~X[,-1],family = binomial(link = "cloglog"))
cloglog_1.5 <- glm(y_p1.5~X[,-1],family = binomial(link = "cloglog"))
cloglog_3 <- glm(y_p3~X[,-1],family = binomial(link = "cloglog"))
cloglog_4 <- glm(y_p4~X[,-1],family = binomial(link = "cloglog"))
cloglog_5 <- glm(y_p5~X[,-1],family = binomial(link = "cloglog"))
cloglog_8 <- glm(y_p8~X[,-1],family = binomial(link = "cloglog"))
#Estimation p-probit model for all scenarios.
p_probit_probit <- Lp_gibbssampler(N_sim=1000,
                burn_in=900, X=X, y=y.probit,true_theta = true_theta,Lp=2,
                initial_theda = MLE_probit_coefficients)
p_probit_logit <- Lp_gibbssampler(N_sim=1000,
                                   burn_in=900, X=X, y=y.logit,true_theta = true_theta,Lp=1,
                                   initial_theda = MLE_logit$coefficients)
p_probit_0.5 <- Lp_gibbssampler(N_sim=1000,
                                   burn_in=900, X=X, y=y_p0.5,true_theta = true_theta,Lp=0.5,
                                   initial_theda = MLE_logit$coefficients)
p_probit_1 <- Lp_gibbssampler(N_sim=1000,
                                   burn_in=900, X=X, y=y_p1,true_theta = true_theta,Lp=1,
                                   initial_theda = MLE_logit$coefficients)
p_probit_1.5 <- Lp_gibbssampler(N_sim=1000,
                                   burn_in=900, X=X, y=y_p1.5,true_theta = true_theta,Lp=1.5,
                                   initial_theda = MLE_logit$coefficients)
p_probit_3 <- Lp_gibbssampler(N_sim=1000,
                                   burn_in=900, X=X, y=y_p3,true_theta = true_theta,Lp=3,
                                   initial_theda = MLE_logit$coefficients)
p_probit_4 <- Lp_gibbssampler(N_sim=1000,
                                   burn_in=900, X=X, y=y_p4,true_theta = true_theta,Lp=4,
                                   initial_theda = MLE_logit$coefficients)
p_probit_5 <- Lp_gibbssampler(N_sim=1000,
                                   burn_in=900, X=X, y=y_p5,true_theta = true_theta,Lp=5,
                                   initial_theda = MLE_logit$coefficients)
p_probit_8 <- Lp_gibbssampler(N_sim=1000,
                                   burn_in=900, X=X, y=y_p8,true_theta = true_theta,Lp=8,
                                   initial_theda = MLE_logit$coefficients)
#Calculate the fitted value of probability for the p-probit model of all scenarios
p_probit_probit_pred <- pgnorm(X%*%p_probit_probit$post.mean,
                               alpha = p_scale(2),beta = 2)
p_probit_logit_pred <- pgnorm(X%*%p_probit_logit$post.mean,
                               alpha = p_scale(1),beta = 1)
p_probit_0.5_pred <- pgnorm(X%*%p_probit_0.5$post.mean,
                              alpha = p_scale(0.5),beta = 0.5)
p_probit_1_pred <- pgnorm(X%*%p_probit_1$post.mean,
                              alpha = p_scale(1),beta = 1)
p_probit_1.5_pred <- pgnorm(X%*%p_probit_1.5$post.mean,
                              alpha = p_scale(1.5),beta = 1.5)
p_probit_3_pred <- pgnorm(X%*%p_probit_3$post.mean,
                              alpha = p_scale(3),beta = 3)
p_probit_4_pred <- pgnorm(X%*%p_probit_4$post.mean,
                              alpha = p_scale(4),beta = 4)
p_probit_5_pred <- pgnorm(X%*%p_probit_5$post.mean,
                              alpha = p_scale(5),beta = 5)
p_probit_8_pred <- pgnorm(X%*%p_probit_8$post.mean,
                              alpha = p_scale(8),beta = 8)

#Calculate the absolute error of all scenarios.

#logit data
sum(abs(pi.logit-p_probit_logit_pred))
sum(abs(pi.logit-probit_logit$fitted.values))
sum(abs(pi.logit-logit_logit$fitted.values))
sum(abs(pi.logit-cloglog_logit$fitted.values))
#probit data
sum(abs(pi.probit-p_probit_probit_pred))
sum(abs(pi.probit-probit_probit$fitted.values))
sum(abs(pi.probit-logit_probit$fitted.values))
sum(abs(pi.probit-cloglog_probit$fitted.values))

#p=0.5
sum(abs(pi_0.5-p_probit_0.5_pred))
sum(abs(pi_0.5-probit0.5$fitted.values))
sum(abs(pi_0.5-logit0.5$fitted.values))
sum(abs(pi_0.5-cloglog_0.5$fitted.values))

#p=1
sum(abs(pi_1-p_probit_1_pred))
sum(abs(pi_1-probit1$fitted.values))
sum(abs(pi_1-logit1$fitted.values))
sum(abs(pi_1-cloglog_1$fitted.values))

#p=1.5
sum(abs(pi_1.5-p_probit_1.5_pred))
sum(abs(pi_1.5-probit1.5$fitted.values))
sum(abs(pi_1.5-logit1.5$fitted.values))
sum(abs(pi_1.5-cloglog_1.5$fitted.values))

#p=3
sum(abs(pi_3-p_probit_3_pred))
sum(abs(pi_3-probit3$fitted.values))
sum(abs(pi_3-logit3$fitted.values))
sum(abs(pi_3-cloglog_3$fitted.values))

#p=4
sum(abs(pi_4-p_probit_4_pred))
sum(abs(pi_4-probit4$fitted.values))
sum(abs(pi_4-logit4$fitted.values))
sum(abs(pi_4-cloglog_4$fitted.values))

#p=5
sum(abs(pi_5-p_probit_5_pred))
sum(abs(pi_5-probit5$fitted.values))
sum(abs(pi_5-logit5$fitted.values))
sum(abs(pi_5-cloglog_5$fitted.values))

#p=8
sum(abs(pi_8-p_probit_8_pred))
sum(abs(pi_8-probit8$fitted.values))
sum(abs(pi_8-logit8$fitted.values))
sum(abs(pi_8-cloglog_8$fitted.values))

#Calculate the likelihood ratio

#logit data
llk(X,true_theta,1,y.logit)/llk(X,p_probit_logit$post.mean,1,y.logit)
llk(X,true_theta,2,y.logit)/logLik(probit_logit)
llk_logit(X,true_theta,y.logit)/logLik(logit_logit)
llk_clog(X,true_theta,y.logit)/logLik(cloglog_logit)
#probit data
llk(X,true_theta,2,y.probit)/llk(X,p_probit_probit$post.mean,2,y.probit)
llk(X,true_theta,2,y.probit)/logLik(probit_probit)
llk_logit(X,true_theta,y.probit)/logLik(logit_probit)
llk_clog(X,true_theta,y.probit)/logLik(cloglog_probit)

#p=0.5
llk(X,true_theta,0.5,y_p0.5)/llk(X,p_probit_0.5$post.mean,0.5,y_p0.5)
llk(X,true_theta,2,y_p0.5)/logLik(probit0.5)
llk_logit(X,true_theta,y_p0.5)/logLik(logit0.5)
llk_clog(X,true_theta,y_p0.5)/logLik(cloglog_0.5)

#p=1
llk(X,true_theta,1,y_p1)/llk(X,p_probit_1$post.mean,1,y_p1)
llk(X,true_theta,2,y_p1)/logLik(probit1)
llk_logit(X,true_theta,y_p1)/logLik(logit1)
llk_clog(X,true_theta,y_p1)/logLik(cloglog_1)

#p=1.5
llk(X,true_theta,1.5,y_p1.5)/llk(X,p_probit_1.5$post.mean,1.5,y_p1.5)
llk(X,true_theta,2,y_p1.5)/logLik(probit1.5)
llk_logit(X,true_theta,y_p1.5)/logLik(logit1.5)
llk_clog(X,true_theta,y_p1.5)/logLik(cloglog_1.5)

#p=3
llk(X,true_theta,3,y_p3)/llk(X,p_probit_3$post.mean,3,y_p3)
llk(X,true_theta,2,y_p3)/logLik(probit3)
llk_logit(X,true_theta,y_p3)/logLik(logit3)
llk_clog(X,true_theta,y_p3)/logLik(cloglog_3)

#p=4
llk(X,true_theta,4,y_p4)/llk(X,p_probit_4$post.mean,4,y_p4)
llk(X,true_theta,2,y_p4)/logLik(probit4)
llk_logit(X,true_theta,y_p4)/logLik(logit4)
llk_clog(X,true_theta,y_p4)/logLik(cloglog_4)


#p=5
llk(X,true_theta,5,y_p5)/llk(X,p_probit_5$post.mean,5,y_p5)
llk(X,true_theta,2,y_p5)/logLik(probit5)
llk_logit(X,true_theta,y_p5)/logLik(logit5)
llk_clog(X,true_theta,y_p5)/logLik(cloglog_5)


#p=8
llk(X,true_theta,8,y_p8)/llk(X,p_probit_8$post.mean,8,y_p8)
llk(X,true_theta,2,y_p8)/logLik(probit8)
llk_logit(X,true_theta,y_p8)/logLik(logit8)
llk_clog(X,true_theta,y_p8)/logLik(cloglog_8)

#MSE

#logit data
sum((y.logit-p_probit_logit_pred)^2)
sum((y.logit-probit_logit$fitted.values)^2)
sum((y.logit-logit_logit$fitted.values)^2)
sum((y.logit-cloglog_logit$fitted.values)^2)

#probit data

sum((y.probit-p_probit_probit_pred)^2)
sum((y.probit-probit_probit$fitted.values)^2)
sum((y.probit-logit_probit$fitted.values)^2)
sum((y.probit-cloglog_probit$fitted.values)^2)

#p=0.5
sum((y_p0.5-p_probit_0.5_pred)^2)
sum((y_p0.5-probit0.5$fitted.values)^2)
sum((y_p0.5-logit0.5$fitted.values)^2)
sum((y_p0.5-cloglog_0.5$fitted.values)^2)

#p=1
sum((y_p1-p_probit_1_pred)^2)
sum((y_p1-probit1$fitted.values)^2)
sum((y_p1-logit1$fitted.values)^2)
sum((y_p1-cloglog_1$fitted.values)^2)

#p=1.5
sum((y_p1.5-p_probit_1.5_pred)^2)
sum((y_p1.5-probit1.5$fitted.values)^2)
sum((y_p1.5-logit1.5$fitted.values)^2)
sum((y_p1.5-cloglog_1.5$fitted.values)^2)

#p=3
sum((y_p3-p_probit_3_pred)^2)
sum((y_p3-probit3$fitted.values)^2)
sum((y_p3-logit3$fitted.values)^2)
sum((y_p3-cloglog_3$fitted.values)^2)

#p=4
sum((y_p4-p_probit_4_pred)^2)
sum((y_p4-probit4$fitted.values)^2)
sum((y_p4-logit4$fitted.values)^2)
sum((y_p4-cloglog_4$fitted.values)^2)

#p=5
sum((y_p5-p_probit_5_pred)^2)
sum((y_p5-probit5$fitted.values)^2)
sum((y_p5-logit5$fitted.values)^2)
sum((y_p5-cloglog_5$fitted.values)^2)

#p=8
sum((y_p8-p_probit_8_pred)^2)
sum((y_p8-probit8$fitted.values)^2)
sum((y_p8-logit8$fitted.values)^2)
sum((y_p8-cloglog_8$fitted.values)^2)




#logit data
sum((pi.logit-p_probit_logit_pred)^2)
sum((pi.logit-probit_logit$fitted.values)^2)
sum((pi.logit-logit_logit$fitted.values)^2)
sum((pi.logit-cloglog_logit$fitted.values)^2)

#probit data

sum((pi.probit-p_probit_probit_pred)^2)
sum((pi.probit-probit_probit$fitted.values)^2)
sum((pi.probit-logit_probit$fitted.values)^2)
sum((pi.probit-cloglog_probit$fitted.values)^2)

#p=0.5
sum((pi_0.5-p_probit_0.5_pred)^2)
sum((pi_0.5-probit0.5$fitted.values)^2)
sum((pi_0.5-logit0.5$fitted.values)^2)
sum((pi_0.5-cloglog_0.5$fitted.values)^2)

#p=1
sum((pi_1-p_probit_1_pred)^2)
sum((pi_1-probit1$fitted.values)^2)
sum((pi_1-logit1$fitted.values)^2)
sum((pi_1-cloglog_1$fitted.values)^2)

#p=1.5
sum((pi_1.5-p_probit_1.5_pred)^2)
sum((pi_1.5-probit1.5$fitted.values)^2)
sum((pi_1.5-logit1.5$fitted.values)^2)
sum((pi_1.5-cloglog_1.5$fitted.values)^2)

#p=3
sum((pi_3-p_probit_3_pred)^2)
sum((pi_3-probit3$fitted.values)^2)
sum((pi_3-logit3$fitted.values)^2)
sum((pi_3-cloglog_3$fitted.values)^2)

#p=4
sum((pi_4-p_probit_4_pred)^2)
sum((pi_4-probit4$fitted.values)^2)
sum((pi_4-logit4$fitted.values)^2)
sum((pi_4-cloglog_4$fitted.values)^2)

#p=5
sum((pi_5-p_probit_5_pred)^2)
sum((pi_5-probit5$fitted.values)^2)
sum((pi_5-logit5$fitted.values)^2)
sum((pi_5-cloglog_5$fitted.values)^2)

#p=8
sum((pi_8-p_probit_8_pred)^2)
sum((pi_8-probit8$fitted.values)^2)
sum((pi_8-logit8$fitted.values)^2)
sum((pi_8-cloglog_8$fitted.values)^2)

















