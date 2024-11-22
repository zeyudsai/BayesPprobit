getwd()
rm(list = ls())
dev.off()



setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("Functions.R")

# setwd("/Users/dingzeyu/Downloads/Fraud Detection Case/")
creditdata <- read.csv("Documents/Fraud Detection Case/creditcard.csv", sep = ",")

## Remove rows that do not have target variable values
creditdata <- creditdata[!(is.na(creditdata$Class)),]

creditdata$Class <- factor(creditdata$Class)
creditdata <- creditdata[,-1]
creditdata$Class <- as.numeric(creditdata$Class)
creditdata$Class[creditdata$Class!=2] <- 0
creditdata$Class[creditdata$Class==2] <- 1
D <- ncol(creditdata)
N <- nrow(creditdata)
X <- creditdata[,1:29]
onevector <- rep(1,N)
X <- as.matrix(cbind(onevector,X))
y <- creditdata$Class
mle_result <- glm(Class~.,data=creditdata,family = binomial(link = "logit"))
mle_result$coefficients
mle_probit <- glm(Class~.,data=creditdata,family = binomial(link = "probit"))
MLE_coefficients <- mle_result$coefficients
probit_coefficients <- mle_probit$coefficients
N_sim <- 1000
burn_in <- N_sim*0.9
credit_p2 <- Lp_gibbssampler(N_sim=N_sim,
                             burn_in=burn_in, X=X, y=y,true_theta = probit_coefficients,Lp=2,initial_theda = probit_coefficients)
credit_p1 <- Lp_gibbssampler(N_sim=N_sim,
                             burn_in=burn_in, X=X, y=y,true_theta = probit_coefficients,Lp=1,initial_theda = MLE_coefficients)

credit_p2$post.mean
mle_probit$coefficients
mle_result$coefficients

sample3 <- multi_chain(N_sim=N_sim,burn_in=burn_in,
                       X=X, y=y,initial_theda=probit_coefficients,
                       true_theta=probit_coefficients,Lp=2,M_iter=100,range=c(0.1,5),
                       step=0.05,times = 1)
sketch_size <- D^2
frac <-seq(0.00017211348,0.05,0.00086056742)

sample3ci<-apply(sample3$Beta_chain[[1]], MARGIN =2, FUN=ci)
sample3_ci_low <- c()
sample3_ci_high <- c()
for (i in 1:D) {
  sample3_ci_low[i] <- sample3ci[[i]]$CI_low
  sample3_ci_high[i] <- sample3ci[[i]]$CI_high
  
}
sample1_ci <-apply(credit_p1$chain, MARGIN =2, FUN=ci)
sample1_ci_low <- c()
sample1_ci_high <- c()
for (i in 1:D) {
  sample1_ci_low[i] <- sample1_ci[[i]]$CI_low
  sample1_ci_high[i] <- sample1_ci[[i]]$CI_high
  
}
sample2_ci <-apply(credit_p2$chain, MARGIN =2, FUN=ci)
sample2_ci_low <- c()
sample2_ci_high <- c()
for (i in 1:D) {
  sample2_ci_low[i] <- sample2_ci[[i]]$CI_low
  sample2_ci_high[i] <- sample2_ci[[i]]$CI_high
  
}
summary(mle_result)
result <-summary(mle_probit)
probit_confidence_low <- result$coefficients[,1]-result$coefficients[,2]*qt(0.95,df=N-2)
probit_confidence_high <- result$coefficients[,1]+result$coefficients[,2]*qt(0.95,df=N-2)
probit_confidence_high <- probit_confidence_high[-1]
probit_confidence_low <- probit_confidence_low[-1]
sample1_ci_low <- sample1_ci_low[-1]
sample1_ci_high <- sample1_ci_high[-1]
sample2_ci_low <- sample2_ci_low[-1]
sample2_ci_high <- sample2_ci_high[-1]
sample3_ci_low <- sample3_ci_low[-1]
sample3_ci_high <- sample3_ci_high[-1]
result <-summary(mle_result)

logit_confidence_low <- result$coefficients[,1]-result$coefficients[,2]*qt(0.95,df=N-2)
logit_confidence_high <- result$coefficients[,1]+result$coefficients[,2]*qt(0.95,df=N-2)
logit_confidence_high <- logit_confidence_high[-1]
logit_confidence_low <- logit_confidence_low[-1]

ci_plot1 <- data.frame("variable"=c(1:15),
                       "sample3_ci_low"=c(sample1_ci_low[1:15],sample2_ci_low[1:15],
                                          logit_confidence_low[1:15],probit_confidence_low[1:15]),
                       "sample3_ci_high"=c(sample1_ci_high[1:15],sample2_ci_high[1:15],
                                           logit_confidence_high[1:15],probit_confidence_high[1:15]),
                       "post_mean"=c(credit_p1$post.mean[2:16],credit_p2$post.mean[2:16],
                                     mle_result$coefficients[2:16],mle_probit$coefficients[2:16]),
                       "Method"=c(rep("Bayesian Credible Interval (P=1)",15),
                                  rep("Bayesian Credible Interval (P=2)",15),
                                  rep("Logit Confidence Interval",15),
                                  rep("Probit Confidence Interval",15)))

ci_plot2 <- data.frame("variable"=c(16:29),
                       "sample3_ci_low"=c(sample1_ci_low[16:29],sample2_ci_low[16:29],
                                          logit_confidence_low[16:29],probit_confidence_low[16:29]),
                       "sample3_ci_high"=c(sample1_ci_high[16:29],sample2_ci_high[16:29],
                                           logit_confidence_high[16:29],probit_confidence_high[16:29]),
                       "post_mean"=c(credit_p1$post.mean[17:30],credit_p2$post.mean[17:30],
                                     mle_result$coefficients[17:30],mle_probit$coefficients[17:30]),
                       "Method"=c(rep("Bayesian Credible Interval (P=1)",14),
                                  rep("Bayesian Credible Interval (P=2)",14),
                                  rep("Logit Confidence Interval",14),
                                  rep("Probit Confidence Interval",14)))



jpeg("credit_ci_1(version2).jpeg", units="in", width=8, height=5, res=300)
ggplot(ci_plot1, aes(variable, post_mean, colour=Method)) +
  geom_point(position=position_dodge(1)) + 
  geom_errorbar(aes(ymin=sample3_ci_low, ymax=sample3_ci_high), width=0.5, position=position_dodge(1)) + 
  theme_bw()+
  labs(x=TeX("$\\hat{\\beta}\\ index"), y=TeX("$\\hat{\\beta}"), colour="Method", fill="Method")+ 
  # ggtitle("Confidence/Credible Interval of Coefficients \nCredit Card Data") +
  theme(legend.position = "none")
dev.off()

jpeg("credit_ci_2(version2).jpeg", units="in", width=8, height=5, res=300)
ggplot(ci_plot2, aes(variable, post_mean, colour=Method)) +
  geom_point(position=position_dodge(1)) + 
  geom_errorbar(aes(ymin=sample3_ci_low, ymax=sample3_ci_high), width=1, position=position_dodge(1)) + 
  theme_bw()+
  labs(x=TeX("$\\hat{\\beta}\\ index"), y=TeX("$\\hat{\\beta}"), colour="Method", fill="Method")+
  # ggtitle("Confidence/Credible Interval of Coefficients \nCredit Card Data") +
  theme(legend.position = "none")

dev.off()

jpeg("credit_p_estimation(version2).jpeg", units="in", width=8, height=5, res=300)
ggplot(data=data.frame(Lp.mean.chain=sample3$Lp.mean.chain[1:N_sim],x=1:length(sample3$Lp.mean.chain)), aes(x=x, y=Lp.mean.chain)) +
  geom_line()+
  annotate("text", x=150, y=2.8,size=6, label= TeX("$\\hat{p} = 2.053", output='character'),parse=TRUE) +
  # annotate("text", x=130, y=4.8,size=6, label= "Estimation of p for credit card data") +
  theme_bw(base_size = 12) +
  ylim(1,5)+
  labs(
     title="Estimation of p for credit card data", 
    # subtitle="Estimated p=2.05", 
    x="MCMC iteration",
       y="p")
dev.off()

plot(1:(D-1),probit_coefficients[-1],type = "l",ylim=c(-1,1),main = "comparison of coefficients",xlab = "Beta-index",
     ylab = "parameters")
Msubtitle <- "credit card data"
mtext(side = 3, line = 0.4, Msubtitle)
legend("topright",legend=c("MLE Probit","Bayesian-1-probit","Bayesian-2-probit"), col=c(1,2,3), 
       bty = 'n', lwd = 2, inset = c(0, 0), lty = 1, cex = 0.73)
lines(credit_p1$post.mean[-1],col=2)
lines(credit_p2$post.mean[-1],col=3)

plot(1:(D-1),probit_coefficients[-1],type = "l",ylim = c(-2,2),main = "comparison of coefficients",xlab = "Beta-index",
     ylab = "parameters")
Msubtitle <- "credit card data"
mtext(side = 3, line = 0.4, Msubtitle)
legend("topright",legend=c("MLE Probit","Bayesian-P-2.15(estimated)-probit"), col=c(1,2,3), 
       bty = 'n', lwd = 2, inset = c(0, 0), lty = 1, cex = 0.73)
lines(sample3$Beta.mean[-1],col=2)

plot(1:(D-1),MLE_coefficients[-1],type = "l",ylim = c(-2,2),main = "comparison of coefficients",xlab = "Beta-index",
     ylab = "parameters")
Msubtitle <- "credit card data"
mtext(side = 3, line = 0.4, Msubtitle)
legend("topright",legend=c("MLE Logit","Bayesian-1-probit","Bayesian-2-probit"), col=c(1,2,3), 
       bty = 'n', lwd = 2, inset = c(0, 0), lty = 1, cex = 0.73)
lines(credit_p1$post.mean[-1],col=2)
lines(credit_p2$post.mean[-1],col=3)

plot(1:(D-1),MLE_coefficients[-1],type = "l",ylim = c(-2,2),main = "comparison of coefficients",xlab = "Beta-index",
     ylab = "parameters")
Msubtitle <- "credit card data"
mtext(side = 3, line = 0.4, Msubtitle)
legend("topright",legend=c("MLE Logit","Bayesian-P-2.15(estimated)-probit"), col=c(1,2,3), 
       bty = 'n', lwd = 2, inset = c(0, 0), lty = 1, cex = 0.73)
lines(sample3$Beta.mean[-1],col=2)


##############p=2######################
for (i in 1:length(frac)) {
  for (j in 1:5) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Lp_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=2)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y[sample_index],
                             initial_theda=probit_coefficients,
                             true_theta=probit_coefficients,Lp=2)
    assign(paste("credit_p_2_rootl2_frac",frac[i],"chain_",j,sep = ""), Model)
  }
}

for (i in 1:length(frac)) {
  for (j in 1:5) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Uniform_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=2)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y[sample_index],
                             initial_theda=probit_coefficients,
                             true_theta=probit_coefficients,Lp=2)
    assign(paste("credit_p_2_uniform_frac",frac[i],"chain_",j,sep = ""), Model)
  }
}
for (i in 1:length(frac)) {
  for (j in 1:5) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Compute_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=2)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y[sample_index],
                             initial_theda=probit_coefficients,
                             true_theta=probit_coefficients,Lp=2)
    assign(paste("credit_p_2_2probit_frac",frac[i],"chain_",j,sep = ""), Model)
  }
}
for (i in 1:length(frac)) {
  for (j in 1:5) {
    coreset_size<- N*frac[i]
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    scores<- one_shot_coreset(sketch_size = sketch_size, coreset_size = coreset_size, X=X, Lp_max = 2.5)
    sample_index <- sample(1:N,size = coreset_size,
                           replace = FALSE, prob = scores)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y[sample_index],
                             initial_theda=probit_coefficients,
                             true_theta=probit_coefficients,Lp=2)
    assign(paste("credit_p_2_one_shot_frac",frac[i],"chain_",j,sep = ""), Model)
  }
}


##############p=1######################
for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Lp_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=1)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y[sample_index],
                             initial_theda=probit_coefficients,
                             true_theta=probit_coefficients,Lp=1)
    assign(paste("credit_p_1_rootl2_frac",frac[i],"chain_",j,sep = ""), Model)
  }
}

for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Uniform_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=1)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y[sample_index],
                             initial_theda=probit_coefficients,
                             true_theta=probit_coefficients,Lp=1)
    assign(paste("credit_p_1_uniform_frac",frac[i],"chain_",j,sep = ""), Model)
  }
}


for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Compute_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=1)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y[sample_index],
                             initial_theda=probit_coefficients,
                             true_theta=probit_coefficients,Lp=1)
    assign(paste("credit_p_1_1probit_frac",frac[i],"chain_",j,sep = ""), Model)
  }
}


for (i in 1:length(frac)) {
  for (j in 1:10) {
    coreset_size<- N*frac[i]
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    scores<- one_shot_coreset(sketch_size = sketch_size, coreset_size = coreset_size, X=X, Lp_max = 1.5)
    sample_index <- sample(1:N,size = coreset_size,
                           replace = FALSE, prob = scores)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y[sample_index],
                             initial_theda=probit_coefficients,
                             true_theta=probit_coefficients,Lp=1)
    assign(paste("credit_p_1_one_shot_frac",frac[i],"chain_",j,sep = ""), Model)
  }
}
############plot of p=2

credit_p_2_one_shot_mean_diff_norm_median <- c()
credit_p_2_one_shot_mean_diff_norm_max <- c()
credit_p_2_one_shot_mean_diff_norm_min <- c()
credit_p_2_one_shot_mean_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    b <- norm(get(paste("credit_p_2_one_shot_frac",frac[i],"chain_",j,sep = ""))$post.mean-credit_p2$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p_2_one_shot_mean_diff_norm_median <- c(credit_p_2_one_shot_mean_diff_norm_median,median_result)
  credit_p_2_one_shot_mean_diff_norm_max <- c(credit_p_2_one_shot_mean_diff_norm_max,maxresult)
  credit_p_2_one_shot_mean_diff_norm_min <- c(credit_p_2_one_shot_mean_diff_norm_min,minresult)
  credit_p_2_one_shot_mean_diff_norm_average <- c(credit_p_2_one_shot_mean_diff_norm_average,meanresult)
}
credit_p_2_rootl2_mean_diff_norm_median <- c()
credit_p_2_rootl2_mean_diff_norm_max <- c()
credit_p_2_rootl2_mean_diff_norm_min <- c()
credit_p_2_rootl2_mean_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    b <- norm(get(paste("credit_p_2_rootl2_frac",frac[i],"chain_",j,sep = ""))$post.mean-credit_p2$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p_2_rootl2_mean_diff_norm_median <- c(credit_p_2_rootl2_mean_diff_norm_median,median_result)
  credit_p_2_rootl2_mean_diff_norm_max <- c(credit_p_2_rootl2_mean_diff_norm_max,maxresult)
  credit_p_2_rootl2_mean_diff_norm_min <- c(credit_p_2_rootl2_mean_diff_norm_min,minresult)
  credit_p_2_rootl2_mean_diff_norm_average <- c(credit_p_2_rootl2_mean_diff_norm_average,meanresult)
}


credit_p_2_2probit_mean_diff_norm_median <- c()
credit_p_2_2probit_mean_diff_norm_max <- c()
credit_p_2_2probit_mean_diff_norm_min <- c()
credit_p_2_2probit_mean_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    b <- norm(get(paste("credit_p_2_2probit_frac",frac[i],"chain_",j,sep = ""))$post.mean-credit_p2$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p_2_2probit_mean_diff_norm_median <- c(credit_p_2_2probit_mean_diff_norm_median,median_result)
  credit_p_2_2probit_mean_diff_norm_max <- c(credit_p_2_2probit_mean_diff_norm_max,maxresult)
  credit_p_2_2probit_mean_diff_norm_min <- c(credit_p_2_2probit_mean_diff_norm_min,minresult)
  credit_p_2_2probit_mean_diff_norm_average <- c(credit_p_2_2probit_mean_diff_norm_average,meanresult)
}


credit_p_2_uniform_mean_diff_norm_median <- c()
credit_p_2_uniform_mean_diff_norm_max <- c()
credit_p_2_uniform_mean_diff_norm_min <- c()
credit_p_2_uniform_mean_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    b <- norm(get(paste("credit_p_2_uniform_frac",frac[i],"chain_",j,sep = ""))$post.mean-credit_p2$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p_2_uniform_mean_diff_norm_median <- c(credit_p_2_uniform_mean_diff_norm_median,median_result)
  credit_p_2_uniform_mean_diff_norm_max <- c(credit_p_2_uniform_mean_diff_norm_max,maxresult)
  credit_p_2_uniform_mean_diff_norm_min <- c(credit_p_2_uniform_mean_diff_norm_min,minresult)
  credit_p_2_uniform_mean_diff_norm_average <- c(credit_p_2_uniform_mean_diff_norm_average,meanresult)
}

p2_meandiff_credit_data <- data.frame(mean=c(credit_p_2_2probit_mean_diff_norm_average,
                                             credit_p_2_one_shot_mean_diff_norm_average,
                                                credit_p_2_uniform_mean_diff_norm_average),
                                         median=c(credit_p_2_2probit_mean_diff_norm_median,
                                                  credit_p_2_one_shot_mean_diff_norm_median,
                                                  credit_p_2_uniform_mean_diff_norm_median),
                                         max=c(credit_p_2_2probit_mean_diff_norm_max,
                                               credit_p_2_one_shot_mean_diff_norm_max,
                                               credit_p_2_uniform_mean_diff_norm_max),
                                         min=c(credit_p_2_2probit_mean_diff_norm_min,
                                               credit_p_2_one_shot_mean_diff_norm_min,
                                               credit_p_2_uniform_mean_diff_norm_min),
                                         label=c(rep("2-probit",length(frac)),rep("one-shot",length(frac)),rep("uniform",length(frac))),
                                         frac=c(rep(frac,3)))
p2_meandiff_credit_data
jpeg("credit_p2_mean.jpeg", units="in", width=8, height=5, res=300)

ggplot(p2_meandiff_credit_data, aes(frac, median, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior mean", colour="method", fill="method") +
  theme_bw()+coord_cartesian(ylim=c(0, 100))+
  ggtitle("Credit card data (p=2)")

dev.off()

############lilelihood ratio###############


credit_p2_2probit_llk_ratio_median <- c()
credit_p2_2probit_llk_ratio_max <- c()
credit_p2_2probit_llk_ratio_min <- c()
credit_p2_2probit_llk_ratio_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    # b <- sum(llk(X=a,get(paste("credit_p_2_2probit_frac",frac[i],"chain_",j,sep = ""))$post.mean,2))/sum(llk(X=a, credit_p2$post.mean,2))
    
    b<- median(abs(Bayes_factor(X, credit_p2$post.mean,2,y)-
                     Bayes_factor(X,get(paste("credit_p_2_2probit_frac",frac[i],"chain_",j,sep = ""))$post.mean,2,y)))
    
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p2_2probit_llk_ratio_median <- c(credit_p2_2probit_llk_ratio_median,median_result)
  credit_p2_2probit_llk_ratio_max <- c(credit_p2_2probit_llk_ratio_max,maxresult)
  credit_p2_2probit_llk_ratio_min <- c(credit_p2_2probit_llk_ratio_min,minresult)
  credit_p2_2probit_llk_ratio_average <- c(credit_p2_2probit_llk_ratio_average,meanresult)
}

credit_p2_uniform_llk_ratio_median <- c()
credit_p2_uniform_llk_ratio_max <- c()
credit_p2_uniform_llk_ratio_min <- c()
credit_p2_uniform_llk_ratio_average <- c()


median(abs(llk(X, credit_p2$post.mean,1,y)-
  llk(X,get(paste("credit_p_2_uniform_frac",frac[i],"chain_",j,sep = ""))$post.mean,1,y)))

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    # b <- sum(llk(a,get(paste("credit_p_2_uniform_frac",frac[i],"chain_",j,sep = ""))$post.mean,2))/sum(llk(a, credit_p2$post.mean,2))
    b<- median(abs(Bayes_factor(X, credit_p2$post.mean,2,y)-
                     Bayes_factor(X,get(paste("credit_p_2_uniform_frac",frac[i],"chain_",j,sep = ""))$post.mean,2,y)))
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p2_uniform_llk_ratio_median <- c(credit_p2_uniform_llk_ratio_median,median_result)
  credit_p2_uniform_llk_ratio_max <- c(credit_p2_uniform_llk_ratio_max,maxresult)
  credit_p2_uniform_llk_ratio_min <- c(credit_p2_uniform_llk_ratio_min,minresult)
  credit_p2_uniform_llk_ratio_average <- c(credit_p2_uniform_llk_ratio_average,meanresult)
}

credit_p2_one_shot_llk_ratio_median <- c()
credit_p2_one_shot_llk_ratio_max <- c()
credit_p2_one_shot_llk_ratio_min <- c()
credit_p2_one_shot_llk_ratio_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    # b <- sum(llk(a,get(paste("credit_p_2_one_shot_frac",frac[i],"chain_",j,sep = ""))$post.mean,2))/sum(llk(a, credit_p2$post.mean,2))
    b<- median(abs(Bayes_factor(X, credit_p2$post.mean,2,y)-
                     Bayes_factor(X,get(paste("credit_p_2_one_shot_frac",frac[i],"chain_",j,sep = ""))$post.mean,2,y)))
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p2_one_shot_llk_ratio_median <- c(credit_p2_one_shot_llk_ratio_median,median_result)
  credit_p2_one_shot_llk_ratio_max <- c(credit_p2_one_shot_llk_ratio_max,maxresult)
  credit_p2_one_shot_llk_ratio_min <- c(credit_p2_one_shot_llk_ratio_min,minresult)
  credit_p2_one_shot_llk_ratio_average <- c(credit_p2_one_shot_llk_ratio_average,meanresult)
}

p2_llk_ratio_credit_data <- data.frame(mean=c(log(credit_p2_2probit_llk_ratio_average),
                                                log(credit_p2_one_shot_llk_ratio_average),
                                                 log(credit_p2_uniform_llk_ratio_average)),
                                          median=c(log(credit_p2_2probit_llk_ratio_median),
                                                   log(credit_p2_one_shot_llk_ratio_median),
                                                   log(credit_p2_uniform_llk_ratio_median)),
                                          max=c(log(credit_p2_2probit_llk_ratio_max),
                                                log(credit_p2_one_shot_llk_ratio_max),
                                                log(credit_p2_uniform_llk_ratio_max)),
                                          min=c(log(credit_p2_2probit_llk_ratio_min),
                                                log(credit_p2_one_shot_llk_ratio_min),
                                                log(credit_p2_uniform_llk_ratio_min)),
                                          label=c(rep("2-probit",length(frac)),rep("one-shot",length(frac)),rep("uniform",length(frac))),
                                          frac=c(rep(frac,3)))
ggplot(p2_llk_ratio_credit_data, aes(frac, mean, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="likelihood ratio", colour="method", fill="method") +
  theme_bw()+coord_cartesian(ylim=c(0, 100))+
  ggtitle("(log)-Likelihood ratio of creditcard data (p=2)")
length(frac)
credit_p2_2probit_llk_ratio_average


p2_llk_ratio_credit_data <- data.frame(mean=c(credit_p2_2probit_llk_ratio_average,
                                              credit_p2_one_shot_llk_ratio_average,
                                              credit_p2_uniform_llk_ratio_average),
                                       median=c(credit_p2_2probit_llk_ratio_median,
                                                credit_p2_one_shot_llk_ratio_median,
                                                credit_p2_uniform_llk_ratio_median),
                                       max=c(credit_p2_2probit_llk_ratio_max,
                                             credit_p2_one_shot_llk_ratio_max,
                                             credit_p2_uniform_llk_ratio_max),
                                       min=c(credit_p2_2probit_llk_ratio_min,
                                             credit_p2_one_shot_llk_ratio_min,
                                             credit_p2_uniform_llk_ratio_min),
                                       label=c(rep("2-probit",length(frac)),rep("one-shot",length(frac)),rep("uniform",length(frac))),
                                       frac=c(rep(frac,3)))

ggplot(p2_llk_ratio_credit_data, aes(frac, median, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +coord_cartesian(ylim=c(0, 50))+
  labs(x="fraction", y="likelihood ratio", colour="method", fill="method") +
  theme_bw()+
  ggtitle("(log)-Likelihood ratio of creditcard data (p=2)")



ggplot(p2_llk_ratio_credit_data, aes(frac, mean, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="(log)-Bayes Factor", colour="method", fill="method") +
  theme_bw()+coord_cartesian(ylim=c(0, 50))+
  ggtitle("(log)-Bayes Factor of creditcard data (p=2)")


######poserior covariance difference

credit_p2_pos_cov_diff_norm_median <- c()
credit_p2_pos_cov_diff_norm_max <- c()
credit_p2_pos_cov_diff_norm_min <- c()
credit_p2_pos_cov_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("credit_p_2_2probit_frac",frac[i],"chain_",j,sep = ""))$chain)-
                cov(credit_p2$chain),type = "2")
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p2_pos_cov_diff_norm_median <- c(credit_p2_pos_cov_diff_norm_median,median_result)
  credit_p2_pos_cov_diff_norm_max <- c(credit_p2_pos_cov_diff_norm_max,maxresult)
  credit_p2_pos_cov_diff_norm_min <- c(credit_p2_pos_cov_diff_norm_min,minresult)
  credit_p2_pos_cov_diff_norm_average <- c(credit_p2_pos_cov_diff_norm_average,averageresult)
  
  
}

credit_p_2_one_shot_cov_diff_norm_median <- c()
credit_p_2_one_shot_cov_diff_norm_max <- c()
credit_p_2_one_shot_cov_diff_norm_min <- c()
credit_p_2_one_shot_cov_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("credit_p_2_one_shot_frac",frac[i],"chain_",j,sep = ""))$chain)-
                cov(credit_p2$chain),type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p_2_one_shot_cov_diff_norm_median <- c(credit_p_2_one_shot_cov_diff_norm_median,median_result)
  credit_p_2_one_shot_cov_diff_norm_max <- c(credit_p_2_one_shot_cov_diff_norm_max,maxresult)
  credit_p_2_one_shot_cov_diff_norm_min <- c(credit_p_2_one_shot_cov_diff_norm_min,minresult)
  credit_p_2_one_shot_cov_diff_norm_average <- c(credit_p_2_one_shot_cov_diff_norm_average,meanresult)
}

credit_p2_lpcoreset_pos_cov_diff_norm_median <- c()
credit_p2_lpcoreset_pos_cov_diff_norm_max <- c()
credit_p2_lpcoreset_pos_cov_diff_norm_min <- c()
credit_p2_lpcoreset_pos_cov_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("credit_p_2_rootl2_frac",frac[i],"chain_",j,sep = ""))$chain)-
                cov(credit_p2$chain),type = "2")
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p2_lpcoreset_pos_cov_diff_norm_median <- c(credit_p2_lpcoreset_pos_cov_diff_norm_median,median_result)
  credit_p2_lpcoreset_pos_cov_diff_norm_max <- c(credit_p2_lpcoreset_pos_cov_diff_norm_max,maxresult)
  credit_p2_lpcoreset_pos_cov_diff_norm_min <- c(credit_p2_lpcoreset_pos_cov_diff_norm_min,minresult)
  credit_p2_lpcoreset_pos_cov_diff_norm_average <- c(credit_p2_lpcoreset_pos_cov_diff_norm_average,averageresult)
  
}

credit_p2_uniform_pos_cov_diff_norm_median <- c()
credit_p2_uniform_pos_cov_diff_norm_max <- c()
credit_p2_uniform_pos_cov_diff_norm_min <- c()
credit_p2_uniform_pos_cov_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("credit_p_2_uniform_frac",frac[i],"chain_",j,sep = ""))$chain)-
                cov(credit_p2$chain),type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p2_uniform_pos_cov_diff_norm_median <- c(credit_p2_uniform_pos_cov_diff_norm_median,median_result)
  credit_p2_uniform_pos_cov_diff_norm_max <- c(credit_p2_uniform_pos_cov_diff_norm_max,maxresult)
  credit_p2_uniform_pos_cov_diff_norm_min <- c(credit_p2_uniform_pos_cov_diff_norm_min,minresult)
  credit_p2_uniform_pos_cov_diff_norm_average <- c(credit_p2_uniform_pos_cov_diff_norm_average,averageresult)
}

p2_covdiff_credit_data <- data.frame(mean=c(credit_p2_pos_cov_diff_norm_median,
                                            credit_p_2_one_shot_cov_diff_norm_average,
                                               credit_p2_uniform_pos_cov_diff_norm_average),
                                        median=c(credit_p2_pos_cov_diff_norm_median,
                                                 credit_p_2_one_shot_cov_diff_norm_median,
                                                 credit_p2_uniform_pos_cov_diff_norm_median),
                                        max=c(credit_p2_pos_cov_diff_norm_max,
                                              credit_p_2_one_shot_cov_diff_norm_max,
                                              credit_p2_uniform_pos_cov_diff_norm_max),
                                        min=c(credit_p2_pos_cov_diff_norm_min,
                                              credit_p_2_one_shot_cov_diff_norm_min,
                                              credit_p2_uniform_pos_cov_diff_norm_min),
                                        label=c(rep("2-probit",length(frac)),rep("one-shot",length(frac)),rep("uniform",length(frac))),
                                        frac=c(rep(frac,3)))
p2_covdiff_credit_data <- p2_covdiff_credit_data[-1,]

ggplot(p2_covdiff_credit_data, aes(frac, median, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior covariance", colour="method", fill="method") +
  theme_bw()+coord_cartesian(ylim=c(0, 100))+
  ggtitle("Credit card data (p=2)")

################plot of p=1


credit_p1_1probit_llk_ratio_median <- c()
credit_p1_1probit_llk_ratio_max <- c()
credit_p1_1probit_llk_ratio_min <- c()
credit_p1_1probit_llk_ratio_average <- c()
sum(llk(X, credit_p1$post.mean,1))
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    b<- median(abs(llk(X, credit_p1$post.mean,1,y)-
                     llk(X,get(paste("credit_p_1_1probit_frac",frac[i],"chain_",j,sep = ""))$post.mean,1,y)))
    
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p1_1probit_llk_ratio_median <- c(credit_p1_1probit_llk_ratio_median,median_result)
  credit_p1_1probit_llk_ratio_max <- c(credit_p1_1probit_llk_ratio_max,maxresult)
  credit_p1_1probit_llk_ratio_min <- c(credit_p1_1probit_llk_ratio_min,minresult)
  credit_p1_1probit_llk_ratio_average <- c(credit_p1_1probit_llk_ratio_average,meanresult)
}

credit_p1_uniform_llk_ratio_median <- c()
credit_p1_uniform_llk_ratio_max <- c()
credit_p1_uniform_llk_ratio_min <- c()
credit_p1_uniform_llk_ratio_average <- c()


for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    # b <- sum(llk(a,get(paste("credit_p_1_uniform_frac",frac[i],"chain_",j,sep = ""))$post.mean,1))/sum(llk(a, credit_p1$post.mean,1))
    b<- median(abs(llk(X, credit_p1$post.mean,1,y)-
                     llk(X,get(paste("credit_p_1_uniform_frac",frac[i],"chain_",j,sep = ""))$post.mean,1,y)))
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p1_uniform_llk_ratio_median <- c(credit_p1_uniform_llk_ratio_median,median_result)
  credit_p1_uniform_llk_ratio_max <- c(credit_p1_uniform_llk_ratio_max,maxresult)
  credit_p1_uniform_llk_ratio_min <- c(credit_p1_uniform_llk_ratio_min,minresult)
  credit_p1_uniform_llk_ratio_average <- c(credit_p1_uniform_llk_ratio_average,meanresult)
}

credit_p1_one_shot_llk_ratio_median <- c()
credit_p1_one_shot_llk_ratio_max <- c()
credit_p1_one_shot_llk_ratio_min <- c()
credit_p1_one_shot_llk_ratio_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    # b <- sum(llk(a,get(paste("credit_p_1_one_shot_frac",frac[i],"chain_",j,sep = ""))$post.mean,1))/sum(llk(a, credit_p1$post.mean,1))
    b<- median(abs(llk(X, credit_p1$post.mean,1,y)-
                     llk(X,get(paste("credit_p_1_one_shot_frac",frac[i],"chain_",j,sep = ""))$post.mean,1,y)))
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p1_one_shot_llk_ratio_median <- c(credit_p1_one_shot_llk_ratio_median,median_result)
  credit_p1_one_shot_llk_ratio_max <- c(credit_p1_one_shot_llk_ratio_max,maxresult)
  credit_p1_one_shot_llk_ratio_min <- c(credit_p1_one_shot_llk_ratio_min,minresult)
  credit_p1_one_shot_llk_ratio_average <- c(credit_p1_one_shot_llk_ratio_average,meanresult)
}

credit_p1_rootl2_llk_ratio_median <- c()
credit_p1_rootl2_llk_ratio_max <- c()
credit_p1_rootl2_llk_ratio_min <- c()
credit_p1_rootl2_llk_ratio_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    # b <- sum(llk(a,get(paste("credit_p_1_one_shot_frac",frac[i],"chain_",j,sep = ""))$post.mean,1))/sum(llk(a, credit_p1$post.mean,1))
    b<- median(abs(llk(X, credit_p1$post.mean,1,y)-
                     llk(X,get(paste("credit_p_1_rootl2_frac",frac[i],"chain_",j,sep = ""))$post.mean,1,y)))
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p1_rootl2_llk_ratio_median <- c(credit_p1_rootl2_llk_ratio_median,median_result)
  credit_p1_rootl2_llk_ratio_max <- c(credit_p1_rootl2_llk_ratio_max,maxresult)
  credit_p1_rootl2_llk_ratio_min <- c(credit_p1_rootl2_llk_ratio_min,minresult)
  credit_p1_rootl2_llk_ratio_average <- c(credit_p1_rootl2_llk_ratio_average,meanresult)
}


p1_llk_ratio_credit_data <- data.frame(mean=c(credit_p1_1probit_llk_ratio_average,
                                              credit_p1_one_shot_llk_ratio_average,
                                              credit_p1_uniform_llk_ratio_average,
                                              credit_p1_rootl2_llk_ratio_average),
                                       median=c(credit_p1_1probit_llk_ratio_median,
                                                credit_p1_one_shot_llk_ratio_median,
                                                credit_p1_uniform_llk_ratio_median,
                                                credit_p1_rootl2_llk_ratio_median),
                                       max=c(credit_p1_1probit_llk_ratio_max,
                                             credit_p1_one_shot_llk_ratio_max,
                                             credit_p1_uniform_llk_ratio_max,
                                             credit_p1_rootl2_llk_ratio_max),
                                       min=c(credit_p1_1probit_llk_ratio_min,
                                             credit_p1_one_shot_llk_ratio_min,
                                             credit_p1_uniform_llk_ratio_min,
                                             credit_p1_rootl2_llk_ratio_min),
                                       label=c(rep("1-probit",length(frac)),
                                               rep("one-shot",length(frac)),
                                               rep("uniform",length(frac)),
                                               rep("rootl2",length(frac))),
                                       frac=c(rep(frac,4)))

p1_llk_ratio_credit_data <- data.frame(mean=c(log(credit_p1_1probit_llk_ratio_average),
                                              log(credit_p1_one_shot_llk_ratio_average),
                                              log(credit_p1_uniform_llk_ratio_average)),
                                       median=c(log(credit_p1_1probit_llk_ratio_median),
                                                log(credit_p1_one_shot_llk_ratio_median),
                                                log(credit_p1_uniform_llk_ratio_median)),
                                       max=c(log(credit_p1_1probit_llk_ratio_max),
                                             log(credit_p1_one_shot_llk_ratio_max),
                                             log(credit_p1_uniform_llk_ratio_max)),
                                       min=c(log(credit_p1_1probit_llk_ratio_min),
                                             log(credit_p1_one_shot_llk_ratio_min),
                                             log(credit_p1_uniform_llk_ratio_min)),
                                       label=c(rep("1-probit",length(frac)),rep("one-shot",length(frac)),rep("uniform",length(frac))),
                                       frac=c(rep(frac,3)))


ggplot(p1_llk_ratio_credit_data, aes(frac, mean, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="(log)-Bayes Factor", colour="method", fill="method") +
  theme_bw()+
  ggtitle("(log)-Bayes Factor of creditcard data (p=1)")



###########post mean

credit_p_1_one_shot_mean_diff_norm_median <- c()
credit_p_1_one_shot_mean_diff_norm_max <- c()
credit_p_1_one_shot_mean_diff_norm_min <- c()
credit_p_1_one_shot_mean_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("credit_p_1_one_shot_frac",frac[i],"chain_",j,sep = ""))$post.mean-credit_p1$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p_1_one_shot_mean_diff_norm_median <- c(credit_p_1_one_shot_mean_diff_norm_median,median_result)
  credit_p_1_one_shot_mean_diff_norm_max <- c(credit_p_1_one_shot_mean_diff_norm_max,maxresult)
  credit_p_1_one_shot_mean_diff_norm_min <- c(credit_p_1_one_shot_mean_diff_norm_min,minresult)
  credit_p_1_one_shot_mean_diff_norm_average <- c(credit_p_1_one_shot_mean_diff_norm_average,meanresult)
}

credit_p_1_rootl2_mean_diff_norm_median <- c()
credit_p_1_rootl2_mean_diff_norm_max <- c()
credit_p_1_rootl2_mean_diff_norm_min <- c()
credit_p_1_rootl2_mean_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    b <- norm(get(paste("credit_p_1_rootl2_frac",frac[i],"chain_",j,sep = ""))$post.mean-credit_p1$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p_1_rootl2_mean_diff_norm_median <- c(credit_p_1_rootl2_mean_diff_norm_median,median_result)
  credit_p_1_rootl2_mean_diff_norm_max <- c(credit_p_1_rootl2_mean_diff_norm_max,maxresult)
  credit_p_1_rootl2_mean_diff_norm_min <- c(credit_p_1_rootl2_mean_diff_norm_min,minresult)
  credit_p_1_rootl2_mean_diff_norm_average <- c(credit_p_1_rootl2_mean_diff_norm_average,meanresult)
}


credit_p_1_2probit_mean_diff_norm_median <- c()
credit_p_1_2probit_mean_diff_norm_max <- c()
credit_p_1_2probit_mean_diff_norm_min <- c()
credit_p_1_2probit_mean_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    b <- norm(get(paste("credit_p_1_2probit_frac",frac[i],"chain_",j,sep = ""))$post.mean-credit_p1$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p_1_2probit_mean_diff_norm_median <- c(credit_p_1_2probit_mean_diff_norm_median,median_result)
  credit_p_1_2probit_mean_diff_norm_max <- c(credit_p_1_2probit_mean_diff_norm_max,maxresult)
  credit_p_1_2probit_mean_diff_norm_min <- c(credit_p_1_2probit_mean_diff_norm_min,minresult)
  credit_p_1_2probit_mean_diff_norm_average <- c(credit_p_1_2probit_mean_diff_norm_average,meanresult)
}


credit_p_1_uniform_mean_diff_norm_median <- c()
credit_p_1_uniform_mean_diff_norm_max <- c()
credit_p_1_uniform_mean_diff_norm_min <- c()
credit_p_1_uniform_mean_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    b <- norm(get(paste("credit_p_1_uniform_frac",frac[i],"chain_",j,sep = ""))$post.mean-credit_p1$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p_1_uniform_mean_diff_norm_median <- c(credit_p_1_uniform_mean_diff_norm_median,median_result)
  credit_p_1_uniform_mean_diff_norm_max <- c(credit_p_1_uniform_mean_diff_norm_max,maxresult)
  credit_p_1_uniform_mean_diff_norm_min <- c(credit_p_1_uniform_mean_diff_norm_min,minresult)
  credit_p_1_uniform_mean_diff_norm_average <- c(credit_p_1_uniform_mean_diff_norm_average,meanresult)
}


p1_meandiff_credit_data <- data.frame(mean=c(credit_p_1_2probit_mean_diff_norm_average,
                                             credit_p_1_one_shot_mean_diff_norm_average,
                                             credit_p_1_uniform_mean_diff_norm_average,
                                             credit_p_1_rootl2_mean_diff_norm_average),
                                      median=c(credit_p_1_2probit_mean_diff_norm_median,
                                               credit_p_1_one_shot_mean_diff_norm_median,
                                               credit_p_1_uniform_mean_diff_norm_median,
                                               credit_p_1_rootl2_mean_diff_norm_median),
                                      max=c(credit_p_1_2probit_mean_diff_norm_max,
                                            credit_p_1_one_shot_mean_diff_norm_max,
                                            credit_p_1_uniform_mean_diff_norm_max,
                                            credit_p_1_rootl2_mean_diff_norm_max),
                                      min=c(credit_p_1_2probit_mean_diff_norm_min,
                                            credit_p_1_one_shot_mean_diff_norm_min,
                                            credit_p_1_uniform_mean_diff_norm_min,
                                            credit_p_1_rootl2_mean_diff_norm_min),
                                      label=c(rep("2-probit",length(frac)),rep("one-shot",length(frac)),rep("uniform",length(frac)),
                                              rep("rootl2",length(frac))),
                                      frac=c(rep(frac,4)))
p2_meandiff_credit_data


ggplot(p1_meandiff_credit_data, aes(frac, median, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior mean", colour="method", fill="method") +
  theme_bw()+
  ggtitle("Credit card data (p=1)")
######poserior covariance difference


credit_p1_pos_cov_diff_norm_median <- c()
credit_p1_pos_cov_diff_norm_max <- c()
credit_p1_pos_cov_diff_norm_min <- c()
credit_p1_pos_cov_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("credit_p_1_2probit_frac",frac[i],"chain_",j,sep = ""))$chain)-
                cov(credit_p1$chain),type = "2")
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p1_pos_cov_diff_norm_median <- c(credit_p1_pos_cov_diff_norm_median,median_result)
  credit_p1_pos_cov_diff_norm_max <- c(credit_p1_pos_cov_diff_norm_max,maxresult)
  credit_p1_pos_cov_diff_norm_min <- c(credit_p1_pos_cov_diff_norm_min,minresult)
  credit_p1_pos_cov_diff_norm_average <- c(credit_p1_pos_cov_diff_norm_average,averageresult)
  
  
}

credit_p_1_one_shot_cov_diff_norm_median <- c()
credit_p_1_one_shot_cov_diff_norm_max <- c()
credit_p_1_one_shot_cov_diff_norm_min <- c()
credit_p_1_one_shot_cov_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    b <- norm(cov(get(paste("credit_p_1_one_shot_frac",frac[i],"chain_",j,sep = ""))$chain)-
                cov(credit_p1$chain),type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p_1_one_shot_cov_diff_norm_median <- c(credit_p_1_one_shot_cov_diff_norm_median,median_result)
  credit_p_1_one_shot_cov_diff_norm_max <- c(credit_p_1_one_shot_cov_diff_norm_max,maxresult)
  credit_p_1_one_shot_cov_diff_norm_min <- c(credit_p_1_one_shot_cov_diff_norm_min,minresult)
  credit_p_1_one_shot_cov_diff_norm_average <- c(credit_p_1_one_shot_cov_diff_norm_average,meanresult)
}

credit_p1_lpcoreset_pos_cov_diff_norm_median <- c()
credit_p1_lpcoreset_pos_cov_diff_norm_max <- c()
credit_p1_lpcoreset_pos_cov_diff_norm_min <- c()
credit_p1_lpcoreset_pos_cov_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("credit_p_1_rootl2_frac",frac[i],"chain_",j,sep = ""))$chain)-
                cov(credit_p1$chain),type = "2")
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p1_lpcoreset_pos_cov_diff_norm_median <- c(credit_p1_lpcoreset_pos_cov_diff_norm_median,median_result)
  credit_p1_lpcoreset_pos_cov_diff_norm_max <- c(credit_p1_lpcoreset_pos_cov_diff_norm_max,maxresult)
  credit_p1_lpcoreset_pos_cov_diff_norm_min <- c(credit_p1_lpcoreset_pos_cov_diff_norm_min,minresult)
  credit_p1_lpcoreset_pos_cov_diff_norm_average <- c(credit_p1_lpcoreset_pos_cov_diff_norm_average,averageresult)
  
}

credit_p1_uniform_pos_cov_diff_norm_median <- c()
credit_p1_uniform_pos_cov_diff_norm_max <- c()
credit_p1_uniform_pos_cov_diff_norm_min <- c()
credit_p1_uniform_pos_cov_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("credit_p_1_uniform_frac",frac[i],"chain_",j,sep = ""))$chain)-
                cov(credit_p1$chain),type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p1_uniform_pos_cov_diff_norm_median <- c(credit_p1_uniform_pos_cov_diff_norm_median,median_result)
  credit_p1_uniform_pos_cov_diff_norm_max <- c(credit_p1_uniform_pos_cov_diff_norm_max,maxresult)
  credit_p1_uniform_pos_cov_diff_norm_min <- c(credit_p1_uniform_pos_cov_diff_norm_min,minresult)
  credit_p1_uniform_pos_cov_diff_norm_average <- c(credit_p1_uniform_pos_cov_diff_norm_average,averageresult)
}

p1_covdiff_credit_data <- data.frame(mean=c(credit_p1_pos_cov_diff_norm_median,
                                            credit_p_1_one_shot_cov_diff_norm_average,
                                            credit_p1_uniform_pos_cov_diff_norm_average,
                                            credit_p1_lpcoreset_pos_cov_diff_norm_average),
                                     median=c(credit_p1_pos_cov_diff_norm_median,
                                              credit_p_1_one_shot_cov_diff_norm_median,
                                              credit_p1_uniform_pos_cov_diff_norm_median,
                                              credit_p1_lpcoreset_pos_cov_diff_norm_median),
                                     max=c(credit_p1_pos_cov_diff_norm_max,
                                           credit_p_1_one_shot_cov_diff_norm_max,
                                           credit_p1_uniform_pos_cov_diff_norm_max,
                                           credit_p1_lpcoreset_pos_cov_diff_norm_max),
                                     min=c(credit_p1_pos_cov_diff_norm_min,
                                           credit_p_1_one_shot_cov_diff_norm_min,
                                           credit_p1_uniform_pos_cov_diff_norm_min,
                                           credit_p1_lpcoreset_pos_cov_diff_norm_min),
                                     label=c(rep("2-probit",length(frac)),rep("one-shot",length(frac)),rep("uniform",length(frac)),
                                             rep("rootl2",length(frac))),
                                     frac=c(rep(frac,4)))


ggplot(p1_covdiff_credit_data, aes(frac, mean, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior covariance", colour="method", fill="method") +
  theme_bw()+
  ggtitle("Credit card data(p=1)")

########fixed p=2#########
#####algorithm 1#####
for (i in 1:length(frac)) {
  
  print(paste("Model_Frac:",i,"chain:",j,sep = ""))
  sample_index<- Compute_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=2)
  Model <- multi_chain(N_sim=N_sim,burn_in=burn_in,
                       X=X[sample_index,], y=y[sample_index],initial_theda=probit_coefficients,
                       true_theta=probit_coefficients,Lp=3,M_iter=100,range=c(0.1,5),step=0.02,times = 1)
  assign(paste("credit_p2_nofix_frac",frac[i],sep = ""), Model)
  
}

credit_p_2probit_summary <- c()
for (i in 1:length(frac)){
  
  b <- get(paste("credit_p2_nofix_frac",frac[i],sep = ""))$Lp.mean
  
  credit_p_2probit_summary<- rbind(credit_p_2probit_summary,b)
  
}
plot(frac,credit_p_2probit_summary,type = "l")

for (i in 1:length(frac)) {
  
  print(paste("Model_Frac:",i,"chain:",j,sep = ""))
  sample_index<- Compute_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=1)
  Model <- multi_chain(N_sim=N_sim,burn_in=burn_in,
                       X=X[sample_index,], y=y[sample_index],initial_theda=probit_coefficients,
                       true_theta=probit_coefficients,Lp=3,M_iter=100,range=c(0.1,5),step=0.02,times = 1)
  assign(paste("credit_p1_nofix_frac",frac[i],sep = ""), Model)
  
}

credit_p_1probit_summary <- c()
for (i in 1:length(frac)){
  
  b <- get(paste("credit_p1_nofix_frac",frac[i],sep = ""))$Lp.mean
  
  credit_p_1probit_summary<- rbind(credit_p_1probit_summary,b)
  
}
plot(frac,credit_p_1probit_summary,type = "l")

#####rootl2#########
for (i in 1:length(frac)) {
  print(paste("Model_Frac:",i,"chain:",j,sep = ""))
  sample_index<- Lp_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=2)
  Model <- multi_chain(N_sim=N_sim,burn_in=burn_in,
                       X=X[sample_index,], y=y[sample_index],initial_theda=probit_coefficients,
                       true_theta=probit_coefficients,Lp=3,M_iter=100,range=c(0.1,5),step=0.02,times = 1)
  assign(paste("credit_p_rootl2_nofix_frac",frac[i],sep = ""), Model)
  
}
credit_p_rootl2_summary <- c()
for (i in 1:length(frac)){
  
  b <- get(paste("credit_p_rootl2_nofix_frac",frac[i],sep = ""))$Lp.mean
  
  credit_p_rootl2_summary<- rbind(credit_p_rootl2_summary,b)
  
}
plot(frac,credit_p_rootl2_summary,type = "l")
#########uniform sampling
for (i in 1:length(frac)) {
  print(paste("Model_Frac:",i,"chain:",j,sep = ""))
  sample_index<- Uniform_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=2)
  Model <- multi_chain(N_sim=N_sim,burn_in=burn_in,
                       X=X[sample_index,], y=y[sample_index],initial_theda=probit_coefficients,
                       true_theta=probit_coefficients,Lp=3,M_iter=100,range=c(0.1,5),step=0.02,times = 1)
  assign(paste("credit_p_uniform_nofix_frac",frac[i],sep = ""), Model)
}
credit_p_uniform_summary <- c()
for (i in 1:length(frac)){
  
  b <- get(paste("credit_p_uniform_nofix_frac",frac[i],sep = ""))$Lp.mean
  
  credit_p_uniform_summary<- rbind(credit_p_uniform_summary,b)
  
}

plot(frac,credit_p_2probit_summary,type = "l",ylim = c(0,3),xlab = "sample fraction", ylab = "estimated P")
lines(frac,credit_p_rootl2_summary,col=2)
lines(frac,credit_p_uniform_summary,col=3)
lines(frac,credit_p_1probit_summary,col=4)

legend("topright",legend=c("2-probit","root l2","uniform"), col=c(1,2,3), 
       bty = 'n', lwd = 2, inset = c(0, 0), lty = 1, cex = 0.73)


length(credit_p_uniform_nofix_frac0.00017211348$Lp.mean.chain)

plot(1:N_sim,get(paste("credit_p_uniform_nofix_frac",frac[58],sep = ""))$Lp.mean.chain,type = "l")
lines(1:N_sim,get(paste("credit_p_rootl2_nofix_frac",frac[58],sep = ""))$Lp.mean.chain,col=2)
lines(1:N_sim,get(paste("credit_p2_nofix_frac",frac[58],sep = ""))$Lp.mean.chain,col=3)

plot(1:N_sim,get(paste("credit_p_uniform_nofix_frac",frac[1],sep = ""))$Lp.mean.chain,type = "l")
lines(1:N_sim,get(paste("credit_p_rootl2_nofix_frac",frac[1],sep = ""))$Lp.mean.chain,col=2)
lines(1:N_sim,get(paste("credit_p2_nofix_frac",frac[1],sep = ""))$Lp.mean.chain,col=3)



credit_uniform_pos_mean_diff_norm <- c()
for (i in 1:length(frac)){
  
  b <- norm(get(paste("credit_p_uniform_nofix_frac",frac[i],sep = ""))$Beta.mean-credit_p2$post.mean,type = "2")
  
  credit_uniform_pos_mean_diff_norm<- rbind(credit_uniform_pos_mean_diff_norm,b)
  
}
credit_uniform_pos_mean_diff_norm
plot(frac,credit_uniform_pos_mean_diff_norm,type = "l")





##############p=1######################
for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Lp_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=1)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y[sample_index],
                             initial_theda=probit_coefficients,
                             true_theta=probit_coefficients,Lp=1)
    assign(paste("credit_p_1_rootl2_frac",frac[i],"chain_",j,sep = ""), Model)
  }
}

for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Uniform_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=1)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y[sample_index],
                             initial_theda=probit_coefficients,
                             true_theta=probit_coefficients,Lp=1)
    assign(paste("credit_p_1_uniform_frac",frac[i],"chain_",j,sep = ""), Model)
  }
}


for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Compute_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=1)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y[sample_index],
                             initial_theda=probit_coefficients,
                             true_theta=probit_coefficients,Lp=1)
    assign(paste("credit_p_1_1probit_frac",frac[i],"chain_",j,sep = ""), Model)
  }
}

credit_p_1_rootl2_mean_diff_norm_median <- c()
credit_p_1_rootl2_mean_diff_norm_max <- c()
credit_p_1_rootl2_mean_diff_norm_min <- c()
credit_p_1_rootl2_mean_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("credit_p_1_rootl2_frac",frac[i],"chain_",j,sep = ""))$post.mean-credit_p1$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p_1_rootl2_mean_diff_norm_median <- c(credit_p_1_rootl2_mean_diff_norm_median,median_result)
  credit_p_1_rootl2_mean_diff_norm_max <- c(credit_p_1_rootl2_mean_diff_norm_max,maxresult)
  credit_p_1_rootl2_mean_diff_norm_min <- c(credit_p_1_rootl2_mean_diff_norm_min,minresult)
  credit_p_1_rootl2_mean_diff_norm_average <- c(credit_p_1_rootl2_mean_diff_norm_average,meanresult)
}


credit_p_1_1probit_mean_diff_norm_median <- c()
credit_p_1_1probit_mean_diff_norm_max <- c()
credit_p_1_1probit_mean_diff_norm_min <- c()
credit_p_1_1probit_mean_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("credit_p_1_1probit_frac",frac[i],"chain_",j,sep = ""))$post.mean-credit_p1$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p_1_1probit_mean_diff_norm_median <- c(credit_p_1_1probit_mean_diff_norm_median,median_result)
  credit_p_1_1probit_mean_diff_norm_max <- c(credit_p_1_1probit_mean_diff_norm_max,maxresult)
  credit_p_1_1probit_mean_diff_norm_min <- c(credit_p_1_1probit_mean_diff_norm_min,minresult)
  credit_p_1_1probit_mean_diff_norm_average <- c(credit_p_1_1probit_mean_diff_norm_average,meanresult)
}


credit_p_1_uniform_mean_diff_norm_median <- c()
credit_p_1_uniform_mean_diff_norm_max <- c()
credit_p_1_uniform_mean_diff_norm_min <- c()
credit_p_1_uniform_mean_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("credit_p_1_uniform_frac",frac[i],"chain_",j,sep = ""))$post.mean-credit_p1$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p_1_uniform_mean_diff_norm_median <- c(credit_p_1_uniform_mean_diff_norm_median,median_result)
  credit_p_1_uniform_mean_diff_norm_max <- c(credit_p_1_uniform_mean_diff_norm_max,maxresult)
  credit_p_1_uniform_mean_diff_norm_min <- c(credit_p_1_uniform_mean_diff_norm_min,minresult)
  credit_p_1_uniform_mean_diff_norm_average <- c(credit_p_1_uniform_mean_diff_norm_average,meanresult)
}
plot(frac,credit_p_1_rootl2_mean_diff_norm_average,type = "l")
lines(frac,credit_p_1_uniform_mean_diff_norm_average,col=2)
p1_meandiff_credit_data <- data.frame(mean=c(credit_p_1_1probit_mean_diff_norm_average,
                                             credit_p_1_rootl2_mean_diff_norm_average,
                                             credit_p_1_uniform_mean_diff_norm_average),
                                      median=c(credit_p_1_1probit_mean_diff_norm_median,
                                               credit_p_1_rootl2_mean_diff_norm_median,
                                               credit_p_1_uniform_mean_diff_norm_median),
                                      max=c(credit_p_1_1probit_mean_diff_norm_max,
                                            credit_p_1_rootl2_mean_diff_norm_max,
                                            credit_p_1_uniform_mean_diff_norm_max),
                                      min=c(credit_p_1_1probit_mean_diff_norm_min,
                                            credit_p_1_rootl2_mean_diff_norm_min,
                                            credit_p_1_uniform_mean_diff_norm_min),
                                      label=c(rep("1-probit",length(frac)),rep("root l2",length(frac)),rep("uniform",length(frac))),
                                      frac=c(rep(frac,3)))
p1_meandiff_credit_data
ggplot(p1_meandiff_credit_data, aes(frac, median, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior mean", colour="method", fill="method") +
  theme_bw()

######poserior covariance difference
credit_p1_pos_cov_diff_norm_median <- c()
credit_p1_pos_cov_diff_norm_max <- c()
credit_p1_pos_cov_diff_norm_min <- c()
credit_p1_pos_cov_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    b <- norm(cov(get(paste("credit_p_1_1probit_frac",frac[i],"chain_",j,sep = ""))$chain)-
                cov(credit_p1$chain),type = "2")
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p1_pos_cov_diff_norm_median <- c(credit_p1_pos_cov_diff_norm_median,median_result)
  credit_p1_pos_cov_diff_norm_max <- c(credit_p1_pos_cov_diff_norm_max,maxresult)
  credit_p1_pos_cov_diff_norm_min <- c(credit_p1_pos_cov_diff_norm_min,minresult)
  credit_p1_pos_cov_diff_norm_average <- c(credit_p1_pos_cov_diff_norm_average,averageresult)
  
  
}

credit_p1_lpcoreset_pos_cov_diff_norm_median <- c()
credit_p1_lpcoreset_pos_cov_diff_norm_max <- c()
credit_p1_lpcoreset_pos_cov_diff_norm_min <- c()
credit_p1_lpcoreset_pos_cov_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    b <- norm(cov(get(paste("credit_p_1_rootl2_frac",frac[i],"chain_",j,sep = ""))$chain)-
                cov(credit_p1$chain),type = "2")
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p1_lpcoreset_pos_cov_diff_norm_median <- c(credit_p1_lpcoreset_pos_cov_diff_norm_median,median_result)
  credit_p1_lpcoreset_pos_cov_diff_norm_max <- c(credit_p1_lpcoreset_pos_cov_diff_norm_max,maxresult)
  credit_p1_lpcoreset_pos_cov_diff_norm_min <- c(credit_p1_lpcoreset_pos_cov_diff_norm_min,minresult)
  credit_p1_lpcoreset_pos_cov_diff_norm_average <- c(credit_p1_lpcoreset_pos_cov_diff_norm_average,averageresult)
  
}

credit_p1_uniform_pos_cov_diff_norm_median <- c()
credit_p1_uniform_pos_cov_diff_norm_max <- c()
credit_p1_uniform_pos_cov_diff_norm_min <- c()
credit_p1_uniform_pos_cov_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    b <- norm(cov(get(paste("credit_p_1_uniform_frac",frac[i],"chain_",j,sep = ""))$chain)-
                cov(credit_p1$chain),type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p1_uniform_pos_cov_diff_norm_median <- c(credit_p1_uniform_pos_cov_diff_norm_median,median_result)
  credit_p1_uniform_pos_cov_diff_norm_max <- c(credit_p1_uniform_pos_cov_diff_norm_max,maxresult)
  credit_p1_uniform_pos_cov_diff_norm_min <- c(credit_p1_uniform_pos_cov_diff_norm_min,minresult)
  credit_p1_uniform_pos_cov_diff_norm_average <- c(credit_p1_uniform_pos_cov_diff_norm_average,averageresult)
}

p1_covdiff_credit_data <- data.frame(mean=c(credit_p1_pos_cov_diff_norm_median,
                                            credit_p1_lpcoreset_pos_cov_diff_norm_average,
                                            credit_p1_uniform_pos_cov_diff_norm_average),
                                     median=c(credit_p1_pos_cov_diff_norm_median,
                                              credit_p1_lpcoreset_pos_cov_diff_norm_median,
                                              credit_p1_uniform_pos_cov_diff_norm_median),
                                     max=c(credit_p1_pos_cov_diff_norm_max,
                                           credit_p1_lpcoreset_pos_cov_diff_norm_max,
                                           credit_p1_uniform_pos_cov_diff_norm_max),
                                     min=c(credit_p1_pos_cov_diff_norm_min,
                                           credit_p1_lpcoreset_pos_cov_diff_norm_min,
                                           credit_p1_uniform_pos_cov_diff_norm_min),
                                     label=c(rep("1-probit",length(frac)),rep("root l2",length(frac)),rep("uniform",length(frac))),
                                     frac=c(rep(frac,3)))
p1_covdiff_credit_data
ggplot(p1_covdiff_credit_data, aes(frac, mean, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior covariance", colour="method", fill="method") +
  theme_bw()



######poserior covariance difference
credit_p1_pos_cov_diff_norm_median <- c()
credit_p1_pos_cov_diff_norm_max <- c()
credit_p1_pos_cov_diff_norm_min <- c()
credit_p1_pos_cov_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    b <- norm(cov(get(paste("credit_p_1_1probit_frac",frac[i],"chain_",j,sep = ""))$chain)-
                cov(credit_p1$chain),type = "2")
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p1_pos_cov_diff_norm_median <- c(credit_p1_pos_cov_diff_norm_median,median_result)
  credit_p1_pos_cov_diff_norm_max <- c(credit_p1_pos_cov_diff_norm_max,maxresult)
  credit_p1_pos_cov_diff_norm_min <- c(credit_p1_pos_cov_diff_norm_min,minresult)
  credit_p1_pos_cov_diff_norm_average <- c(credit_p1_pos_cov_diff_norm_average,averageresult)
  
  
}

credit_p1_lpcoreset_pos_cov_diff_norm_median <- c()
credit_p1_lpcoreset_pos_cov_diff_norm_max <- c()
credit_p1_lpcoreset_pos_cov_diff_norm_min <- c()
credit_p1_lpcoreset_pos_cov_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    b <- norm(cov(get(paste("credit_p_1_rootl2_frac",frac[i],"chain_",j,sep = ""))$chain)-
                cov(credit_p1$chain),type = "2")
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p1_lpcoreset_pos_cov_diff_norm_median <- c(credit_p1_lpcoreset_pos_cov_diff_norm_median,median_result)
  credit_p1_lpcoreset_pos_cov_diff_norm_max <- c(credit_p1_lpcoreset_pos_cov_diff_norm_max,maxresult)
  credit_p1_lpcoreset_pos_cov_diff_norm_min <- c(credit_p1_lpcoreset_pos_cov_diff_norm_min,minresult)
  credit_p1_lpcoreset_pos_cov_diff_norm_average <- c(credit_p1_lpcoreset_pos_cov_diff_norm_average,averageresult)
  
}

credit_p1_uniform_pos_cov_diff_norm_median <- c()
credit_p1_uniform_pos_cov_diff_norm_max <- c()
credit_p1_uniform_pos_cov_diff_norm_min <- c()
credit_p1_uniform_pos_cov_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    b <- norm(cov(get(paste("credit_p_1_uniform_frac",frac[i],"chain_",j,sep = ""))$chain)-
                cov(credit_p1$chain),type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN = 2, mean)
  credit_p1_uniform_pos_cov_diff_norm_median <- c(credit_p1_uniform_pos_cov_diff_norm_median,median_result)
  credit_p1_uniform_pos_cov_diff_norm_max <- c(credit_p1_uniform_pos_cov_diff_norm_max,maxresult)
  credit_p1_uniform_pos_cov_diff_norm_min <- c(credit_p1_uniform_pos_cov_diff_norm_min,minresult)
  credit_p1_uniform_pos_cov_diff_norm_average <- c(credit_p1_uniform_pos_cov_diff_norm_average,averageresult)
}

p1_covdiff_credit_data <- data.frame(mean=c(credit_p1_pos_cov_diff_norm_median,
                                            credit_p1_lpcoreset_pos_cov_diff_norm_average,
                                            credit_p1_uniform_pos_cov_diff_norm_average),
                                     median=c(credit_p1_pos_cov_diff_norm_median,
                                              credit_p1_lpcoreset_pos_cov_diff_norm_median,
                                              credit_p1_uniform_pos_cov_diff_norm_median),
                                     max=c(credit_p1_pos_cov_diff_norm_max,
                                           credit_p1_lpcoreset_pos_cov_diff_norm_max,
                                           credit_p1_uniform_pos_cov_diff_norm_max),
                                     min=c(credit_p1_pos_cov_diff_norm_min,
                                           credit_p1_lpcoreset_pos_cov_diff_norm_min,
                                           credit_p1_uniform_pos_cov_diff_norm_min),
                                     label=c(rep("1-probit",length(frac)),rep("root l2",length(frac)),rep("uniform",length(frac))),
                                     frac=c(rep(frac,3)))
p1_covdiff_credit_data
ggplot(p1_covdiff_credit_data, aes(frac, mean, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior covariance", colour="method", fill="method") +
  theme_bw()


############one-shot-comparison#########
times <- seq(1,15,2)

for (i in 1:length(frac)) {
  for (j in 1:5) {
    for (k in times) {
    coreset_size<- N*frac[i]
    print(paste("Model_Frac:",i,"chain:",j,"time:",k,sep = ""))
    scores<- one_shot_coreset_times(sketch_size = sketch_size, coreset_size = coreset_size, X=X, Lp_max = 2.5,
                                    times=k)
    sample_index <- sample(1:N,size = coreset_size,
                           replace = FALSE, prob = scores$score)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y[sample_index],
                             initial_theda=probit_coefficients,
                             true_theta=probit_coefficients,Lp=2)
    assign(paste("credit_p_2_one_shot",frac[i],"chain_",j,"time_",k,sep = ""), Model)
    assign(paste("credit_p_2_one_shot",frac[i],"chain_",j,"coreset_time",k,sep = ""), scores$modeling.time)
  }
  }
}


one_shot_time_average1 <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    b <- norm(get(paste("credit_p_2_one_shot",frac[i],"chain_",j,"time_",1,sep = ""))$post.mean-credit_p2$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  one_shot_time_average1 <- c(one_shot_time_average1,meanresult)
}

one_shot_time_average3 <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    b <- norm(get(paste("credit_p_2_one_shot",frac[i],"chain_",j,"time_",3,sep = ""))$post.mean-credit_p2$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  one_shot_time_average3 <- c(one_shot_time_average3,meanresult)
}
one_shot_time_average5 <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    b <- norm(get(paste("credit_p_2_one_shot",frac[i],"chain_",j,"time_",5,sep = ""))$post.mean-credit_p2$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  one_shot_time_average5 <- c(one_shot_time_average5,meanresult)
}
one_shot_time_average7 <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    b <- norm(get(paste("credit_p_2_one_shot",frac[i],"chain_",j,"time_",7,sep = ""))$post.mean-credit_p2$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  one_shot_time_average7 <- c(one_shot_time_average7,meanresult)
}
one_shot_time_average9 <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    b <- norm(get(paste("credit_p_2_one_shot",frac[i],"chain_",j,"time_",9,sep = ""))$post.mean-credit_p2$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  one_shot_time_average9 <- c(one_shot_time_average9,meanresult)
}
one_shot_time_average11 <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    b <- norm(get(paste("credit_p_2_one_shot",frac[i],"chain_",j,"time_",11,sep = ""))$post.mean-credit_p2$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  one_shot_time_average11 <- c(one_shot_time_average11,meanresult)
}
one_shot_time_average13 <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    b <- norm(get(paste("credit_p_2_one_shot",frac[i],"chain_",j,"time_",13,sep = ""))$post.mean-credit_p2$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  one_shot_time_average13 <- c(one_shot_time_average13,meanresult)
}
one_shot_time_average15 <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    b <- norm(get(paste("credit_p_2_one_shot",frac[i],"chain_",j,"time_",15,sep = ""))$post.mean-credit_p2$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  one_shot_time_average15 <- c(one_shot_time_average15,meanresult)
}

one_shot_plot_data <- data.frame(mean=c(one_shot_time_average1,
                                        one_shot_time_average3,
                                        one_shot_time_average5,
                                        one_shot_time_average7,
                                        one_shot_time_average9,
                                        one_shot_time_average11,
                                        one_shot_time_average13,
                                        one_shot_time_average15,
                                        credit_p_2_one_shot_mean_diff_norm_average),
                                     label=c(rep("grid-1",length(frac)),
                                             rep("grid-3",length(frac)),
                                             rep("grid-5",length(frac)),
                                             rep("grid-7",length(frac)),
                                             rep("grid-9",length(frac)),
                                             rep("grid-11",length(frac)),
                                             rep("grid-13",length(frac)),
                                             rep("grid-15",length(frac)),
                                             rep("one-shot",length(frac))),
                                     frac=c(rep(frac,9)))
one_shot_plot_data <- data.frame(mean=c(one_shot_time_average1,
                                        one_shot_time_average7,
                                        one_shot_time_average15,
                                        credit_p_2_one_shot_mean_diff_norm_average),
                                 label=c(rep("grid-1",length(frac)),
                                         rep("grid-7",length(frac)),
                                         rep("grid-15",length(frac)),
                                         rep("one-shot",length(frac))),
                                 frac=c(rep(frac,4)))




probit2_time_count <- c()

for (i in 1:length(frac)){
  medianchain <- c()

    b <- get(paste("credit_p_2_2probit_frac",frac[i],"chain_",1,sep = ""))$modeling.time[3]
    medianchain<- rbind(medianchain,b)
    probit2_time_count <- c(probit2_time_count,medianchain)
}
probit2_time_count

oneshot_time_count <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  
  b <- get(paste("credit_p_2_one_shot_frac",frac[i],"chain_",1,sep = ""))$modeling.time[3]
  medianchain<- rbind(medianchain,b)
  oneshot_time_count <- c(oneshot_time_count,medianchain)
}

rootl2_time_count <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  
  b <- get(paste("credit_p_2_rootl2_frac",frac[i],"chain_",1,sep = ""))$modeling.time[3]
  medianchain<- rbind(medianchain,b)
  rootl2_time_count <- c(rootl2_time_count,medianchain)
}

uniform_time_count <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  
  b <- get(paste("credit_p_2_uniform_frac",frac[i],"chain_",1,sep = ""))$modeling.time[3]
  medianchain<- rbind(medianchain,b)
  uniform_time_count <- c(uniform_time_count,medianchain)
}

modeling_time_plot <- data.frame(mean=c(probit2_time_count,
                                        oneshot_time_count,
                                        uniform_time_count,
                                        rootl2_time_count),
                                 label=c(rep("2-probit",length(frac)),rep("one-shot",length(frac)),
                                         rep("uniform",length(frac)),
                                         rep("rootl2",length(frac))),
                                 frac=c(rep(frac,4)))



ggplot(modeling_time_plot, aes(frac, mean, fill=factor(label), colour=factor(label))) +
  # geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  # geom_point(aes(x=0.05,y=500),colour="red")+
  labs(x="fraction", y="Seconds", colour="method", fill="method") +
  theme_bw()+ggtitle("Computation time of MCMC chain on different coresets")


one_shot_time_average1 <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  
  b <- get(paste("credit_p_2_one_shot",frac[i],"chain_",j,"coreset_time",1,sep = ""))[3]

  medianchain<- rbind(medianchain,b)
  one_shot_time_average1 <- c(one_shot_time_average1,medianchain)
}


one_shot_time_average7 <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  
  b <- get(paste("credit_p_2_one_shot",frac[i],"chain_",j,"coreset_time",7,sep = ""))[3]

  medianchain<- rbind(medianchain,b)
  one_shot_time_average7 <- c(one_shot_time_average7,medianchain)
}

one_shot_time_average15 <- c()

for (i in times){
  medianchain <- c()
  
  b <- get(paste("credit_p_2_one_shot",frac[1],"chain_",1,"coreset_time",i,sep = ""))[3]

  medianchain<- rbind(medianchain,b)
  one_shot_time_average15 <- c(one_shot_time_average15,medianchain)
}
one_shot_time_average15
one_shot_plot_data <- data.frame(mean=c(one_shot_time_average1,
                                        one_shot_time_average7,
                                        one_shot_time_average15),
                                 label=c(rep("grid-1",length(times)),
                                         rep("grid-7",length(times)),
                                         rep("grid-15",length(times))),
                                 frac=c(rep(times,3)))

jpeg("one_shot_time.jpg", units="in", width=8, height=5, res=300)

ggplot(one_shot_plot_data, aes(frac, mean, fill=factor(label), colour=factor(label))) +
  # geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  # geom_point(aes(x=0.05,y=500),colour="red")+
  labs(x="fraction", y="Seconds", colour="method", fill="method") +
  theme_bw()+ggtitle("Computation time of MCMC chain on different coresets")
