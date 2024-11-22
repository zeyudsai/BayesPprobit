rm(list = ls())
dev.off()


options(scipen=100, digits=8)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("Functions.R")



d<- fread("/Real_world_data/heart_disease_health_indicators_BRFSS2015.csv")

table(d$HeartDiseaseorAttack)

#sample_index <- sample(1:nrow(d),nrow(d)/10,replace = FALSE)
#d <- d[sample_index,]
D <- ncol(d)
N <- nrow(d)
X <- d[,2:22]
d <- d[,-23]
#scale.d <- scale(X)

#scale.d <- data.frame(scale.d,"Spam"= d[,34])
#scale.d
onevector <- rep(1,N)

#X <- as.matrix(cbind(onevector,scale.d))
X <- as.matrix(cbind(onevector,X))
y <- d$HeartDiseaseorAttack

mle_result <- glm(HeartDiseaseorAttack~.,data=d,family = binomial(link = "logit"))
MLE_coefficients <- mle_result$coefficients

mle_probit <- glm(HeartDiseaseorAttack~.,data=d,family = binomial(link = "probit"))
probit_coefficients <- mle_probit$coefficients


N_sim <- 1000
burn_in <- N_sim*0.9
heartattack_2 <- Lp_gibbssampler(N_sim=N_sim,
                               burn_in=burn_in, X=X, y=y,true_theta = probit_coefficients,
                               Lp=2,initial_theda = probit_coefficients)
heartattack_2
heartattack_1 <- Lp_gibbssampler(N_sim=N_sim,
                                burn_in=burn_in, X=X, y=y,true_theta = probit_coefficients,
                                Lp=1,initial_theda = MLE_coefficients)
sample3 <- multi_chain(N_sim=N_sim,burn_in=burn_in,
                       X=X, y=y,initial_theda=probit_coefficients,
                       true_theta=probit_coefficients,Lp=3,M_iter=300,range=c(0.1,5),step=0.02,times = 1)

sample3ci<-apply(sample3$Beta_chain[[1]], MARGIN =2, FUN=ci)
sample3_ci_low <- c()
sample3_ci_high <- c()
for (i in 1:D) {
  sample3_ci_low[i] <- sample3ci[[i]]$CI_low
  sample3_ci_high[i] <- sample3ci[[i]]$CI_high
  
}
sample1_ci <-apply(heartattack_1$chain, MARGIN =2, FUN=ci)
sample1_ci_low <- c()
sample1_ci_high <- c()
for (i in 1:D) {
  sample1_ci_low[i] <- sample1_ci[[i]]$CI_low
  sample1_ci_high[i] <- sample1_ci[[i]]$CI_high
  
}
sample2_ci <-apply(heartattack_2$chain, MARGIN =2, FUN=ci)
sample2_ci_low <- c()
sample2_ci_high <- c()
for (i in 1:D) {
  sample2_ci_low[i] <- sample2_ci[[i]]$CI_low
  sample2_ci_high[i] <- sample2_ci[[i]]$CI_high
  
}
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

ci_plot1 <- data.frame("variable"=c(1:11),
                       "sample3_ci_low"=c(sample1_ci_low[1:11],sample2_ci_low[1:11],
                                          logit_confidence_low[1:11],probit_confidence_low[1:11]),
                       "sample3_ci_high"=c(sample1_ci_high[1:11],sample2_ci_high[1:11],
                                           logit_confidence_high[1:11],probit_confidence_high[1:11]),
                       "post_mean"=c(heartattack_1$post.mean[2:12],heartattack_2$post.mean[2:12],
                                     mle_result$coefficients[2:12],mle_probit$coefficients[2:12]),
                       "Method"=c(rep("Bayesian Credible Interval (P=1)",11),
                                  rep("Bayesian Credible Interval (P=2)",11),
                                  rep("Logit Confidence Interval",11),
                                  rep("Probit Confidence Interval",11)))
ci_plot2 <- data.frame("variable"=c(12:21),
                       "sample3_ci_low"=c(sample1_ci_low[12:21],sample2_ci_low[12:21],
                                          logit_confidence_low[12:21],probit_confidence_low[12:21]),
                       "sample3_ci_high"=c(sample1_ci_high[12:21],sample2_ci_high[12:21],
                                           logit_confidence_high[12:21],probit_confidence_high[12:21]),
                       "post_mean"=c(heartattack_1$post.mean[13:22],heartattack_2$post.mean[13:22],
                                     mle_result$coefficients[13:22],mle_probit$coefficients[13:22]),
                       "Method"=c(rep("Bayesian Credible Interval (P=1)",10),
                                  rep("Bayesian Credible Interval (P=2)",10),
                                  rep("Logit Confidence Interval",10),
                                  rep("Probit Confidence Interval",10)))




# jpeg("heartattack_ci_1(version2).jpeg", units="in", width=8, height=5, res=300)
ggplot(ci_plot1, aes(variable, post_mean, colour=Method)) +
  geom_point(position=position_dodge(1)) + 
  geom_errorbar(aes(ymin=sample3_ci_low, ymax=sample3_ci_high), width=0.5, position=position_dodge(1)) + 
  theme_bw()+
  labs(x=TeX("$\\hat{\\beta}\\ index"), y=TeX("$\\hat{\\beta}"), colour="Method", fill="Method")+ 
  # ggtitle("Confidence/Credible Interval of Coefficients \nHeart Disease Data") +
  theme(legend.position = "none")
dev.off()

# jpeg("heartattack_ci_2(version2).jpeg", units="in", width=8, height=5, res=300)
ggplot(ci_plot2, aes(variable, post_mean, colour=Method)) +
  geom_point(position=position_dodge(1)) + 
  geom_errorbar(aes(ymin=sample3_ci_low, ymax=sample3_ci_high), width=1, position=position_dodge(1)) + 
  theme_bw()+
  labs(x=TeX("$\\hat{\\beta}\\ index"), y=TeX("$\\hat{\\beta}"), colour="Method", fill="Method")
  # +ggtitle("Confidence/Credible Interval of Coefficients \nHeart Disease Data") 
dev.off()

# jpeg("heartattack_p_estimation(version2).jpeg", units="in", width=8, height=5, res=300)
ggplot(data=data.frame(Lp.mean.chain=sample3$Lp.mean.chain[seq(250,1000,5)],x=seq(250,1000,5)), aes(x=x, y=Lp.mean.chain)) +
  geom_line()+
  theme_bw(base_size = 12) +
  ylim(2,3)+
  labs(title="Estimation of p for heart disease data",  
       # subtitle="Estimated p=2.21", 
       x="MCMC iteration",
       y="p")  +
annotate("text", x=800, y=2.4,size=6, label= TeX("$\\hat{p} = 2.212", output='character'),parse=TRUE)

# insert ggplot code
dev.off()



plot(1:(D-1),probit_coefficients[-1],type = "l",ylim = c(-1,1),main = "comparison of coefficients",xlab = "Beta-index",
     ylab = "parameters")
Msubtitle <- "Heart Disease data"
mtext(side = 3, line = 0.4, Msubtitle)
legend("topright",legend=c("MLE Probit","Bayesian-1-probit","Bayesian-2-probit"), col=c(1,2,3), 
       bty = 'n', lwd = 2, inset = c(0, 0), lty = 1, cex = 0.73)
lines(heartattack_1$post.mean[-1],col=2)
lines(heartattack_2$post.mean[-1],col=3)

plot(1:(D-1),probit_coefficients[-1],type = "l",ylim = c(-1,1),main = "comparison of coefficients",xlab = "Beta-index",
     ylab = "parameters")
Msubtitle <- "Heart Disease data"
mtext(side = 3, line = 0.4, Msubtitle)
legend("topright",legend=c("MLE Probit","Bayesian-2.26(estimated)-probit"), col=c(1,2), 
       bty = 'n', lwd = 2, inset = c(0, 0), lty = 1, cex = 0.73)
lines(sample3$Beta.mean[-1],col=2)

plot(1:(D-1),MLE_coefficients[-1],type = "l",ylim = c(-1,1.5),main = "comparison of coefficients",xlab = "Beta-index",
     ylab = "parameters")
Msubtitle <- "Heart Disease data"
mtext(side = 3, line = 0.4, Msubtitle)
legend("topright",legend=c("MLE Logit","Bayesian-1-probit","Bayesian-2-probit"), col=c(1,2,3), 
       bty = 'n', lwd = 2, inset = c(0, 0), lty = 1, cex = 0.73)
lines(heartattack_1$post.mean[-1],col=2)
lines(heartattack_2$post.mean[-1],col=3)

plot(1:(D-1),MLE_coefficients[-1],type = "l",ylim = c(-1,1.5),main = "comparison of coefficients",xlab = "Beta-index",
     ylab = "parameters")
Msubtitle <- "Heart Disease data"
mtext(side = 3, line = 0.4, Msubtitle)
legend("topright",legend=c("MLE Logit","Bayesian-2.26(estimated)-probit"), col=c(1,2), 
       bty = 'n', lwd = 2, inset = c(0, 0), lty = 1, cex = 0.73)
lines(sample3$Beta.mean[-1],col=2)



plot(10:N_sim, sample3$Lp.mean.chain[10:2000],ylim=c(2,2.5),type = "l",xlab = "MCMC iteration",ylab = "estimated P",
     main = "The estimated P of the heart disease data")

Msubtitle <- "P is estimated as 2.26"
mtext(side = 3, line = 0.4, Msubtitle)


sketch_size <- D^2

frac <- seq(0.0001,0.005,0.0001)


########fixed p=2#########
#####algorithm 1#####
for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Compute_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=2)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y[sample_index],
                             initial_theda=probit_coefficients,
                             true_theta=probit_coefficients,Lp=2)
    assign(paste("heartattack_p_2_frac",frac[i],"chain_",j,sep = ""), Model)
  }
}

#####LP sampling#########
for (i in 1:length(frac)) {
  for (j in 1:5) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Lp_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=2)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y[sample_index],
                             initial_theda=probit_coefficients,
                             true_theta=probit_coefficients,Lp=2)
    assign(paste("heartattack_p_2_lpcoreset_frac",frac[i],"chain_",j,sep = ""), Model)
  }
}
#########uniform sampling
for (i in 1:length(frac)) {
  for (j in 1:5) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Uniform_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=2)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y[sample_index],
                             initial_theda=probit_coefficients,
                             true_theta=probit_coefficients,Lp=2)
    assign(paste("heartattack_p_2_uniform_frac",frac[i],"chain_",j,sep = ""), Model)
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
    assign(paste("heartattack_p_2_one_shot_frac",frac[i],"chain_",j,sep = ""), Model)
  }
}
#############post mean difference######


heartattack_p_2_one_shot_mean_diff_norm_median <- c()
heartattack_p_2_one_shot_mean_diff_norm_max <- c()
heartattack_p_2_one_shot_mean_diff_norm_min <- c()
heartattack_p_2_one_shot_mean_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    b <- norm(get(paste("heartattack_p_2_one_shot_frac",frac[i],"chain_",j,sep = ""))$post.mean-heartattack_2$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  heartattack_p_2_one_shot_mean_diff_norm_median <- c(heartattack_p_2_one_shot_mean_diff_norm_median,median_result)
  heartattack_p_2_one_shot_mean_diff_norm_max <- c(heartattack_p_2_one_shot_mean_diff_norm_max,maxresult)
  heartattack_p_2_one_shot_mean_diff_norm_min <- c(heartattack_p_2_one_shot_mean_diff_norm_min,minresult)
  heartattack_p_2_one_shot_mean_diff_norm_average <- c(heartattack_p_2_one_shot_mean_diff_norm_average,meanresult)
}
heartattack_p_2_frac_mean_diff_norm_median <- c()
heartattack_p_2_frac_mean_diff_norm_max <- c()
heartattack_p_2_frac_mean_diff_norm_min <- c()
heartattack_p_2_frac_mean_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    b <- norm(get(paste("heartattack_p_2_frac",frac[i],"chain_",j,sep = ""))$post.mean-heartattack$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  heartattack_p_2_frac_mean_diff_norm_median <- c(heartattack_p_2_frac_mean_diff_norm_median,median_result)
  heartattack_p_2_frac_mean_diff_norm_max <- c(heartattack_p_2_frac_mean_diff_norm_max,maxresult)
  heartattack_p_2_frac_mean_diff_norm_min <- c(heartattack_p_2_frac_mean_diff_norm_min,minresult)
  heartattack_p_2_frac_mean_diff_norm_average <- c(heartattack_p_2_frac_mean_diff_norm_average,meanresult)
}

heartattack_p_2_lpcoreset_frac_mean_diff_norm_median <- c()
heartattack_p_2_lpcoreset_frac_mean_diff_norm_max <- c()
heartattack_p_2_lpcoreset_frac_mean_diff_norm_min <- c()
heartattack_p_2_lpcoreset_frac_mean_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    b <- norm(get(paste("heartattack_p_2_lpcoreset_frac",frac[i],"chain_",j,sep = ""))$post.mean-heartattack$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  heartattack_p_2_lpcoreset_frac_mean_diff_norm_median <- c(heartattack_p_2_lpcoreset_frac_mean_diff_norm_median,median_result)
  heartattack_p_2_lpcoreset_frac_mean_diff_norm_max <- c(heartattack_p_2_lpcoreset_frac_mean_diff_norm_max,maxresult)
  heartattack_p_2_lpcoreset_frac_mean_diff_norm_min <- c(heartattack_p_2_lpcoreset_frac_mean_diff_norm_min,minresult)
  heartattack_p_2_lpcoreset_frac_mean_diff_norm_average <- c(heartattack_p_2_lpcoreset_frac_mean_diff_norm_average,meanresult)
  
}
heartattack_p_2_uniform_frac_mean_diff_norm_median <- c()
heartattack_p_2_uniform_frac_mean_diff_norm_max <- c()
heartattack_p_2_uniform_frac_mean_diff_norm_min <- c()
heartattack_p_2_uniform_frac_mean_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    b <- norm(get(paste("heartattack_p_2_uniform_frac",frac[i],"chain_",j,sep = ""))$post.mean-heartattack$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  heartattack_p_2_uniform_frac_mean_diff_norm_median <- c(heartattack_p_2_uniform_frac_mean_diff_norm_median,median_result)
  heartattack_p_2_uniform_frac_mean_diff_norm_max <- c(heartattack_p_2_uniform_frac_mean_diff_norm_max,maxresult)
  heartattack_p_2_uniform_frac_mean_diff_norm_min <- c(heartattack_p_2_uniform_frac_mean_diff_norm_min,minresult)
  heartattack_p_2_uniform_frac_mean_diff_norm_average <- c(heartattack_p_2_uniform_frac_mean_diff_norm_average,meanresult)
  
}

p2_meandiff_heartattack_data <- data.frame(mean=c(heartattack_p_2_frac_mean_diff_norm_average,
                                                  heartattack_p_2_one_shot_mean_diff_norm_median,
                                                heartattack_p_2_uniform_frac_mean_diff_norm_average),
                                         median=c(heartattack_p_2_frac_mean_diff_norm_median,
                                                  heartattack_p_2_one_shot_mean_diff_norm_median,
                                                  heartattack_p_2_uniform_frac_mean_diff_norm_median),
                                         max=c(heartattack_p_2_frac_mean_diff_norm_max,
                                               heartattack_p_2_one_shot_mean_diff_norm_median,
                                               heartattack_p_2_uniform_frac_mean_diff_norm_max),
                                         min=c(heartattack_p_2_frac_mean_diff_norm_min,
                                               heartattack_p_2_one_shot_mean_diff_norm_median,
                                               heartattack_p_2_uniform_frac_mean_diff_norm_min),
                                         label=c(rep("2-probit",length(frac)),rep("one-shot",length(frac)),rep("uniform",length(frac))),
                                         frac=c(rep(frac,3)))
jpeg("heartdisease_p2_mean.jpeg", units="in", width=8, height=5, res=300)
ggplot(p2_meandiff_heartattack_data, aes(frac, mean, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior mean", colour="method", fill="method") +
  theme_bw()+
  ggtitle("Heart disease data (p=2)")
dev.off()
######poserior covariance difference
heartattack_p2_pos_cov_diff_norm_median <- c()
heartattack_p2_pos_cov_diff_norm_max <- c()
heartattack_p2_pos_cov_diff_norm_min <- c()
heartattack_p2_pos_cov_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("heartattack_p_2_frac",frac[i],"chain_",j,sep = ""))$chain)-
                cov(heartattack$chain),type = "2")
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN = 2, mean)
  heartattack_p2_pos_cov_diff_norm_median <- c(heartattack_p2_pos_cov_diff_norm_median,median_result)
  heartattack_p2_pos_cov_diff_norm_max <- c(heartattack_p2_pos_cov_diff_norm_max,maxresult)
  heartattack_p2_pos_cov_diff_norm_min <- c(heartattack_p2_pos_cov_diff_norm_min,minresult)
  heartattack_p2_pos_cov_diff_norm_average <- c(heartattack_p2_pos_cov_diff_norm_average,averageresult)
  
  
}

heartattack_p_2_one_shot_cov_diff_norm_median <- c()
heartattack_p_2_one_shot_cov_diff_norm_max <- c()
heartattack_p_2_one_shot_cov_diff_norm_min <- c()
heartattack_p_2_one_shot_cov_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("heartattack_p_2_one_shot_frac",frac[i],"chain_",j,sep = ""))$chain)-
                cov(heartattack_2$chain),type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  heartattack_p_2_one_shot_cov_diff_norm_median <- c(heartattack_p_2_one_shot_cov_diff_norm_median,median_result)
  heartattack_p_2_one_shot_cov_diff_norm_max <- c(heartattack_p_2_one_shot_cov_diff_norm_max,maxresult)
  heartattack_p_2_one_shot_cov_diff_norm_min <- c(heartattack_p_2_one_shot_cov_diff_norm_min,minresult)
  heartattack_p_2_one_shot_cov_diff_norm_average <- c(heartattack_p_2_one_shot_cov_diff_norm_average,meanresult)
}

heartattack_p2_lpcoreset_pos_cov_diff_norm_median <- c()
heartattack_p2_lpcoreset_pos_cov_diff_norm_max <- c()
heartattack_p2_lpcoreset_pos_cov_diff_norm_min <- c()
heartattack_p2_lpcoreset_pos_cov_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("heartattack_p_2_lpcoreset_frac",frac[i],"chain_",j,sep = ""))$chain)-
                cov(heartattack$chain),type = "2")
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN = 2, mean)
  heartattack_p2_lpcoreset_pos_cov_diff_norm_median <- c(heartattack_p2_lpcoreset_pos_cov_diff_norm_median,median_result)
  heartattack_p2_lpcoreset_pos_cov_diff_norm_max <- c(heartattack_p2_lpcoreset_pos_cov_diff_norm_max,maxresult)
  heartattack_p2_lpcoreset_pos_cov_diff_norm_min <- c(heartattack_p2_lpcoreset_pos_cov_diff_norm_min,minresult)
  heartattack_p2_lpcoreset_pos_cov_diff_norm_average <- c(heartattack_p2_lpcoreset_pos_cov_diff_norm_average,averageresult)
  
}

heartattack_p2_uniform_pos_cov_diff_norm_median <- c()
heartattack_p2_uniform_pos_cov_diff_norm_max <- c()
heartattack_p2_uniform_pos_cov_diff_norm_min <- c()
heartattack_p2_uniform_pos_cov_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("heartattack_p_2_uniform_frac",frac[i],"chain_",j,sep = ""))$chain)-
                cov(heartattack$chain),type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN = 2, mean)
  heartattack_p2_uniform_pos_cov_diff_norm_median <- c(heartattack_p2_uniform_pos_cov_diff_norm_median,median_result)
  heartattack_p2_uniform_pos_cov_diff_norm_max <- c(heartattack_p2_uniform_pos_cov_diff_norm_max,maxresult)
  heartattack_p2_uniform_pos_cov_diff_norm_min <- c(heartattack_p2_uniform_pos_cov_diff_norm_min,minresult)
  heartattack_p2_uniform_pos_cov_diff_norm_average <- c(heartattack_p2_uniform_pos_cov_diff_norm_average,averageresult)
}

p2_covdiff_heartattack_data <- data.frame(mean=c(heartattack_p2_pos_cov_diff_norm_average,
                                                 heartattack_p_2_one_shot_cov_diff_norm_average,
                                               heartattack_p2_uniform_pos_cov_diff_norm_average),
                                        median=c(heartattack_p2_pos_cov_diff_norm_median,
                                                 heartattack_p_2_one_shot_cov_diff_norm_median,
                                                 heartattack_p2_uniform_pos_cov_diff_norm_median),
                                        max=c(heartattack_p2_pos_cov_diff_norm_max,
                                              heartattack_p_2_one_shot_cov_diff_norm_max,
                                              heartattack_p2_uniform_pos_cov_diff_norm_max),
                                        min=c(heartattack_p2_pos_cov_diff_norm_min,
                                              heartattack_p_2_one_shot_cov_diff_norm_min,
                                              heartattack_p2_uniform_pos_cov_diff_norm_min),
                                        label=c(rep("2-probit",length(frac)),rep("one-shot",length(frac)),rep("uniform",length(frac))),
                                        frac=c(rep(frac,3)))
p2_covdiff_heartattack_data

jpeg("heartdisease_p2_cov.jpeg", units="in", width=8, height=5, res=300)
ggplot(p2_covdiff_heartattack_data, aes(frac, mean, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior covariance", colour="method", fill="method") +
  theme_bw()+
  ggtitle("Heart disease data (p=2)")
dev.off()
  
# 
# #####MMD#######
# heartattack_p_2_MMD_median <- c()
# heartattack_p_2_MMD_max <- c()
# heartattack_p_2_MMD_min <- c()
# heartattack_p_2_MMD_average <- c()
# 
# for (i in 1:length(frac)){
#   print(paste("2-probit_MMD_Frac:",i,sep = ""))
#   medianchain <- c()
#   for (j in 1:10) {
#     b <- mmd(sample_x=tail(get(paste("heartattack_p_2_frac",frac[i],"chain_",j,sep = ""))$chain,100),
#              sample_y = tail(cover_type_2$chain,100))
#     medianchain<- rbind(medianchain,b)
#   }
#   median_result <- apply(medianchain, MARGIN = 2, median)
#   maxresult <- apply(medianchain, MARGIN = 2, max)
#   minresult <- apply(medianchain, MARGIN = 2, min)
#   averageresult <- apply(medianchain, MARGIN=2, mean)
#   heartattack_p_2_MMD_median <- c(heartattack_p_2_MMD_median,median_result)
#   heartattack_p_2_MMD_max <- c(heartattack_p_2_MMD_max,maxresult)
#   heartattack_p_2_MMD_min <- c(heartattack_p_2_MMD_min,minresult)
#   heartattack_p_2_MMD_average <- c(heartattack_p_2_MMD_average,averageresult)
# }
# 
# heartattack_p_2_lpcoreset_MMD_median <- c()
# heartattack_p_2_lpcoreset_MMD_max <- c()
# heartattack_p_2_lpcoreset_MMD_min <- c()
# heartattack_p_2_lpcoreset_MMD_average <- c()
# for (i in 1:length(frac)){
#   print(paste("rootl2_MMD_Frac:",i,sep = ""))
#   medianchain <- c()
#   for (j in 1:10) {
#     b <- mmd(sample_x=tail(get(paste("heartattack_p_2_lpcoreset_frac",frac[i],"chain_",j,sep = ""))$chain,100),
#              sample_y = tail(cover_type_2$chain,100))
#     medianchain<- rbind(medianchain,b)
#   }
#   median_result <- apply(medianchain, MARGIN = 2, median)
#   maxresult <- apply(medianchain, MARGIN = 2, max)
#   minresult <- apply(medianchain, MARGIN = 2, min)
#   averageresult <- apply(medianchain, MARGIN=2, mean)
#   heartattack_p_2_lpcoreset_MMD_median <- c(heartattack_p_2_lpcoreset_MMD_median,median_result)
#   heartattack_p_2_lpcoreset_MMD_max <- c(heartattack_p_2_lpcoreset_MMD_max,maxresult)
#   heartattack_p_2_lpcoreset_MMD_min <- c(heartattack_p_2_lpcoreset_MMD_min,minresult)
#   heartattack_p_2_lpcoreset_MMD_average <- c(heartattack_p_2_lpcoreset_MMD_average,averageresult)
# }
# 
# heartattack_p_2_uniform_MMD_median <- c()
# heartattack_p_2_uniform_MMD_max <- c()
# heartattack_p_2_uniform_MMD_min <- c()
# heartattack_p_2_uniform_MMD_average <- c()
# for (i in 1:length(frac)){
#   print(paste("uniform_MMD_Frac:",i,sep = ""))
#   medianchain <- c()
#   for (j in 1:5) {
#     b <- mmd(sample_x=tail(get(paste("heartattack_p_2_uniform_frac",frac[i],"chain_",j,sep = ""))$chain,100),
#              sample_y = tail(cover_type_2$chain,100))
#     medianchain<- rbind(medianchain,b)
#   }
#   median_result <- apply(medianchain, MARGIN = 2, median)
#   maxresult <- apply(medianchain, MARGIN = 2, max)
#   minresult <- apply(medianchain, MARGIN = 2, min)
#   averageresult <- apply(medianchain, MARGIN=2, mean)
#   heartattack_p_2_uniform_MMD_median <- c(heartattack_p_2_uniform_MMD_median,median_result)
#   heartattack_p_2_uniform_MMD_max <- c(heartattack_p_2_uniform_MMD_max,maxresult)
#   heartattack_p_2_uniform_MMD_min <- c(heartattack_p_2_uniform_MMD_min,minresult)
#   heartattack_p_2_uniform_MMD_average <- c(heartattack_p_2_uniform_MMD_average,averageresult)
# }
# 
# p2_mmd_heartattack_data <- data.frame(mean=c(heartattack_p_2_MMD_average,
#                                            heartattack_p_2_lpcoreset_MMD_average,
#                                            heartattack_p_2_uniform_MMD_average),
#                                     median=c(heartattack_p_2_MMD_median,
#                                              heartattack_p_2_lpcoreset_MMD_median,
#                                              heartattack_p_2_uniform_MMD_median),
#                                     max=c(heartattack_p_2_MMD_max,
#                                           heartattack_p_2_lpcoreset_MMD_max,
#                                           heartattack_p_2_uniform_MMD_max),
#                                     min=c(heartattack_p_2_MMD_min,
#                                           heartattack_p_2_lpcoreset_MMD_min,
#                                           heartattack_p_2_uniform_MMD_min),
#                                     label=c(rep("2-probit",length(frac)),rep("root l2",length(frac)),rep("uniform",length(frac))),
#                                     frac=c(rep(frac,3)))
# p2_mmd_heartattack_data
# ggplot(p2_mmd_heartattack_data, aes(frac, mean, fill=factor(label), colour=factor(label))) +
#   geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
#   geom_line() +
#   labs(x="fraction", y="MMD", colour="method", fill="method") +
#   theme_bw()
###########
###########


############lilelihood ratio###############


heartattack_p2_2probit_llk_ratio_median <- c()
heartattack_p2_2probit_llk_ratio_max <- c()
heartattack_p2_2probit_llk_ratio_min <- c()
heartattack_p2_2probit_llk_ratio_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    # b <- sum(llk(X=a,get(paste("heartattack_p_2_frac",frac[i],"chain_",j,sep = ""))$post.mean,2))/sum(llk(X=a, heartattack_2$post.mean,2))
    
    b<- median(abs(llk(X, heartattack_2$post.mean,2)-
                     llk(X,get(paste("heartattack_p_2_frac",frac[i],"chain_",j,sep = ""))$post.mean,2)))
    
     medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  heartattack_p2_2probit_llk_ratio_median <- c(heartattack_p2_2probit_llk_ratio_median,median_result)
  heartattack_p2_2probit_llk_ratio_max <- c(heartattack_p2_2probit_llk_ratio_max,maxresult)
  heartattack_p2_2probit_llk_ratio_min <- c(heartattack_p2_2probit_llk_ratio_min,minresult)
  heartattack_p2_2probit_llk_ratio_average <- c(heartattack_p2_2probit_llk_ratio_average,meanresult)
}

heartattack_p2_uniform_llk_ratio_median <- c()
heartattack_p2_uniform_llk_ratio_max <- c()
heartattack_p2_uniform_llk_ratio_min <- c()
heartattack_p2_uniform_llk_ratio_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    # b <- sum(llk(a,get(paste("heartattack_p_2_uniform_frac",frac[i],"chain_",j,sep = ""))$post.mean,2))/sum(llk(a, heartattack_2$post.mean,2))
    b<- median(abs(llk(X, heartattack_2$post.mean,2)-
                     llk(X,get(paste("heartattack_p_2_uniform_frac",frac[i],"chain_",j,sep = ""))$post.mean,2)))
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  heartattack_p2_uniform_llk_ratio_median <- c(heartattack_p2_uniform_llk_ratio_median,median_result)
  heartattack_p2_uniform_llk_ratio_max <- c(heartattack_p2_uniform_llk_ratio_max,maxresult)
  heartattack_p2_uniform_llk_ratio_min <- c(heartattack_p2_uniform_llk_ratio_min,minresult)
  heartattack_p2_uniform_llk_ratio_average <- c(heartattack_p2_uniform_llk_ratio_average,meanresult)
}

heartattack_p2_one_shot_llk_ratio_median <- c()
heartattack_p2_one_shot_llk_ratio_max <- c()
heartattack_p2_one_shot_llk_ratio_min <- c()
heartattack_p2_one_shot_llk_ratio_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    # b <- sum(llk(a,get(paste("heartattack_p_2_one_shot_frac",frac[i],"chain_",j,sep = ""))$post.mean,2))/sum(llk(a, heartattack_2$post.mean,2))
    b<- median(abs(llk(X, heartattack_2$post.mean,2)-
                     llk(X,get(paste("heartattack_p_2_one_shot_frac",frac[i],"chain_",j,sep = ""))$post.mean,2)))
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  heartattack_p2_one_shot_llk_ratio_median <- c(heartattack_p2_one_shot_llk_ratio_median,median_result)
  heartattack_p2_one_shot_llk_ratio_max <- c(heartattack_p2_one_shot_llk_ratio_max,maxresult)
  heartattack_p2_one_shot_llk_ratio_min <- c(heartattack_p2_one_shot_llk_ratio_min,minresult)
  heartattack_p2_one_shot_llk_ratio_average <- c(heartattack_p2_one_shot_llk_ratio_average,meanresult)
}

p2_llk_ratio_heartattack_data <- data.frame(mean=c(log(heartattack_p2_2probit_llk_ratio_average),
                                              log(heartattack_p2_one_shot_llk_ratio_average),
                                              log(heartattack_p2_uniform_llk_ratio_average)),
                                       median=c(log(heartattack_p2_2probit_llk_ratio_median),
                                                log(heartattack_p2_one_shot_llk_ratio_median),
                                                log(heartattack_p2_uniform_llk_ratio_median)),
                                       max=c(log(heartattack_p2_2probit_llk_ratio_max),
                                             log(heartattack_p2_one_shot_llk_ratio_max),
                                             log(heartattack_p2_uniform_llk_ratio_max)),
                                       min=c(log(heartattack_p2_2probit_llk_ratio_min),
                                             log(heartattack_p2_one_shot_llk_ratio_min),
                                             log(heartattack_p2_uniform_llk_ratio_min)),
                                       label=c(rep("2-probit",length(frac)),rep("one-shot",length(frac)),rep("uniform",length(frac))),
                                       frac=c(rep(frac,3)))

p2_llk_ratio_heartattack_data <- data.frame(mean=c(heartattack_p2_2probit_llk_ratio_average,
                                              heartattack_p2_one_shot_llk_ratio_average,
                                              heartattack_p2_uniform_llk_ratio_average),
                                       median=c(heartattack_p2_2probit_llk_ratio_median,
                                                heartattack_p2_one_shot_llk_ratio_median,
                                                heartattack_p2_uniform_llk_ratio_median),
                                       max=c(heartattack_p2_2probit_llk_ratio_max,
                                             heartattack_p2_one_shot_llk_ratio_max,
                                             heartattack_p2_uniform_llk_ratio_max),
                                       min=c(heartattack_p2_2probit_llk_ratio_min,
                                             heartattack_p2_one_shot_llk_ratio_min,
                                             heartattack_p2_uniform_llk_ratio_min),
                                       label=c(rep("2-probit",length(frac)),rep("one-shot",length(frac)),rep("uniform",length(frac))),
                                       frac=c(rep(frac,3)))

ggplot(p2_llk_ratio_heartattack_data, aes(frac, mean, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="(log) Bayes Factor", colour="method", fill="method") +
  theme_bw()+
  ggtitle("(log)-Bayes Factor of heartattackcard data (p=2)")
dev.off()
########fixed p=1#########
for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Compute_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=1)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y[sample_index],
                             initial_theda=probit_coefficients,
                             true_theta=probit_coefficients,Lp=1)
    assign(paste("heartattack_p_1_frac",frac[i],"chain_",j,sep = ""), Model)
  }
}

#####LP sampling#########
for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Lp_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=1)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y[sample_index],
                             initial_theda=probit_coefficients,
                             true_theta=probit_coefficients,Lp=1)
    assign(paste("heartattack_p_1_lpcoreset_frac",frac[i],"chain_",j,sep = ""), Model)
  }
}
#########uniform sampling##############
for (i in 1:length(frac)) {
  for (j in 1:10) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Uniform_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=1)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y[sample_index],
                             initial_theda=probit_coefficients,
                             true_theta=probit_coefficients,Lp=1)
    assign(paste("heartattack_p_1_uniform_frac",frac[i],"chain_",j,sep = ""), Model)
  }
}

###############one-shot##########
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
    assign(paste("heartattack_p_1_one_shot_frac",frac[i],"chain_",j,sep = ""), Model)
  }
}
########################################
############lilelihood ratio###############


heartattack_p1_1probit_llk_ratio_median <- c()
heartattack_p1_1probit_llk_ratio_max <- c()
heartattack_p1_1probit_llk_ratio_min <- c()
heartattack_p1_1probit_llk_ratio_average <- c()
sum(llk(X, heartattack_2$post.mean,1))
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    # b <- sum(llk(X=a,get(paste("heartattack_p_1_frac",frac[i],"chain_",j,sep = ""))$post.mean,1))/sum(llk(X=a, heartattack_1$post.mean,1))
    b<- median(abs(llk(X, heartattack_1$post.mean,1,y)-
                     llk(X,get(paste("heartattack_p_1_frac",frac[i],"chain_",j,sep = ""))$post.mean,1,y)))
    
     medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  heartattack_p1_1probit_llk_ratio_median <- c(heartattack_p1_1probit_llk_ratio_median,median_result)
  heartattack_p1_1probit_llk_ratio_max <- c(heartattack_p1_1probit_llk_ratio_max,maxresult)
  heartattack_p1_1probit_llk_ratio_min <- c(heartattack_p1_1probit_llk_ratio_min,minresult)
  heartattack_p1_1probit_llk_ratio_average <- c(heartattack_p1_1probit_llk_ratio_average,meanresult)
}

heartattack_p1_uniform_llk_ratio_median <- c()
heartattack_p1_uniform_llk_ratio_max <- c()
heartattack_p1_uniform_llk_ratio_min <- c()
heartattack_p1_uniform_llk_ratio_average <- c()


for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    # b <- sum(llk(a,get(paste("heartattack_p_1_uniform_frac",frac[i],"chain_",j,sep = ""))$post.mean,1))/sum(llk(a, heartattack_1$post.mean,1))
    b<- median(abs(llk(X, heartattack_1$post.mean,1,y)-
                     llk(X,get(paste("heartattack_p_1_uniform_frac",frac[i],"chain_",j,sep = ""))$post.mean,1,y)))
    
     medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  heartattack_p1_uniform_llk_ratio_median <- c(heartattack_p1_uniform_llk_ratio_median,median_result)
  heartattack_p1_uniform_llk_ratio_max <- c(heartattack_p1_uniform_llk_ratio_max,maxresult)
  heartattack_p1_uniform_llk_ratio_min <- c(heartattack_p1_uniform_llk_ratio_min,minresult)
  heartattack_p1_uniform_llk_ratio_average <- c(heartattack_p1_uniform_llk_ratio_average,meanresult)
}

heartattack_p1_one_shot_llk_ratio_median <- c()
heartattack_p1_one_shot_llk_ratio_max <- c()
heartattack_p1_one_shot_llk_ratio_min <- c()
heartattack_p1_one_shot_llk_ratio_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    # b <- sum(llk(a,get(paste("heartattack_p_1_one_shot_frac",frac[i],"chain_",j,sep = ""))$post.mean,1))/sum(llk(a, heartattack_1$post.mean,1))
    b<- median(abs(llk(X, heartattack_1$post.mean,1,y)-
                     llk(X,get(paste("heartattack_p_1_one_shot_frac",frac[i],"chain_",j,sep = ""))$post.mean,1,y)))
    
     medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  heartattack_p1_one_shot_llk_ratio_median <- c(heartattack_p1_one_shot_llk_ratio_median,median_result)
  heartattack_p1_one_shot_llk_ratio_max <- c(heartattack_p1_one_shot_llk_ratio_max,maxresult)
  heartattack_p1_one_shot_llk_ratio_min <- c(heartattack_p1_one_shot_llk_ratio_min,minresult)
  heartattack_p1_one_shot_llk_ratio_average <- c(heartattack_p1_one_shot_llk_ratio_average,meanresult)
}

heartattack_p1_rootl2_llk_ratio_median <- c()
heartattack_p1_rootl2_llk_ratio_max <- c()
heartattack_p1_rootl2_llk_ratio_min <- c()
heartattack_p1_rootl2_llk_ratio_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    # b <- sum(llk(a,get(paste("heartattack_p_1_one_shot_frac",frac[i],"chain_",j,sep = ""))$post.mean,1))/sum(llk(a, heartattack_1$post.mean,1))
    b<- median(abs(llk(X, heartattack_1$post.mean,1,y)-
                     llk(X,get(paste("heartattack_p_1_lpcoreset_frac",frac[i],"chain_",j,sep = ""))$post.mean,1,y)))
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  heartattack_p1_rootl2_llk_ratio_median <- c(heartattack_p1_rootl2_llk_ratio_median,median_result)
  heartattack_p1_rootl2_llk_ratio_max <- c(heartattack_p1_rootl2_llk_ratio_max,maxresult)
  heartattack_p1_rootl2_llk_ratio_min <- c(heartattack_p1_rootl2_llk_ratio_min,minresult)
  heartattack_p1_rootl2_llk_ratio_average <- c(heartattack_p1_rootl2_llk_ratio_average,meanresult)
}

p1_llk_ratio_heartattack_data <- data.frame(mean=c(log(heartattack_p1_1probit_llk_ratio_average),
                                                   log(heartattack_p1_one_shot_llk_ratio_average),
                                                   log(heartattack_p1_uniform_llk_ratio_average)),
                                            median=c(log(heartattack_p1_1probit_llk_ratio_median),
                                                     log(heartattack_p1_one_shot_llk_ratio_median),
                                                     log(heartattack_p1_uniform_llk_ratio_median)),
                                            max=c(log(heartattack_p1_1probit_llk_ratio_max),
                                                  log(heartattack_p1_one_shot_llk_ratio_max),
                                                  log(heartattack_p1_uniform_llk_ratio_max)),
                                            min=c(log(heartattack_p1_1probit_llk_ratio_min),
                                                  log(heartattack_p1_one_shot_llk_ratio_min),
                                                  log(heartattack_p1_uniform_llk_ratio_min)),
                                            label=c(rep("1-probit",length(frac)),rep("one-shot",length(frac)),rep("uniform",length(frac))),
                                            frac=c(rep(frac,3)))

p1_llk_ratio_heartattack_data <- data.frame(mean=c(heartattack_p1_1probit_llk_ratio_average,
                                                   heartattack_p1_one_shot_llk_ratio_average,
                                                   heartattack_p1_uniform_llk_ratio_average,
                                                   heartattack_p1_rootl2_llk_ratio_average),
                                            median=c(heartattack_p1_1probit_llk_ratio_median,
                                                     heartattack_p1_one_shot_llk_ratio_median,
                                                     heartattack_p1_uniform_llk_ratio_median,
                                                     heartattack_p1_rootl2_llk_ratio_median),
                                            max=c(heartattack_p1_1probit_llk_ratio_max,
                                                  heartattack_p1_one_shot_llk_ratio_max,
                                                  heartattack_p1_uniform_llk_ratio_max,
                                                  heartattack_p1_rootl2_llk_ratio_max),
                                            min=c(heartattack_p1_1probit_llk_ratio_min,
                                                  heartattack_p1_one_shot_llk_ratio_min,
                                                  heartattack_p1_uniform_llk_ratio_min,
                                                  heartattack_p1_rootl2_llk_ratio_min),
                                            label=c(rep("1-probit",length(frac)),
                                                    rep("one-shot",length(frac)),
                                                    rep("uniform",length(frac)),
                                                    rep("rootl2",length(frac))),
                                            frac=c(rep(frac,4)))


ggplot(p1_llk_ratio_heartattack_data, aes(frac, mean, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="(log)-Bayes Factor", colour="method", fill="method") +
  theme_bw()+
  ggtitle("(log)-Bayes Factor of heartattackcard data (p=1)")
dev.off()
#############post mean difference######

heartattack_p_1_one_shot_mean_diff_norm_median <- c()
heartattack_p_1_one_shot_mean_diff_norm_max <- c()
heartattack_p_1_one_shot_mean_diff_norm_min <- c()
heartattack_p_1_one_shot_mean_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("heartattack_p_1_one_shot_frac",frac[i],"chain_",j,sep = ""))$post.mean-heartattack_1$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  heartattack_p_1_one_shot_mean_diff_norm_median <- c(heartattack_p_1_one_shot_mean_diff_norm_median,median_result)
  heartattack_p_1_one_shot_mean_diff_norm_max <- c(heartattack_p_1_one_shot_mean_diff_norm_max,maxresult)
  heartattack_p_1_one_shot_mean_diff_norm_min <- c(heartattack_p_1_one_shot_mean_diff_norm_min,minresult)
  heartattack_p_1_one_shot_mean_diff_norm_average <- c(heartattack_p_1_one_shot_mean_diff_norm_average,meanresult)
}
heartattack_p_1_frac_mean_diff_norm_median <- c()
heartattack_p_1_frac_mean_diff_norm_max <- c()
heartattack_p_1_frac_mean_diff_norm_min <- c()
heartattack_p_1_frac_mean_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("heartattack_p_1_frac",frac[i],"chain_",j,sep = ""))$post.mean-heartattack_1$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  heartattack_p_1_frac_mean_diff_norm_median <- c(heartattack_p_1_frac_mean_diff_norm_median,median_result)
  heartattack_p_1_frac_mean_diff_norm_max <- c(heartattack_p_1_frac_mean_diff_norm_max,maxresult)
  heartattack_p_1_frac_mean_diff_norm_min <- c(heartattack_p_1_frac_mean_diff_norm_min,minresult)
  heartattack_p_1_frac_mean_diff_norm_average <- c(heartattack_p_1_frac_mean_diff_norm_average,meanresult)
  
}
##########lpcoreset######
heartattack_p_1_lpcoreset_frac_mean_diff_norm_median <- c()
heartattack_p_1_lpcoreset_frac_mean_diff_norm_max <- c()
heartattack_p_1_lpcoreset_frac_mean_diff_norm_min <- c()
heartattack_p_1_lpcoreset_frac_mean_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("heartattack_p_1_lpcoreset_frac",frac[i],"chain_",j,sep = ""))$post.mean-heartattack_1$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  heartattack_p_1_lpcoreset_frac_mean_diff_norm_median <- c(heartattack_p_1_lpcoreset_frac_mean_diff_norm_median,median_result)
  heartattack_p_1_lpcoreset_frac_mean_diff_norm_max <- c(heartattack_p_1_lpcoreset_frac_mean_diff_norm_max,maxresult)
  heartattack_p_1_lpcoreset_frac_mean_diff_norm_min <- c(heartattack_p_1_lpcoreset_frac_mean_diff_norm_min,minresult)
  heartattack_p_1_lpcoreset_frac_mean_diff_norm_average <- c(heartattack_p_1_lpcoreset_frac_mean_diff_norm_average,meanresult)
  
}
#######uniform sampling###########
heartattack_p_1_uniform_frac_mean_diff_norm_median <- c()
heartattack_p_1_uniform_frac_mean_diff_norm_max <- c()
heartattack_p_1_uniform_frac_mean_diff_norm_min <- c()
heartattack_p_1_uniform_frac_mean_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    
    b <- norm(get(paste("heartattack_p_1_uniform_frac",frac[i],"chain_",j,sep = ""))$post.mean-heartattack_1$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  heartattack_p_1_uniform_frac_mean_diff_norm_median <- c(heartattack_p_1_uniform_frac_mean_diff_norm_median,median_result)
  heartattack_p_1_uniform_frac_mean_diff_norm_max <- c(heartattack_p_1_uniform_frac_mean_diff_norm_max,maxresult)
  heartattack_p_1_uniform_frac_mean_diff_norm_min <- c(heartattack_p_1_uniform_frac_mean_diff_norm_min,minresult)
  heartattack_p_1_uniform_frac_mean_diff_norm_average <- c(heartattack_p_1_uniform_frac_mean_diff_norm_average,meanresult)
  
}

p1_meandiff_heartattack_data <- data.frame(mean=c(log(heartattack_p_1_frac_mean_diff_norm_average),
                                                 log(heartattack_p_1_lpcoreset_frac_mean_diff_norm_average),
                                                 log(heartattack_p_1_uniform_frac_mean_diff_norm_average)),
                                          median=c(log(heartattack_p_1_frac_mean_diff_norm_median),
                                                   log(heartattack_p_1_lpcoreset_frac_mean_diff_norm_median),
                                                   log(heartattack_p_1_uniform_frac_mean_diff_norm_median)),
                                          max=c(log(heartattack_p_1_frac_mean_diff_norm_max),
                                                log(heartattack_p_1_lpcoreset_frac_mean_diff_norm_max),
                                                log(heartattack_p_1_uniform_frac_mean_diff_norm_max)),
                                          min=c(log(heartattack_p_1_frac_mean_diff_norm_min),
                                                log(heartattack_p_1_lpcoreset_frac_mean_diff_norm_min),
                                                log(heartattack_p_1_uniform_frac_mean_diff_norm_min)),
                                          label=c(rep("1-probit",length(frac)),rep("root l2",length(frac)),rep("uniform",length(frac))),
                                          frac=c(rep(frac,3)))


p1_meandiff_heartattack_data <- data.frame(mean=c(heartattack_p_1_frac_mean_diff_norm_average,
                                                heartattack_p_1_lpcoreset_frac_mean_diff_norm_average,
                                                heartattack_p_1_uniform_frac_mean_diff_norm_average,
                                                heartattack_p_1_one_shot_mean_diff_norm_average),
                                         median=c(heartattack_p_1_frac_mean_diff_norm_median,
                                                  heartattack_p_1_lpcoreset_frac_mean_diff_norm_median,
                                                  heartattack_p_1_uniform_frac_mean_diff_norm_median,
                                                  heartattack_p_1_one_shot_mean_diff_norm_median),
                                         max=c(heartattack_p_1_frac_mean_diff_norm_max,
                                               heartattack_p_1_lpcoreset_frac_mean_diff_norm_max,
                                               heartattack_p_1_uniform_frac_mean_diff_norm_max,
                                               heartattack_p_1_one_shot_mean_diff_norm_max),
                                         min=c(heartattack_p_1_frac_mean_diff_norm_min,
                                               heartattack_p_1_lpcoreset_frac_mean_diff_norm_min,
                                               heartattack_p_1_uniform_frac_mean_diff_norm_min,
                                               heartattack_p_1_one_shot_mean_diff_norm_min),
                                         label=c(rep("1-probit",length(frac)),rep("root l2",length(frac)),rep("uniform",length(frac)),
                                                 rep("one-shot",length(frac))),
                                         frac=c(rep(frac,4)))



jpeg("heartdisease_p1_mean.jpeg", units="in", width=8, height=5, res=300)
ggplot(p1_meandiff_heartattack_data, aes(frac, median, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior mean", colour="method", fill="method") +
  theme_bw()+
  ggtitle("Heart disease data (p=1)")
dev.off()

###############post cov diff##########
######poserior covariance difference
heartattack_p1_pos_cov_diff_norm_median <- c()
heartattack_p1_pos_cov_diff_norm_max <- c()
heartattack_p1_pos_cov_diff_norm_min <- c()
heartattack_p1_pos_cov_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    b <- norm(cov(get(paste("heartattack_p_1_frac",frac[i],"chain_",j,sep = ""))$chain)-
                cov(heartattack_1$chain),type = "2")
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN = 2, mean)
  heartattack_p1_pos_cov_diff_norm_median <- c(heartattack_p1_pos_cov_diff_norm_median,median_result)
  heartattack_p1_pos_cov_diff_norm_max <- c(heartattack_p1_pos_cov_diff_norm_max,maxresult)
  heartattack_p1_pos_cov_diff_norm_min <- c(heartattack_p1_pos_cov_diff_norm_min,minresult)
  heartattack_p1_pos_cov_diff_norm_average <- c(heartattack_p1_pos_cov_diff_norm_average,averageresult)
  
  
}
####lpcoreset######
heartattack_p1_lpcoreset_pos_cov_diff_norm_median <- c()
heartattack_p1_lpcoreset_pos_cov_diff_norm_max <- c()
heartattack_p1_lpcoreset_pos_cov_diff_norm_min <- c()
heartattack_p1_lpcoreset_pos_cov_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    b <- norm(cov(get(paste("heartattack_p_1_lpcoreset_frac",frac[i],"chain_",j,sep = ""))$chain)-
                cov(heartattack_1$chain),type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN = 2, mean)
  heartattack_p1_lpcoreset_pos_cov_diff_norm_median <- c(heartattack_p1_lpcoreset_pos_cov_diff_norm_median,median_result)
  heartattack_p1_lpcoreset_pos_cov_diff_norm_max <- c(heartattack_p1_lpcoreset_pos_cov_diff_norm_max,maxresult)
  heartattack_p1_lpcoreset_pos_cov_diff_norm_min <- c(heartattack_p1_lpcoreset_pos_cov_diff_norm_min,minresult)
  heartattack_p1_lpcoreset_pos_cov_diff_norm_average <- c(heartattack_p1_lpcoreset_pos_cov_diff_norm_average,averageresult)
  
}

heartattack_p_1_one_shot_cov_diff_norm_median <- c()
heartattack_p_1_one_shot_cov_diff_norm_max <- c()
heartattack_p_1_one_shot_cov_diff_norm_min <- c()
heartattack_p_1_one_shot_cov_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    b <- norm(cov(get(paste("heartattack_p_1_one_shot_frac",frac[i],"chain_",j,sep = ""))$chain)-
                cov(heartattack_1$chain),type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  heartattack_p_1_one_shot_cov_diff_norm_median <- c(heartattack_p_1_one_shot_cov_diff_norm_median,median_result)
  heartattack_p_1_one_shot_cov_diff_norm_max <- c(heartattack_p_1_one_shot_cov_diff_norm_max,maxresult)
  heartattack_p_1_one_shot_cov_diff_norm_min <- c(heartattack_p_1_one_shot_cov_diff_norm_min,minresult)
  heartattack_p_1_one_shot_cov_diff_norm_average <- c(heartattack_p_1_one_shot_cov_diff_norm_average,meanresult)
}


heartattack_p1_uniform_pos_cov_diff_norm_median <- c()
heartattack_p1_uniform_pos_cov_diff_norm_max <- c()
heartattack_p1_uniform_pos_cov_diff_norm_min <- c()
heartattack_p1_uniform_pos_cov_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:10) {
    b <- norm(cov(get(paste("heartattack_p_1_uniform_frac",frac[i],"chain_",j,sep = ""))$chain)-
                cov(heartattack_1$chain),type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN = 2, mean)
  heartattack_p1_uniform_pos_cov_diff_norm_median <- c(heartattack_p1_uniform_pos_cov_diff_norm_median,median_result)
  heartattack_p1_uniform_pos_cov_diff_norm_max <- c(heartattack_p1_uniform_pos_cov_diff_norm_max,maxresult)
  heartattack_p1_uniform_pos_cov_diff_norm_min <- c(heartattack_p1_uniform_pos_cov_diff_norm_min,minresult)
  heartattack_p1_uniform_pos_cov_diff_norm_average <- c(heartattack_p1_uniform_pos_cov_diff_norm_average,averageresult)
}

p1_covdiff_heartattack_data <- data.frame(mean=c(log(heartattack_p1_pos_cov_diff_norm_average),
                                             log(heartattack_p1_lpcoreset_pos_cov_diff_norm_average),
                                             log(heartattack_p1_uniform_pos_cov_diff_norm_average),
                                             log(heartattack_p1_pos_cov_diff_norm_average)),
                                      median=c(log(heartattack_p1_pos_cov_diff_norm_median),
                                               log(heartattack_p1_lpcoreset_pos_cov_diff_norm_median),
                                               log(heartattack_p1_uniform_pos_cov_diff_norm_median),
                                               log(heartattack_p1_pos_cov_diff_norm_median)),
                                      max=c(log(heartattack_p1_pos_cov_diff_norm_max),
                                            log(heartattack_p1_lpcoreset_pos_cov_diff_norm_max),
                                            log(heartattack_p1_uniform_pos_cov_diff_norm_max),
                                            log(heartattack_p1_pos_cov_diff_norm_max)),
                                      min=c(log(heartattack_p1_pos_cov_diff_norm_min),
                                            log(heartattack_p1_lpcoreset_pos_cov_diff_norm_min),
                                            log(heartattack_p1_uniform_pos_cov_diff_norm_min),
                                            log(heartattack_p1_pos_cov_diff_norm_min)),
                                      label=c(rep("1-probit",length(frac)),rep("root l2",length(frac)),rep("uniform",length(frac)),
                                              rep("one-shot",length(frac))),
                                      frac=c(rep(frac,4)))


p1_covdiff_heartattack_data
ggplot(p1_covdiff_heartattack_data, aes(frac, median, fill=label, colour=label)) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior covariance", colour="method", fill="method") +
  theme_bw()
#####MMD#######
heartattack_p_1_MMD_median <- c()
heartattack_p_1_MMD_max <- c()
heartattack_p_1_MMD_min <- c()
heartattack_p_1_MMD_average <- c()

for (i in 1:length(frac)){
  print(paste("1-probit_MMD_Frac:",i,sep = ""))
  medianchain <- c()
  for (j in 1:10) {
    b <- mmd(sample_x=tail(get(paste("heartattack_p_1_frac",frac[i],"chain_",j,sep = ""))$chain,100),
             sample_y = tail(cover_type_1$chain,100))
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN=2, mean)
  heartattack_p_1_MMD_median <- c(heartattack_p_1_MMD_median,median_result)
  heartattack_p_1_MMD_max <- c(heartattack_p_1_MMD_max,maxresult)
  heartattack_p_1_MMD_min <- c(heartattack_p_1_MMD_min,minresult)
  heartattack_p_1_MMD_average <- c(heartattack_p_1_MMD_average,averageresult)
}

heartattack_p_1_lpcoreset_MMD_median <- c()
heartattack_p_1_lpcoreset_MMD_max <- c()
heartattack_p_1_lpcoreset_MMD_min <- c()
heartattack_p_1_lpcoreset_MMD_average <- c()
for (i in 1:length(frac)){
  print(paste("rootl2_MMD_Frac:",i,sep = ""))
  medianchain <- c()
  for (j in 1:10) {
    b <- mmd(sample_x=tail(get(paste("heartattack_p_1_lpcoreset_frac",frac[i],"chain_",j,sep = ""))$chain,100),
             sample_y = tail(cover_type_1$chain,100))
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN=2, mean)
  heartattack_p_1_lpcoreset_MMD_median <- c(heartattack_p_1_lpcoreset_MMD_median,median_result)
  heartattack_p_1_lpcoreset_MMD_max <- c(heartattack_p_1_lpcoreset_MMD_max,maxresult)
  heartattack_p_1_lpcoreset_MMD_min <- c(heartattack_p_1_lpcoreset_MMD_min,minresult)
  heartattack_p_1_lpcoreset_MMD_average <- c(heartattack_p_1_lpcoreset_MMD_average,averageresult)
}

heartattack_p_1_uniform_MMD_median <- c()
heartattack_p_1_uniform_MMD_max <- c()
heartattack_p_1_uniform_MMD_min <- c()
heartattack_p_1_uniform_MMD_average <- c()
for (i in 1:length(frac)){
  print(paste("uniform_MMD_Frac:",i,sep = ""))
  medianchain <- c()
  for (j in 1:5) {
    b <- mmd(sample_x=tail(get(paste("heartattack_p_1_uniform_frac",frac[i],"chain_",j,sep = ""))$chain,100),
             sample_y = tail(cover_type_1$chain,100))
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN=2, mean)
  heartattack_p_1_uniform_MMD_median <- c(heartattack_p_1_uniform_MMD_median,median_result)
  heartattack_p_1_uniform_MMD_max <- c(heartattack_p_1_uniform_MMD_max,maxresult)
  heartattack_p_1_uniform_MMD_min <- c(heartattack_p_1_uniform_MMD_min,minresult)
  heartattack_p_1_uniform_MMD_average <- c(heartattack_p_1_uniform_MMD_average,averageresult)
}

p1_mmd_heartattack_data <- data.frame(mean=c(heartattack_p_1_MMD_average,
                                           heartattack_p_1_lpcoreset_MMD_average,
                                           heartattack_p_1_uniform_MMD_average),
                                    median=c(heartattack_p_1_MMD_median,
                                             heartattack_p_1_lpcoreset_MMD_median,
                                             heartattack_p_1_uniform_MMD_median),
                                    max=c(heartattack_p_1_MMD_max,
                                          heartattack_p_1_lpcoreset_MMD_max,
                                          heartattack_p_1_uniform_MMD_max),
                                    min=c(heartattack_p_1_MMD_min,
                                          heartattack_p_1_lpcoreset_MMD_min,
                                          heartattack_p_1_uniform_MMD_min),
                                    label=c(rep("1-probit",length(frac)),rep("root l2",length(frac)),rep("uniform",length(frac))),
                                    frac=c(rep(frac,3)))


p1_mmd_heartattack_data <- data.frame(mean=c(log(heartattack_p_1_MMD_average),
                                             log(heartattack_p_1_lpcoreset_MMD_average),
                                             log(heartattack_p_1_uniform_MMD_average)),
                                      median=c(log(heartattack_p_1_MMD_median),
                                               log(heartattack_p_1_lpcoreset_MMD_median),
                                               log(heartattack_p_1_uniform_MMD_median)),
                                      max=c(log(heartattack_p_1_MMD_max),
                                            log(heartattack_p_1_lpcoreset_MMD_max),
                                            log(heartattack_p_1_uniform_MMD_max)),
                                      min=c(log(heartattack_p_1_MMD_min),
                                            log(heartattack_p_1_lpcoreset_MMD_min),
                                            log(heartattack_p_1_uniform_MMD_min)),
                                      label=c(rep("1-probit",length(frac)),rep("root l2",length(frac)),rep("uniform",length(frac))),
                                      frac=c(rep(frac,3)))
p1_mmd_heartattack_data
ggplot(p1_meandiff_heartattack_data, aes(frac, median, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior mean", colour="method", fill="method") +
  theme_bw()
ggplot(p1_covdiff_heartattack_data, aes(frac, median, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior covariance", colour="method", fill="method") +
  theme_bw()
ggplot(p1_mmd_heartattack_data, aes(frac, median, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="MMD", colour="method", fill="method") +
  theme_bw()
######################

########fixed p=2#########
#####algorithm 1#####
for (i in 1:length(frac)) {
  
  print(paste("Model_Frac:",i,"chain:",j,sep = ""))
  sample_index<- Compute_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=2)
  Model <- multi_chain(N_sim=N_sim,burn_in=burn_in,
                       X=X[sample_index,], y=y[sample_index],initial_theda=probit_coefficients,
                       true_theta=probit_coefficients,Lp=3,M_iter=100,range=c(0.1,5),step=0.02,times = 1)
  assign(paste("heartattack_p2_nofix_frac",frac[i],sep = ""), Model)
  
}



for (i in 1:length(frac)) {
  
  print(paste("Model_Frac:",i,"chain:",j,sep = ""))
  sample_index<- Compute_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=1)
  Model <- multi_chain(N_sim=N_sim,burn_in=burn_in,
                       X=X[sample_index,], y=y[sample_index],initial_theda=probit_coefficients,
                       true_theta=probit_coefficients,Lp=3,M_iter=100,range=c(0.1,5),step=0.02,times = 1)
  assign(paste("heartattack_p1_nofix_frac",frac[i],sep = ""), Model)
  
}



#####LP sampling#########
for (i in 1:length(frac)) {
  print(paste("Model_Frac:",i,"chain:",j,sep = ""))
  sample_index<- Lp_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=2)
  Model <- multi_chain(N_sim=N_sim,burn_in=burn_in,
                       X=X[sample_index,], y=y[sample_index],initial_theda=probit_coefficients,
                       true_theta=probit_coefficients,Lp=3,M_iter=100,range=c(0.1,5),step=0.02,times = 3)
  assign(paste("heartattack_p_lpcoreset_nofix_frac",frac[i],sep = ""), Model)
  
}


#########uniform sampling
for (i in 1:length(frac)) {
  print(paste("Model_Frac:",i,"chain:",j,sep = ""))
  sample_index<- Uniform_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=2)
  Model <- multi_chain(N_sim=N_sim,burn_in=burn_in,
                       X=X[sample_index,], y=y[sample_index],initial_theda=probit_coefficients,
                       true_theta=probit_coefficients,Lp=3,M_iter=100,range=c(0.1,5),step=0.02,times = 3)
  assign(paste("heartattack_p_uniform_nofix_frac",frac[i],sep = ""), Model)
}

for (i in 1:length(frac)) {
  
  assign(paste("heartattack_p_uniform_nofix_frac",frac[i],sep = ""), get(paste("hearattack_p_uniform_nofix_frac",frac[i],sep = "")))
  
}

heartattack_p_summary1 <- c()
for (i in 1:length(frac)){
  
  b <- get(paste("heartattack_p2_nofix_frac",frac[i],sep = ""))$Lp.mean
  
  heartattack_p_summary1<- rbind(heartattack_p_summary1,b)
  
}
plot(frac,heartattack_p_summary1,type = "l")

heartattack_p_summary2 <- c()
for (i in 1:length(frac)){
  
  b <- get(paste("heartattack_p_lpcoreset_nofix_frac",frac[i],sep = ""))$Lp.mean
  
  heartattack_p_summary2<- rbind(heartattack_p_summary2,b)
  
}
plot(frac,heartattack_p_summary2,type = "l")

heartattack_p_summary3 <- c()

for (i in 1:length(frac)){
  
  b <- get(paste("heartattack_p_uniform_nofix_frac",frac[i],sep = ""))$Lp.mean
  
  heartattack_p_summary3<- rbind(heartattack_p_summary3,b)
  
}
plot(frac,heartattack_p_summary2,type = "l")
lines(frac,heartattack_p_summary3,col=2)

plot(frac,heartattack_p_summary1,type = "l",ylim = c(0,3),xlab = "sample fraction", ylab = "estimated P")
lines(frac,heartattack_p_summary2,col=2)
lines(frac,heartattack_p_summary3,col=3)
legend("topright",legend=c("2-probit","root l2","uniform"), col=c(1,2,3), 
       bty = 'n', lwd = 2, inset = c(0, 0), lty = 1, cex = 0.73)

