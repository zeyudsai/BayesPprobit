rm(list = ls())
dev.off()
options(scipen=100, digits=8)


setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("Functions.R")


d<- fread("Real_world_data/covtype.data.gz")

Wilderness_Area<-c()
Soil_Type <- c()
for (i in 1:4) {
  Wilderness_Area[i]<-paste("Wilderness_Area.",i,sep = "")
  
}
for (i in 1:40) {
  Soil_Type[i]<-paste("Soil_Type.",i,sep = "")
  
}
names(d)<-c("Elevation","Aspect","Slope","Horizontal_Distance_To_Hydrology",
            "Vertical_Distance_To_Hydrology","Horizontal_Distance_To_Roadways",
            "Hillshade_9am","Hillshade_Noon","Hillshade_3pm","Horizontal_Distance_To_Fire_Points",
            Wilderness_Area,Soil_Type,"Cover_Type")
Wilderness_Area
as.data.frame(table(d$Wilderness_Area.4))
as.data.frame(table(d$Soil_Type.40))
hist(d$Cover_Type)
d$Cover_Type[d$Cover_Type!=2] <- 0
d$Cover_Type[d$Cover_Type==2] <- 1
d <- d[,c(-14,-54)]
D <- ncol(d)
N <- nrow(d)
X <- d[,1:52]
ncol(X)
head(X)
onevector <- rep(1,N)

X <- as.matrix(cbind(onevector,X))

y <- d$Cover_Type

mle_result <- glm(Cover_Type~.,data=d,family = binomial(link = "logit"))
mle_probit <- glm(Cover_Type~.,data=d,family = binomial(link = probit))
MLE_coefficients <- mle_result$coefficients
probit_coefficients <- mle_probit$coefficients
MLE_coefficients[is.na(MLE_coefficients)] <- 0
probit_coefficients[is.na(probit_coefficients)] <- 0

N_sim <- 1000
burn_in <- N_sim*0.5
cover_type_2 <- Lp_gibbssampler(N_sim=N_sim,
                                burn_in=burn_in, X=X, y=y,true_theta = probit_coefficients,Lp=2,initial_theda = probit_coefficients)
MLE_coefficients
cover_type_2$post.mean
as.vector(MLE_coefficients)
as.vector(probit_coefficients)
norm(na.omit(cover_type_2$post.mean-MLE_coefficients),type = "2")
norm(na.omit(cover_type_2$post.mean-probit_coefficients),type = "2")

cover_type_1 <- Lp_gibbssampler(N_sim=N_sim,
                burn_in=burn_in, X=X, y=y,true_theta = probit_coefficients,Lp=1,initial_theda = probit_coefficients)
cover_type_1$post.mean

as.vector(MLE_coefficients)
as.vector(probit_coefficients)
norm(na.omit(cover_type_1$post.mean-MLE_coefficients),type = "2")
norm(na.omit(cover_type_1$post.mean-probit_coefficients),type = "2")

sample3 <- multi_chain(N_sim=N_sim,burn_in=burn_in,
                       X=X, y=y,initial_theda=probit_coefficients,
                       true_theta=probit_coefficients,Lp=3,M_iter=300,range=c(1,5),step=0.05,times = 1)

sample4 <- multi_chain(N_sim=N_sim,burn_in=burn_in,
                       X=X, y=y,initial_theda=probit_coefficients,
                       true_theta=probit_coefficients,Lp=3,M_iter=300,range=c(0,5),step=0.05,times = 1)



N_sim
sample3$Lp.mean.chain
sample3$Beta.mean
cover_type_1


plot(100:N_sim, sample3$Lp.mean.chain[100:N_sim],ylim=c(0.999,1.01),type = "l",xlab = "MCMC iteration",ylab = "estimated P",
     main = "The estimated P of the Covertype data")

Msubtitle <- "P is restricted between (1,5)"
mtext(side = 3, line = 0.4, Msubtitle)
sample4

plot(100:N_sim, sample4$Lp.mean.chain[100:2000],type = "l",xlab = "MCMC iteration",ylab = "estimated P",
     main = "The estimated P of the Covertype data")
Msubtitle <- "Estimated P 0.67"
mtext(side = 3, line = 0.4, Msubtitle)

head(as.data.frame(sample4$Lp.mean.chain))

setwd("/Users/dingzeyu/Documents/文稿 - MacBook Pro/Research Project/Generalized Normal Distribution using MCMC/plots")
jpeg("cover_type_p1.jpeg", units="in", width=8, height=5, res=300)
ggplot(data=data.frame(Lp.mean.chain=sample3$Lp.mean.chain[seq(250,1000,3)],x=seq(250,1000,3)), aes(x=x, y=Lp.mean.chain)) +
  geom_line()+
  theme_bw(base_size = 12) +  
  annotate("text", x=750, y=1.005,size=6, label= TeX("$\\hat{p} = 1.001", output='character'),parse=TRUE)+
  labs(title="Estimation of p for covertype data", 
       # subtitle="Estimated p=0.673", 
       x="MCMC iteration",
       y="p")+ylim(c(1,1.01))
# insert ggplot code
dev.off()

jpeg("cover_type_p2.jpeg", units="in", width=8, height=5, res=300)

ggplot(data=data.frame(Lp.mean.chain=sample4$Lp.mean.chain[seq(250,1000,5)],x=seq(250,1000,5)), aes(x=x, y=Lp.mean.chain)) +
  geom_line()+
  theme_bw(base_size = 12) +  
  annotate("text", x=750, y=1.008,size=6, label= TeX("$\\hat{p} = 0.673", output='character'),parse=TRUE)+
  labs(title="Estimation of p for covertype data", 
       # subtitle="Estimated p=0.673", 
       x="MCMC iteration",
       y="p")+ylim(c(0.5,1.5))
dev.off()

# insert ggplot code


jpeg("cover_type_p2.jpeg", units="in", width=8, height=5, res=300)
ggplot(data=data.frame(Lp.mean.chain=sample4$Lp.mean.chain[250:N_sim],x=250:length(sample3$Lp.mean.chain)), aes(x=x, y=Lp.mean.chain)) +
  geom_line()+
  theme_bw(base_size = 12) +  
  annotate("text", x=800, y=0.85,size=6, label= TeX("$\\hat{p} = 0.673", output='character'),parse=TRUE)+
  labs(title="Estimation of p for covertype data", 
       # subtitle="Estimated p=0.673", 
       x="MCMC iteration",
       y="p")+ylim(0.5,1.5)
dev.off()



sample3ci<-apply(sample3$Beta_chain[[1]], MARGIN =2, FUN=ci)
sample3_ci_low <- c()
sample3_ci_high <- c()
for (i in 1:D) {
  sample3_ci_low[i] <- sample3ci[[i]]$CI_low
  sample3_ci_high[i] <- sample3ci[[i]]$CI_high

}
covertype_1_ci<-apply(cover_type_1$chain[800:1000,], MARGIN =2, FUN=ci)
covertype_1_ci
covertype_1_ci_low <- c()
covertype_1_ci_high <- c()
for (i in 1:D) {
  covertype_1_ci_low[i] <- covertype_1_ci[[i]]$CI_low
  covertype_1_ci_high[i] <- covertype_1_ci[[i]]$CI_high
  
}
cover_type_2$chain
covertype_2_ci<-apply(cover_type_2$chain[800:1000,], MARGIN =2, FUN=ci)
covertype_2_ci
covertype_2_ci_low <- c()
covertype_2_ci_high <- c()
for (i in 1:D) {
  covertype_2_ci_low[i] <- covertype_2_ci[[i]]$CI_low
  covertype_2_ci_high[i] <- covertype_2_ci[[i]]$CI_high
  
}

result <-summary(mle_probit)
probit_confidence_low <- result$coefficients[,1]-result$coefficients[,2]*qt(0.95,df=N-2)
probit_confidence_high <- result$coefficients[,1]+result$coefficients[,2]*qt(0.95,df=N-2)
probit_confidence_high <- probit_confidence_high[c(-15,-19,-21,-28,-29,-51)]
probit_confidence_low <- probit_confidence_low[c(-15,-19,-21,-28,-29,-51)]
covertype_1_ci_low <- covertype_1_ci_low[c(-15,-19,-21,-28,-29,-51)]
covertype_1_ci_high <- covertype_1_ci_high[c(-15,-19,-21,-28,-29,-51)]
covertype_2_ci_low <- covertype_2_ci_low[c(-15,-19,-21,-28,-29,-51)]
covertype_2_ci_high <- covertype_2_ci_high[c(-15,-19,-21,-28,-29,-51)]
result <-summary(mle_result)
logit_confidence_low <- result$coefficients[,1]-result$coefficients[,2]*qt(0.95,df=N-2)
logit_confidence_high <- result$coefficients[,1]+result$coefficients[,2]*qt(0.95,df=N-2)
logit_confidence_high <- logit_confidence_high[c(-15,-19,-21,-28,-29,-51)]
logit_confidence_low <- logit_confidence_low[c(-15,-19,-21,-28,-29,-51)]

ci_plot1 <- data.frame("variable"=c(1:(D-6)),
  "sample3_ci_low"=c(covertype_1_ci_low,covertype_2_ci_low,
                     probit_confidence_low,logit_confidence_low),
                       "sample3_ci_high"=c(covertype_1_ci_high,covertype_2_ci_high,
                                           probit_confidence_high,logit_confidence_high),
                       "post_mean"=c(cover_type_1$post.mean[c(-15,-19,-21,-28,-29,-51)],
                                     cover_type_2$post.mean[c(-15,-19,-21,-28,-29,-51)],
                                     mle_probit$coefficients[c(-15,-19,-21,-28,-29,-51)],
                                     mle_result$coefficients[c(-15,-19,-21,-28,-29,-51)]),
                       "Method"=c(rep("Bayesian(P=1) Credible Interval",(D-6)),
                                  rep("Bayesian(P=2) Credible Interval",(D-6)),
                                  rep("Probit Confidence Interval",(D-6)),
                                  rep("Logit Confidence Interval",(D-6))))
ci_plot1 <- subset(ci_plot1,ci_plot1$variable>12)
ci_plot1

jpeg("cover_type_ci1.jpeg", units="in", width=8, height=5, res=300)
ggplot(ci_plot1, aes(variable, post_mean, colour=Method)) +
  geom_point(position=position_dodge(1)) + 
  geom_errorbar(aes(ymin=sample3_ci_low, ymax=sample3_ci_high), width=0.5, position=position_dodge(1)) + 
  theme_bw()+
  labs(x=TeX("$\\hat{\\beta}\\ index"), y=TeX("$\\hat{\\beta}"), colour="Method", fill="Method")+ 
  ggtitle("Covertype Data") 
dev.off()


plot(1:D,probit_coefficients,type = "l",ylim = c(-10,10),main = "comparison of coefficients",xlab = "Beta-index",
     ylab = "parameters")
Msubtitle <- "covertype data"
mtext(side = 3, line = 0.4, Msubtitle)
legend("topright",legend=c("MLE Probit","Bayesian-1-probit","Bayesian-2-probit"), col=c(1,2,3), 
       bty = 'n', lwd = 2, inset = c(0, 0), lty = 1, cex = 0.73)
lines(cover_type_1$post.mean,col=2)
lines(cover_type_2$post.mean,col=3)

plot(1:D,probit_coefficients,type = "l",ylim = c(-20,15),main = "comparison of coefficients",xlab = "Beta-index",
     ylab = "parameters")
Msubtitle <- "covertype data"
mtext(side = 3, line = 0.4, Msubtitle)
legend("topright",legend=c("MLE Probit","Bayesian-1(estimated)-probit","Bayesian-0.6(estimated)-probit"), col=c(1,2,3), 
       bty = 'n', lwd = 2, inset = c(0, 0), lty = 1, cex = 0.73)
lines(sample3$Beta.mean,col=2)
lines(sample4$Beta.mean,col=3)

plot(1:D,MLE_coefficients,type = "l",ylim = c(-20,20),main = "comparison of coefficients",xlab = "Beta-index",
     ylab = "parameters")
Msubtitle <- "covertype data"
mtext(side = 3, line = 0.4, Msubtitle)
legend("topright",legend=c("MLE Logit","Bayesian-1-probit","Bayesian-2-probit"), col=c(1,2,3), 
       bty = 'n', lwd = 2, inset = c(0, 0), lty = 1, cex = 0.73)
lines(cover_type_1$post.mean,col=2)
lines(cover_type_2$post.mean,col=3)

plot(1:D,MLE_coefficients,type = "l",ylim = c(-20,20),main = "comparison of coefficients",xlab = "Beta-index",
     ylab = "parameters")
Msubtitle <- "covertype data"
mtext(side = 3, line = 0.4, Msubtitle)
legend("topright",legend=c("MLE Logit","Bayesian-1(estimated)-probit","Bayesian-0.6(estimated)-probit"), col=c(1,2,3), 
       bty = 'n', lwd = 2, inset = c(0, 0), lty = 1, cex = 0.73)
lines(sample3$Beta.mean,col=2)
lines(sample4$Beta.mean,col=3)

sketch_size <- D^2
frac <-seq(0.00017211348,0.05,0.00086056742)


########fixed p=2#########
#####algorithm 1#####
for (i in 1:length(frac)) {
  for (j in 1:5) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Compute_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=2)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y[sample_index],
                             initial_theda=probit_coefficients,
                             true_theta=probit_coefficients,Lp=2)
    assign(paste("covertype_p_2_frac",frac[i],"chain_",j,sep = ""), Model)
  }
}
#####LP sampling#########
# for (i in 1:length(frac)) {
#   for (j in 1:10) {
#     print(paste("Model_Frac:",i,"chain:",j,sep = ""))
#     sample_index<- Lp_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=2)
#     Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
#                              X=X[sample_index,], y=y[sample_index],
#                              initial_theda=probit_coefficients,
#                              true_theta=probit_coefficients,Lp=2)
#     assign(paste("covertype_p_2_lpcoreset_frac",frac[i],"chain_",j,sep = ""), Model)
#   }
# }
#########uniform sampling
for (i in 1:length(frac)) {
  for (j in 1:5) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Uniform_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=2)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y[sample_index],
                             initial_theda=probit_coefficients,
                             true_theta=probit_coefficients,Lp=2)
    assign(paste("covertype_p_2_uniform_frac",frac[i],"chain_",j,sep = ""), Model)
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
    assign(paste("covertype_p_2_one_shot",frac[i],"chain_",j,sep = ""), Model)
  }
}


#############post mean difference######
covertype_p_2_frac_mean_diff_norm_median <- c()
covertype_p_2_frac_mean_diff_norm_max <- c()
covertype_p_2_frac_mean_diff_norm_min <- c()
covertype_p_2_frac_mean_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    b <- norm(get(paste("covertype_p_2_frac",frac[i],"chain_",j,sep = ""))$post.mean-cover_type_2$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  covertype_p_2_frac_mean_diff_norm_median <- c(covertype_p_2_frac_mean_diff_norm_median,median_result)
  covertype_p_2_frac_mean_diff_norm_max <- c(covertype_p_2_frac_mean_diff_norm_max,maxresult)
  covertype_p_2_frac_mean_diff_norm_min <- c(covertype_p_2_frac_mean_diff_norm_min,minresult)
  covertype_p_2_frac_mean_diff_norm_average <- c(covertype_p_2_frac_mean_diff_norm_average,meanresult)
}

# covertype_p_2_lpcoreset_frac_mean_diff_norm_median <- c()
# covertype_p_2_lpcoreset_frac_mean_diff_norm_max <- c()
# covertype_p_2_lpcoreset_frac_mean_diff_norm_min <- c()
# covertype_p_2_lpcoreset_frac_mean_diff_norm_average <- c()
# for (i in 1:length(frac)){
#   medianchain <- c()
#   for (j in 1:10) {
#     
#     b <- norm(get(paste("covertype_p_2_lpcoreset_frac",frac[i],"chain_",j,sep = ""))$post.mean-cover_type_2$post.mean,type = "2")
#     medianchain<- rbind(medianchain,b)
#   }
#   median_result <- apply(medianchain, MARGIN = 2, median)
#   maxresult <- apply(medianchain, MARGIN = 2, max)
#   minresult <- apply(medianchain, MARGIN = 2, min)
#   meanresult <- apply(medianchain, MARGIN = 2, mean)
#   covertype_p_2_lpcoreset_frac_mean_diff_norm_median <- c(covertype_p_2_lpcoreset_frac_mean_diff_norm_median,median_result)
#   covertype_p_2_lpcoreset_frac_mean_diff_norm_max <- c(covertype_p_2_lpcoreset_frac_mean_diff_norm_max,maxresult)
#   covertype_p_2_lpcoreset_frac_mean_diff_norm_min <- c(covertype_p_2_lpcoreset_frac_mean_diff_norm_min,minresult)
#   covertype_p_2_lpcoreset_frac_mean_diff_norm_average <- c(covertype_p_2_lpcoreset_frac_mean_diff_norm_average,meanresult)
#   
# }

covertype_p_2_one_shot_mean_diff_norm_median <- c()
covertype_p_2_one_shot_mean_diff_norm_max <- c()
covertype_p_2_one_shot_mean_diff_norm_min <- c()
covertype_p_2_one_shot_mean_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    b <- norm(get(paste("covertype_p_2_one_shot",frac[i],"chain_",j,sep = ""))$post.mean-cover_type_2$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  covertype_p_2_one_shot_mean_diff_norm_median <- c(covertype_p_2_one_shot_mean_diff_norm_median,median_result)
  covertype_p_2_one_shot_mean_diff_norm_max <- c(covertype_p_2_one_shot_mean_diff_norm_max,maxresult)
  covertype_p_2_one_shot_mean_diff_norm_min <- c(covertype_p_2_one_shot_mean_diff_norm_min,minresult)
  covertype_p_2_one_shot_mean_diff_norm_average <- c(covertype_p_2_one_shot_mean_diff_norm_average,meanresult)
}

covertype_p_2_uniform_frac_mean_diff_norm_median <- c()
covertype_p_2_uniform_frac_mean_diff_norm_max <- c()
covertype_p_2_uniform_frac_mean_diff_norm_min <- c()
covertype_p_2_uniform_frac_mean_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    b <- norm(get(paste("covertype_p_2_uniform_frac",frac[i],"chain_",j,sep = ""))$post.mean-cover_type_2$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  covertype_p_2_uniform_frac_mean_diff_norm_median <- c(covertype_p_2_uniform_frac_mean_diff_norm_median,median_result)
  covertype_p_2_uniform_frac_mean_diff_norm_max <- c(covertype_p_2_uniform_frac_mean_diff_norm_max,maxresult)
  covertype_p_2_uniform_frac_mean_diff_norm_min <- c(covertype_p_2_uniform_frac_mean_diff_norm_min,minresult)
  covertype_p_2_uniform_frac_mean_diff_norm_average <- c(covertype_p_2_uniform_frac_mean_diff_norm_average,meanresult)
  
}

p2_meandiff_covertype_data <- data.frame(mean=c(covertype_p_2_frac_mean_diff_norm_average,
                                                covertype_p_2_one_shot_mean_diff_norm_average,
                                     covertype_p_2_uniform_frac_mean_diff_norm_average),
                            median=c(covertype_p_2_frac_mean_diff_norm_median,
                                     covertype_p_2_one_shot_mean_diff_norm_median,
                                   covertype_p_2_uniform_frac_mean_diff_norm_median),
                            max=c(covertype_p_2_frac_mean_diff_norm_max,
                                  covertype_p_2_one_shot_mean_diff_norm_max,
                                     covertype_p_2_uniform_frac_mean_diff_norm_max),
                            min=c(covertype_p_2_frac_mean_diff_norm_min,
                                  covertype_p_2_one_shot_mean_diff_norm_min,
                                  covertype_p_2_uniform_frac_mean_diff_norm_min),
                            label=c(rep("2-probit",length(frac)),rep("one-shot",length(frac)),rep("uniform",length(frac))),
                            frac=c(rep(frac,3)))

jpeg("cover_type_p2_mean.jpeg", units="in", width=8, height=5, res=300)

ggplot(p2_meandiff_covertype_data, aes(frac, mean, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior mean", colour="method", fill="method") +
  theme_bw()+
  ggtitle("Covertype data (p=2)")


dev.off()


######poserior covariance difference
covertype_p2_pos_cov_diff_norm_median <- c()
covertype_p2_pos_cov_diff_norm_max <- c()
covertype_p2_pos_cov_diff_norm_min <- c()
covertype_p2_pos_cov_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("covertype_p_2_frac",frac[i],"chain_",j,sep = ""))$chain)-
                cov(cover_type_2$chain),type = "2")
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN = 2, mean)
  covertype_p2_pos_cov_diff_norm_median <- c(covertype_p2_pos_cov_diff_norm_median,median_result)
  covertype_p2_pos_cov_diff_norm_max <- c(covertype_p2_pos_cov_diff_norm_max,maxresult)
  covertype_p2_pos_cov_diff_norm_min <- c(covertype_p2_pos_cov_diff_norm_min,minresult)
  covertype_p2_pos_cov_diff_norm_average <- c(covertype_p2_pos_cov_diff_norm_average,averageresult)
  
  
}

covertype_p2_lpcoreset_pos_cov_diff_norm_median <- c()
covertype_p2_lpcoreset_pos_cov_diff_norm_max <- c()
covertype_p2_lpcoreset_pos_cov_diff_norm_min <- c()
covertype_p2_lpcoreset_pos_cov_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("covertype_p_2_lpcoreset_frac",frac[i],"chain_",j,sep = ""))$chain)-
                cov(cover_type_2$chain),type = "2")
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN = 2, mean)
  covertype_p2_lpcoreset_pos_cov_diff_norm_median <- c(covertype_p2_lpcoreset_pos_cov_diff_norm_median,median_result)
  covertype_p2_lpcoreset_pos_cov_diff_norm_max <- c(covertype_p2_lpcoreset_pos_cov_diff_norm_max,maxresult)
  covertype_p2_lpcoreset_pos_cov_diff_norm_min <- c(covertype_p2_lpcoreset_pos_cov_diff_norm_min,minresult)
  covertype_p2_lpcoreset_pos_cov_diff_norm_average <- c(covertype_p2_lpcoreset_pos_cov_diff_norm_average,averageresult)
  
}


covertype_p_2_one_shot_cov_diff_norm_median <- c()
covertype_p_2_one_shot_cov_diff_norm_max <- c()
covertype_p_2_one_shot_cov_diff_norm_min <- c()
covertype_p_2_one_shot_cov_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("covertype_p_2_one_shot",frac[i],"chain_",j,sep = ""))$chain)-
                cov(cover_type_2$chain),type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  covertype_p_2_one_shot_cov_diff_norm_median <- c(covertype_p_2_one_shot_cov_diff_norm_median,median_result)
  covertype_p_2_one_shot_cov_diff_norm_max <- c(covertype_p_2_one_shot_cov_diff_norm_max,maxresult)
  covertype_p_2_one_shot_cov_diff_norm_min <- c(covertype_p_2_one_shot_cov_diff_norm_min,minresult)
  covertype_p_2_one_shot_cov_diff_norm_average <- c(covertype_p_2_one_shot_cov_diff_norm_average,meanresult)
}

covertype_p2_uniform_pos_cov_diff_norm_median <- c()
covertype_p2_uniform_pos_cov_diff_norm_max <- c()
covertype_p2_uniform_pos_cov_diff_norm_min <- c()
covertype_p2_uniform_pos_cov_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("covertype_p_2_uniform_frac",frac[i],"chain_",j,sep = ""))$chain)-
                cov(cover_type_2$chain),type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN = 2, mean)
  covertype_p2_uniform_pos_cov_diff_norm_median <- c(covertype_p2_uniform_pos_cov_diff_norm_median,median_result)
  covertype_p2_uniform_pos_cov_diff_norm_max <- c(covertype_p2_uniform_pos_cov_diff_norm_max,maxresult)
  covertype_p2_uniform_pos_cov_diff_norm_min <- c(covertype_p2_uniform_pos_cov_diff_norm_min,minresult)
  covertype_p2_uniform_pos_cov_diff_norm_average <- c(covertype_p2_uniform_pos_cov_diff_norm_average,averageresult)
}

p2_covdiff_covertype_data <- data.frame(mean=c(covertype_p2_pos_cov_diff_norm_median,
                                               covertype_p_2_one_shot_cov_diff_norm_average,
                                                covertype_p2_uniform_pos_cov_diff_norm_average),
                                         median=c(covertype_p2_pos_cov_diff_norm_median,
                                                  covertype_p_2_one_shot_cov_diff_norm_median,
                                                  covertype_p2_uniform_pos_cov_diff_norm_median),
                                         max=c(covertype_p2_pos_cov_diff_norm_max,
                                               covertype_p_2_one_shot_cov_diff_norm_max,
                                               covertype_p2_uniform_pos_cov_diff_norm_max),
                                         min=c(covertype_p2_pos_cov_diff_norm_min,
                                               covertype_p_2_one_shot_cov_diff_norm_min,
                                               covertype_p2_uniform_pos_cov_diff_norm_min),
                                         label=c(rep("2-probit",length(frac)),rep("one-shot",length(frac)),rep("uniform",length(frac))),
                                         frac=c(rep(frac,3)))
p2_covdiff_covertype_data

jpeg("cover_type_p2_cov.jpeg", units="in", width=8, height=5, res=300)

ggplot(p2_covdiff_covertype_data, aes(frac, mean, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior covariance", colour="method", fill="method") +
  theme_bw()+
  ggtitle("Covertype data (p=2)")
dev.off()

############lilelihood ratio###############


covertype_p2_2probit_llk_ratio_median <- c()
covertype_p2_2probit_llk_ratio_max <- c()
covertype_p2_2probit_llk_ratio_min <- c()
covertype_p2_2probit_llk_ratio_average <- c()




for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    b<- median(abs(llk(X, cover_type_2$post.mean,2,y)/
                     llk(X,get(paste("covertype_p_2_frac",frac[i],"chain_",j,sep = ""))$post.mean,2,y)))
    
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  covertype_p2_2probit_llk_ratio_median <- c(covertype_p2_2probit_llk_ratio_median,median_result)
  covertype_p2_2probit_llk_ratio_max <- c(covertype_p2_2probit_llk_ratio_max,maxresult)
  covertype_p2_2probit_llk_ratio_min <- c(covertype_p2_2probit_llk_ratio_min,minresult)
  covertype_p2_2probit_llk_ratio_average <- c(covertype_p2_2probit_llk_ratio_average,meanresult)
}

covertype_p2_uniform_llk_ratio_median <- c()
covertype_p2_uniform_llk_ratio_max <- c()
covertype_p2_uniform_llk_ratio_min <- c()
covertype_p2_uniform_llk_ratio_average <- c()



for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    

    b<- median(abs(llk(X, cover_type_2$post.mean,2,y)/
                     llk(X,get(paste("covertype_p_2_uniform_frac",frac[i],"chain_",j,sep = ""))$post.mean,2,y)))
    
  
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  covertype_p2_uniform_llk_ratio_median <- c(covertype_p2_uniform_llk_ratio_median,median_result)
  covertype_p2_uniform_llk_ratio_max <- c(covertype_p2_uniform_llk_ratio_max,maxresult)
  covertype_p2_uniform_llk_ratio_min <- c(covertype_p2_uniform_llk_ratio_min,minresult)
  covertype_p2_uniform_llk_ratio_average <- c(covertype_p2_uniform_llk_ratio_average,meanresult)
}

covertype_p2_one_shot_llk_ratio_median <- c()
covertype_p2_one_shot_llk_ratio_max <- c()
covertype_p2_one_shot_llk_ratio_min <- c()
covertype_p2_one_shot_llk_ratio_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
 
    b<- median(abs(llk(X, cover_type_2$post.mean,2,y)/
                     llk(X,get(paste("covertype_p_2_one_shot",frac[i],"chain_",j,sep = ""))$post.mean,2,y)))
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  covertype_p2_one_shot_llk_ratio_median <- c(covertype_p2_one_shot_llk_ratio_median,median_result)
  covertype_p2_one_shot_llk_ratio_max <- c(covertype_p2_one_shot_llk_ratio_max,maxresult)
  covertype_p2_one_shot_llk_ratio_min <- c(covertype_p2_one_shot_llk_ratio_min,minresult)
  covertype_p2_one_shot_llk_ratio_average <- c(covertype_p2_one_shot_llk_ratio_average,meanresult)
}

p2_llk_ratio_covertype_data <- data.frame(mean=c(covertype_p2_2probit_llk_ratio_average,
                                                 covertype_p2_one_shot_llk_ratio_average,
                                                 covertype_p2_uniform_llk_ratio_average),
                                        median=c(covertype_p2_2probit_llk_ratio_median,
                                                 covertype_p2_one_shot_llk_ratio_median,
                                                 covertype_p2_uniform_llk_ratio_median),
                                        max=c(covertype_p2_2probit_llk_ratio_max,
                                              covertype_p2_one_shot_llk_ratio_max,
                                              covertype_p2_uniform_llk_ratio_max),
                                        min=c(covertype_p2_2probit_llk_ratio_min,
                                              covertype_p2_one_shot_llk_ratio_min,
                                              covertype_p2_uniform_llk_ratio_min),
                                        label=c(rep("2-probit",length(frac)),rep("one-shot",length(frac)),rep("uniform",length(frac))),
                                        frac=c(rep(frac,3)))

p2_llk_ratio_covertype_data <- data.frame(mean=c(log(covertype_p2_2probit_llk_ratio_average),
                                                 log(covertype_p2_one_shot_llk_ratio_average),
                                                 log(covertype_p2_uniform_llk_ratio_average)),
                                          median=c(log(covertype_p2_2probit_llk_ratio_median),
                                                   log(covertype_p2_one_shot_llk_ratio_median),
                                                   log(covertype_p2_uniform_llk_ratio_median)),
                                          max=c(log(covertype_p2_2probit_llk_ratio_max),
                                                log(covertype_p2_one_shot_llk_ratio_max),
                                                log(covertype_p2_uniform_llk_ratio_max)),
                                          min=c(log(covertype_p2_2probit_llk_ratio_min),
                                                log(covertype_p2_one_shot_llk_ratio_min),
                                                log(covertype_p2_uniform_llk_ratio_min)),
                                          label=c(rep("2-probit",length(frac)),rep("one-shot",length(frac)),rep("uniform",length(frac))),
                                          frac=c(rep(frac,3)))


jpeg("cover_type_p2_llk_ratio.jpeg", units="in", width=8, height=5, res=300)

ggplot(p2_llk_ratio_covertype_data, aes(frac, median, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="likelihood ratio", colour="method", fill="method") +
  theme_bw()

dev.off()


#####MMD#######
covertype_p_2_MMD_median <- c()
covertype_p_2_MMD_max <- c()
covertype_p_2_MMD_min <- c()
covertype_p_2_MMD_average <- c()

for (i in 1:length(frac)){
  print(paste("2-probit_MMD_Frac:",i,sep = ""))
  medianchain <- c()
  for (j in 1:10) {
    b <- mmd(sample_x=tail(get(paste("covertype_p_2_frac",frac[i],"chain_",j,sep = ""))$chain,100),
             sample_y = tail(cover_type_2$chain,100))
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN=2, mean)
  covertype_p_2_MMD_median <- c(covertype_p_2_MMD_median,median_result)
  covertype_p_2_MMD_max <- c(covertype_p_2_MMD_max,maxresult)
  covertype_p_2_MMD_min <- c(covertype_p_2_MMD_min,minresult)
  covertype_p_2_MMD_average <- c(covertype_p_2_MMD_average,averageresult)
}

covertype_p_2_lpcoreset_MMD_median <- c()
covertype_p_2_lpcoreset_MMD_max <- c()
covertype_p_2_lpcoreset_MMD_min <- c()
covertype_p_2_lpcoreset_MMD_average <- c()
for (i in 1:length(frac)){
  print(paste("rootl2_MMD_Frac:",i,sep = ""))
  medianchain <- c()
  for (j in 1:10) {
    b <- mmd(sample_x=tail(get(paste("covertype_p_2_lpcoreset_frac",frac[i],"chain_",j,sep = ""))$chain,100),
             sample_y = tail(cover_type_2$chain,100))
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN=2, mean)
  covertype_p_2_lpcoreset_MMD_median <- c(covertype_p_2_lpcoreset_MMD_median,median_result)
  covertype_p_2_lpcoreset_MMD_max <- c(covertype_p_2_lpcoreset_MMD_max,maxresult)
  covertype_p_2_lpcoreset_MMD_min <- c(covertype_p_2_lpcoreset_MMD_min,minresult)
  covertype_p_2_lpcoreset_MMD_average <- c(covertype_p_2_lpcoreset_MMD_average,averageresult)
}

covertype_p_2_uniform_MMD_median <- c()
covertype_p_2_uniform_MMD_max <- c()
covertype_p_2_uniform_MMD_min <- c()
covertype_p_2_uniform_MMD_average <- c()
for (i in 1:length(frac)){
  print(paste("uniform_MMD_Frac:",i,sep = ""))
  medianchain <- c()
  for (j in 1:5) {
    b <- mmd(sample_x=tail(get(paste("covertype_p_2_uniform_frac",frac[i],"chain_",j,sep = ""))$chain,100),
             sample_y = tail(cover_type_2$chain,100))
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN=2, mean)
  covertype_p_2_uniform_MMD_median <- c(covertype_p_2_uniform_MMD_median,median_result)
  covertype_p_2_uniform_MMD_max <- c(covertype_p_2_uniform_MMD_max,maxresult)
  covertype_p_2_uniform_MMD_min <- c(covertype_p_2_uniform_MMD_min,minresult)
  covertype_p_2_uniform_MMD_average <- c(covertype_p_2_uniform_MMD_average,averageresult)
}

p2_mmd_covertype_data <- data.frame(mean=c(covertype_p_2_MMD_average,
                                           covertype_p_2_lpcoreset_MMD_average,
                                           covertype_p_2_uniform_MMD_average),
                                         median=c(covertype_p_2_MMD_median,
                                                  covertype_p_2_lpcoreset_MMD_median,
                                                  covertype_p_2_uniform_MMD_median),
                                         max=c(covertype_p_2_MMD_max,
                                               covertype_p_2_lpcoreset_MMD_max,
                                               covertype_p_2_uniform_MMD_max),
                                         min=c(covertype_p_2_MMD_min,
                                               covertype_p_2_lpcoreset_MMD_min,
                                               covertype_p_2_uniform_MMD_min),
                                         label=c(rep("2-probit",length(frac)),rep("root l2",length(frac)),rep("uniform",length(frac))),
                                         frac=c(rep(frac,3)))
p2_mmd_covertype_data
ggplot(p2_mmd_covertype_data, aes(frac, mean, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="MMD", colour="method", fill="method") +
  theme_bw()
###########
###########
########fixed p=1#########
for (i in 1:length(frac)) {
  for (j in 1:5) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Compute_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=1)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y[sample_index],
                             initial_theda=probit_coefficients,
                             true_theta=probit_coefficients,Lp=1)
    assign(paste("covertype_p_1_frac",frac[i],"chain_",j,sep = ""), Model)
  }
}
#####LP sampling#########
for (i in 1:length(frac)) {
  for (j in 1:5) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Lp_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=1)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y[sample_index],
                             initial_theda=probit_coefficients,
                             true_theta=probit_coefficients,Lp=1)
    assign(paste("covertype_p_1_lpcoreset_frac",frac[i],"chain_",j,sep = ""), Model)
  }
}
#########uniform sampling##############
for (i in 1:length(frac)) {
  for (j in 1:5) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Uniform_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=1)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y[sample_index],
                             initial_theda=probit_coefficients,
                             true_theta=probit_coefficients,Lp=1)
    assign(paste("covertype_p_1_uniform_frac",frac[i],"chain_",j,sep = ""), Model)
  }
}

for (i in 1:length(frac)) {
  for (j in 1:5) {
    coreset_size<- N*frac[i]
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    scores<- one_shot_coreset(sketch_size = sketch_size, coreset_size = coreset_size, X=X, Lp_max = 2)
    sample_index <- sample(1:N,size = coreset_size,
                           replace = FALSE, prob = scores)
    Model <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                             X=X[sample_index,], y=y[sample_index],
                             initial_theda=probit_coefficients,
                             true_theta=probit_coefficients,Lp=1)
    assign(paste("covertype_p_1_one_shot",frac[i],"chain_",j,sep = ""), Model)
  }
}

########################################

############lilelihood ratio###############


covertype_p1_1probit_llk_ratio_median <- c()
covertype_p1_1probit_llk_ratio_max <- c()
covertype_p1_1probit_llk_ratio_min <- c()
covertype_p1_1probit_llk_ratio_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    b<- median(abs(llk(X, cover_type_1$post.mean,1)/
                  llk(X,get(paste("covertype_p_1_frac",frac[i],"chain_",j,sep = ""))$post.mean,1)))
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  covertype_p1_1probit_llk_ratio_median <- c(covertype_p1_1probit_llk_ratio_median,median_result)
  covertype_p1_1probit_llk_ratio_max <- c(covertype_p1_1probit_llk_ratio_max,maxresult)
  covertype_p1_1probit_llk_ratio_min <- c(covertype_p1_1probit_llk_ratio_min,minresult)
  covertype_p1_1probit_llk_ratio_average <- c(covertype_p1_1probit_llk_ratio_average,meanresult)
}

covertype_p1_uniform_llk_ratio_median <- c()
covertype_p1_uniform_llk_ratio_max <- c()
covertype_p1_uniform_llk_ratio_min <- c()
covertype_p1_uniform_llk_ratio_average <- c()


for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b<- median(abs(llk(X, cover_type_1$post.mean,1)/
                  llk(X,get(paste("covertype_p_1_uniform_frac",frac[i],"chain_",j,sep = ""))$post.mean,1)))
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  covertype_p1_uniform_llk_ratio_median <- c(covertype_p1_uniform_llk_ratio_median,median_result)
  covertype_p1_uniform_llk_ratio_max <- c(covertype_p1_uniform_llk_ratio_max,maxresult)
  covertype_p1_uniform_llk_ratio_min <- c(covertype_p1_uniform_llk_ratio_min,minresult)
  covertype_p1_uniform_llk_ratio_average <- c(covertype_p1_uniform_llk_ratio_average,meanresult)
}

covertype_p1_one_shot_llk_ratio_median <- c()
covertype_p1_one_shot_llk_ratio_max <- c()
covertype_p1_one_shot_llk_ratio_min <- c()
covertype_p1_one_shot_llk_ratio_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    b<- median(abs(llk(X, cover_type_1$post.mean,1)/
                  llk(X,get(paste("covertype_p_1_one_shot",frac[i],"chain_",j,sep = ""))$post.mean,1)))
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  covertype_p1_one_shot_llk_ratio_median <- c(covertype_p1_one_shot_llk_ratio_median,median_result)
  covertype_p1_one_shot_llk_ratio_max <- c(covertype_p1_one_shot_llk_ratio_max,maxresult)
  covertype_p1_one_shot_llk_ratio_min <- c(covertype_p1_one_shot_llk_ratio_min,minresult)
  covertype_p1_one_shot_llk_ratio_average <- c(covertype_p1_one_shot_llk_ratio_average,meanresult)
}

covertype_p1_rootl2_llk_ratio_median <- c()
covertype_p1_rootl2_llk_ratio_max <- c()
covertype_p1_rootl2_llk_ratio_min <- c()
covertype_p1_rootl2_llk_ratio_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    b<- median(abs(llk(X, cover_type_1$post.mean,1)-
              llk(X,get(paste("covertype_p_1_lpcoreset_frac",frac[i],"chain_",j,sep = ""))$post.mean,1)))
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  covertype_p1_rootl2_llk_ratio_median <- c(covertype_p1_rootl2_llk_ratio_median,median_result)
  covertype_p1_rootl2_llk_ratio_max <- c(covertype_p1_rootl2_llk_ratio_max,maxresult)
  covertype_p1_rootl2_llk_ratio_min <- c(covertype_p1_rootl2_llk_ratio_min,minresult)
  covertype_p1_rootl2_llk_ratio_average <- c(covertype_p1_rootl2_llk_ratio_average,meanresult)
}

p1_llk_ratio_covertype_data <- data.frame(mean=c(covertype_p1_1probit_llk_ratio_average,
                                                 covertype_p1_one_shot_llk_ratio_average,
                                                 covertype_p1_uniform_llk_ratio_average,
                                                 covertype_p1_rootl2_llk_ratio_average),
                                          median=c(covertype_p1_1probit_llk_ratio_median,
                                                   covertype_p1_one_shot_llk_ratio_median,
                                                   covertype_p1_uniform_llk_ratio_median,
                                                   covertype_p1_rootl2_llk_ratio_median),
                                          max=c(covertype_p1_1probit_llk_ratio_max,
                                                covertype_p1_one_shot_llk_ratio_max,
                                                covertype_p1_uniform_llk_ratio_max,
                                                covertype_p1_rootl2_llk_ratio_max),
                                          min=c(covertype_p1_1probit_llk_ratio_min,
                                                covertype_p1_one_shot_llk_ratio_min,
                                                covertype_p1_uniform_llk_ratio_min,
                                                covertype_p1_rootl2_llk_ratio_min),
                                          label=c(rep("1-probit",length(frac)),rep("one-shot",length(frac)),rep("uniform",length(frac)),
                                                  rep("rootl2",length(frac))),
                                          frac=c(rep(frac,4)))

p1_llk_ratio_covertype_data <- data.frame(mean=c(log(covertype_p1_1probit_llk_ratio_average),
                                                 log(covertype_p1_one_shot_llk_ratio_average),
                                                 log(covertype_p1_uniform_llk_ratio_average)),
                                          median=c(log(covertype_p1_1probit_llk_ratio_median),
                                                   log(covertype_p1_one_shot_llk_ratio_median),
                                                   log(covertype_p1_uniform_llk_ratio_median)),
                                          max=c(log(covertype_p1_1probit_llk_ratio_max),
                                                log(covertype_p1_one_shot_llk_ratio_max),
                                                log(covertype_p1_uniform_llk_ratio_max)),
                                          min=c(log(covertype_p1_1probit_llk_ratio_min),
                                                log(covertype_p1_one_shot_llk_ratio_min),
                                                log(covertype_p1_uniform_llk_ratio_min)),
                                          label=c(rep("1-probit",length(frac)),rep("one-shot",length(frac)),rep("uniform",length(frac))),
                                          frac=c(rep(frac,3)))


p1_llk_ratio_covertype_data <- data.frame(mean=c(covertype_p1_1probit_llk_ratio_average,
                                                 covertype_p1_one_shot_llk_ratio_average,
                                                 covertype_p1_uniform_llk_ratio_average),
                                          median=c(covertype_p1_1probit_llk_ratio_median,
                                                   covertype_p1_one_shot_llk_ratio_median,
                                                   covertype_p1_uniform_llk_ratio_median),
                                          max=c(covertype_p1_1probit_llk_ratio_max,
                                                covertype_p1_one_shot_llk_ratio_max,
                                                covertype_p1_uniform_llk_ratio_max),
                                          min=c(covertype_p1_1probit_llk_ratio_min,
                                                covertype_p1_one_shot_llk_ratio_min,
                                                covertype_p1_uniform_llk_ratio_min),
                                          label=c(rep("1-probit",length(frac)),rep("one-shot",length(frac)),rep("uniform",length(frac))),
                                          frac=c(rep(frac,3)))

# Xb <- X%*%cover_type_1$post.mean
# (2*y-1)*Xb==a%*%cover_type_1$post.mean
# 
# a%*%cover_type_1$post.mean==X%*%cover_type_1$post.mean
# y[990:1000]
# getwd()
jpeg("cover_type_p1_llk_ratio.jpeg", units="in", width=8, height=5, res=300)

ggplot(p1_llk_ratio_covertype_data, aes(frac, median, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="Bayes factor", colour="method", fill="method") +
  ggtitle("(log)-Bayes Factor of covertype data (p=1)")+coord_cartesian(ylim=c(0, 2.5))+
  theme_bw()

 dev.off()

#############post mean difference######
covertype_p_1_frac_mean_diff_norm_median <- c()
covertype_p_1_frac_mean_diff_norm_max <- c()
covertype_p_1_frac_mean_diff_norm_min <- c()
covertype_p_1_frac_mean_diff_norm_average <- c()
cover_type_1$post.mean
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    b <- norm(get(paste("covertype_p_1_frac",frac[i],"chain_",j,sep = ""))$post.mean-cover_type_1$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  covertype_p_1_frac_mean_diff_norm_median <- c(covertype_p_1_frac_mean_diff_norm_median,median_result)
  covertype_p_1_frac_mean_diff_norm_max <- c(covertype_p_1_frac_mean_diff_norm_max,maxresult)
  covertype_p_1_frac_mean_diff_norm_min <- c(covertype_p_1_frac_mean_diff_norm_min,minresult)
  covertype_p_1_frac_mean_diff_norm_average <- c(covertype_p_1_frac_mean_diff_norm_average,meanresult)
  
}
##########lpcoreset######
covertype_p_1_lpcoreset_frac_mean_diff_norm_median <- c()
covertype_p_1_lpcoreset_frac_mean_diff_norm_max <- c()
covertype_p_1_lpcoreset_frac_mean_diff_norm_min <- c()
covertype_p_1_lpcoreset_frac_mean_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    b <- norm(get(paste("covertype_p_2_lpcoreset_frac",frac[i],"chain_",j,sep = ""))$post.mean-cover_type_1$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  covertype_p_1_lpcoreset_frac_mean_diff_norm_median <- c(covertype_p_1_lpcoreset_frac_mean_diff_norm_median,median_result)
  covertype_p_1_lpcoreset_frac_mean_diff_norm_max <- c(covertype_p_1_lpcoreset_frac_mean_diff_norm_max,maxresult)
  covertype_p_1_lpcoreset_frac_mean_diff_norm_min <- c(covertype_p_1_lpcoreset_frac_mean_diff_norm_min,minresult)
  covertype_p_1_lpcoreset_frac_mean_diff_norm_average <- c(covertype_p_1_lpcoreset_frac_mean_diff_norm_average,meanresult)
  
}
#######uniform sampling###########
covertype_p_1_uniform_frac_mean_diff_norm_median <- c()
covertype_p_1_uniform_frac_mean_diff_norm_max <- c()
covertype_p_1_uniform_frac_mean_diff_norm_min <- c()
covertype_p_1_uniform_frac_mean_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    b <- norm(get(paste("covertype_p_1_uniform_frac",frac[i],"chain_",j,sep = ""))$post.mean-cover_type_1$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  covertype_p_1_uniform_frac_mean_diff_norm_median <- c(covertype_p_1_uniform_frac_mean_diff_norm_median,median_result)
  covertype_p_1_uniform_frac_mean_diff_norm_max <- c(covertype_p_1_uniform_frac_mean_diff_norm_max,maxresult)
  covertype_p_1_uniform_frac_mean_diff_norm_min <- c(covertype_p_1_uniform_frac_mean_diff_norm_min,minresult)
  covertype_p_1_uniform_frac_mean_diff_norm_average <- c(covertype_p_1_uniform_frac_mean_diff_norm_average,meanresult)
  
}

covertype_p_1_one_shot_mean_diff_norm_median <- c()
covertype_p_1_one_shot_mean_diff_norm_max <- c()
covertype_p_1_one_shot_mean_diff_norm_min <- c()
covertype_p_1_one_shot_mean_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    
    b <- norm(get(paste("covertype_p_1_one_shot",frac[i],"chain_",j,sep = ""))$post.mean-cover_type_1$post.mean,type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  covertype_p_1_one_shot_mean_diff_norm_median <- c(covertype_p_1_one_shot_mean_diff_norm_median,median_result)
  covertype_p_1_one_shot_mean_diff_norm_max <- c(covertype_p_1_one_shot_mean_diff_norm_max,maxresult)
  covertype_p_1_one_shot_mean_diff_norm_min <- c(covertype_p_1_one_shot_mean_diff_norm_min,minresult)
  covertype_p_1_one_shot_mean_diff_norm_average <- c(covertype_p_1_one_shot_mean_diff_norm_average,meanresult)
}
p1_meandiff_covertype_data
p1_meandiff_covertype_data <- data.frame(mean=c(covertype_p_1_frac_mean_diff_norm_average,
                                               covertype_p_1_lpcoreset_frac_mean_diff_norm_average,
                                               covertype_p_1_uniform_frac_mean_diff_norm_average,
                                               covertype_p_1_one_shot_mean_diff_norm_average),
                                        median=c(covertype_p_1_frac_mean_diff_norm_median,
                                                 covertype_p_1_lpcoreset_frac_mean_diff_norm_median,
                                                 covertype_p_1_uniform_frac_mean_diff_norm_median,
                                                 covertype_p_1_one_shot_mean_diff_norm_median),
                                        max=c(covertype_p_1_frac_mean_diff_norm_max,
                                              covertype_p_1_lpcoreset_frac_mean_diff_norm_max,
                                              covertype_p_1_uniform_frac_mean_diff_norm_max,
                                              covertype_p_1_one_shot_mean_diff_norm_max),
                                        min=c(covertype_p_1_frac_mean_diff_norm_min,
                                              covertype_p_1_lpcoreset_frac_mean_diff_norm_min,
                                              covertype_p_1_uniform_frac_mean_diff_norm_min,
                                              covertype_p_1_one_shot_mean_diff_norm_min),
                                        label=c(rep("1-probit",length(frac)),rep("root l2",length(frac)),rep("uniform",length(frac)),
                                                rep("one-shot",length(frac))),
                                        frac=c(rep(frac,4)))
p1_meandiff_covertype_data$label
jpeg("cover_type_p1_mean.jpeg", units="in", width=8, height=5, res=300)

ggplot(p1_meandiff_covertype_data, aes(frac, mean, fill=label, colour=label)) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior mean", colour="method", fill="method") +
  theme_bw()+
  ggtitle("Covertype data (p=1)")
dev.off()
###############post cov diff##########
######poserior covariance difference
covertype_p1_pos_cov_diff_norm_median <- c()
covertype_p1_pos_cov_diff_norm_max <- c()
covertype_p1_pos_cov_diff_norm_min <- c()
covertype_p1_pos_cov_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("covertype_p_1_frac",frac[i],"chain_",j,sep = ""))$chain)-
                cov(cover_type_1$chain),type = "2")
    
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN = 2, mean)
  covertype_p1_pos_cov_diff_norm_median <- c(covertype_p1_pos_cov_diff_norm_median,median_result)
  covertype_p1_pos_cov_diff_norm_max <- c(covertype_p1_pos_cov_diff_norm_max,maxresult)
  covertype_p1_pos_cov_diff_norm_min <- c(covertype_p1_pos_cov_diff_norm_min,minresult)
  covertype_p1_pos_cov_diff_norm_average <- c(covertype_p1_pos_cov_diff_norm_average,averageresult)
  
  
}
####lpcoreset######

covertype_p1_lpcoreset_pos_cov_diff_norm_median <- c()
covertype_p1_lpcoreset_pos_cov_diff_norm_max <- c()
covertype_p1_lpcoreset_pos_cov_diff_norm_min <- c()
covertype_p1_lpcoreset_pos_cov_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("covertype_p_1_lpcoreset_frac",frac[i],"chain_",j,sep = ""))$chain)-
                cov(cover_type_1$chain),type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN = 2, mean)
  covertype_p1_lpcoreset_pos_cov_diff_norm_median <- c(covertype_p1_lpcoreset_pos_cov_diff_norm_median,median_result)
  covertype_p1_lpcoreset_pos_cov_diff_norm_max <- c(covertype_p1_lpcoreset_pos_cov_diff_norm_max,maxresult)
  covertype_p1_lpcoreset_pos_cov_diff_norm_min <- c(covertype_p1_lpcoreset_pos_cov_diff_norm_min,minresult)
  covertype_p1_lpcoreset_pos_cov_diff_norm_average <- c(covertype_p1_lpcoreset_pos_cov_diff_norm_average,averageresult)
  
}

covertype_p1_uniform_pos_cov_diff_norm_median <- c()
covertype_p1_uniform_pos_cov_diff_norm_max <- c()
covertype_p1_uniform_pos_cov_diff_norm_min <- c()
covertype_p1_uniform_pos_cov_diff_norm_average <- c()

for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("covertype_p_1_uniform_frac",frac[i],"chain_",j,sep = ""))$chain)-
                cov(cover_type_1$chain),type = "2")
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN = 2, mean)
  covertype_p1_uniform_pos_cov_diff_norm_median <- c(covertype_p1_uniform_pos_cov_diff_norm_median,median_result)
  covertype_p1_uniform_pos_cov_diff_norm_max <- c(covertype_p1_uniform_pos_cov_diff_norm_max,maxresult)
  covertype_p1_uniform_pos_cov_diff_norm_min <- c(covertype_p1_uniform_pos_cov_diff_norm_min,minresult)
  covertype_p1_uniform_pos_cov_diff_norm_average <- c(covertype_p1_uniform_pos_cov_diff_norm_average,averageresult)
}


covertype_p_1_one_shot_cov_diff_norm_median <- c()
covertype_p_1_one_shot_cov_diff_norm_max <- c()
covertype_p_1_one_shot_cov_diff_norm_min <- c()
covertype_p_1_one_shot_cov_diff_norm_average <- c()
for (i in 1:length(frac)){
  medianchain <- c()
  for (j in 1:5) {
    b <- norm(cov(get(paste("covertype_p_1_one_shot",frac[i],"chain_",j,sep = ""))$chain)-
                cov(cover_type_1$chain),type = "2")
        medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  meanresult <- apply(medianchain, MARGIN = 2, mean)
  covertype_p_1_one_shot_cov_diff_norm_median <- c(covertype_p_1_one_shot_cov_diff_norm_median,median_result)
  covertype_p_1_one_shot_cov_diff_norm_max <- c(covertype_p_1_one_shot_cov_diff_norm_max,maxresult)
  covertype_p_1_one_shot_cov_diff_norm_min <- c(covertype_p_1_one_shot_cov_diff_norm_min,minresult)
  covertype_p_1_one_shot_cov_diff_norm_average <- c(covertype_p_1_one_shot_cov_diff_norm_average,meanresult)
}

p1_covdiff_covertype_data <- data.frame(mean=c(covertype_p1_pos_cov_diff_norm_average,
                                                covertype_p1_lpcoreset_pos_cov_diff_norm_average,
                                                covertype_p1_uniform_pos_cov_diff_norm_average,
                                               covertype_p_1_one_shot_cov_diff_norm_average),
                                         median=c(covertype_p1_pos_cov_diff_norm_median,
                                                  covertype_p1_lpcoreset_pos_cov_diff_norm_median,
                                                  covertype_p1_uniform_pos_cov_diff_norm_median,
                                                  covertype_p_1_one_shot_cov_diff_norm_median),
                                         max=c(covertype_p1_pos_cov_diff_norm_max,
                                               covertype_p1_lpcoreset_pos_cov_diff_norm_max,
                                               covertype_p1_uniform_pos_cov_diff_norm_max,
                                               covertype_p_1_one_shot_cov_diff_norm_max),
                                         min=c(covertype_p1_pos_cov_diff_norm_min,
                                               covertype_p1_lpcoreset_pos_cov_diff_norm_min,
                                               covertype_p1_uniform_pos_cov_diff_norm_min,
                                               covertype_p_1_one_shot_cov_diff_norm_min),
                                         label=c(rep("1-probit",length(frac)),rep("root l2",length(frac)),rep("uniform",length(frac)),
                                                 rep("one-shot",length(frac))),
                                         frac=c(rep(frac,4)))

jpeg("cover_type_p1_cov.jpeg", units="in", width=8, height=5, res=300)

ggplot(p1_covdiff_covertype_data, aes(frac, median, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior covariance", colour="method", fill="method") +
  theme_bw()+
  ggtitle("Covertype data (p=1)")
dev.off()

#####MMD#######
covertype_p_1_MMD_median <- c()
covertype_p_1_MMD_max <- c()
covertype_p_1_MMD_min <- c()
covertype_p_1_MMD_average <- c()

for (i in 1:length(frac)){
  print(paste("1-probit_MMD_Frac:",i,sep = ""))
  medianchain <- c()
  for (j in 1:10) {
    b <- mmd(sample_x=tail(get(paste("covertype_p_1_frac",frac[i],"chain_",j,sep = ""))$chain,100),
             sample_y = tail(cover_type_1$chain,100))
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN=2, mean)
  covertype_p_1_MMD_median <- c(covertype_p_1_MMD_median,median_result)
  covertype_p_1_MMD_max <- c(covertype_p_1_MMD_max,maxresult)
  covertype_p_1_MMD_min <- c(covertype_p_1_MMD_min,minresult)
  covertype_p_1_MMD_average <- c(covertype_p_1_MMD_average,averageresult)
}

covertype_p_1_lpcoreset_MMD_median <- c()
covertype_p_1_lpcoreset_MMD_max <- c()
covertype_p_1_lpcoreset_MMD_min <- c()
covertype_p_1_lpcoreset_MMD_average <- c()
for (i in 1:length(frac)){
  print(paste("rootl2_MMD_Frac:",i,sep = ""))
  medianchain <- c()
  for (j in 1:10) {
    b <- mmd(sample_x=tail(get(paste("covertype_p_1_lpcoreset_frac",frac[i],"chain_",j,sep = ""))$chain,100),
             sample_y = tail(cover_type_1$chain,100))
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN=2, mean)
  covertype_p_1_lpcoreset_MMD_median <- c(covertype_p_1_lpcoreset_MMD_median,median_result)
  covertype_p_1_lpcoreset_MMD_max <- c(covertype_p_1_lpcoreset_MMD_max,maxresult)
  covertype_p_1_lpcoreset_MMD_min <- c(covertype_p_1_lpcoreset_MMD_min,minresult)
  covertype_p_1_lpcoreset_MMD_average <- c(covertype_p_1_lpcoreset_MMD_average,averageresult)
}

covertype_p_1_uniform_MMD_median <- c()
covertype_p_1_uniform_MMD_max <- c()
covertype_p_1_uniform_MMD_min <- c()
covertype_p_1_uniform_MMD_average <- c()
for (i in 1:length(frac)){
  print(paste("uniform_MMD_Frac:",i,sep = ""))
  medianchain <- c()
  for (j in 1:5) {
    b <- mmd(sample_x=tail(get(paste("covertype_p_1_uniform_frac",frac[i],"chain_",j,sep = ""))$chain,100),
             sample_y = tail(cover_type_1$chain,100))
    medianchain<- rbind(medianchain,b)
  }
  median_result <- apply(medianchain, MARGIN = 2, median)
  maxresult <- apply(medianchain, MARGIN = 2, max)
  minresult <- apply(medianchain, MARGIN = 2, min)
  averageresult <- apply(medianchain, MARGIN=2, mean)
  covertype_p_1_uniform_MMD_median <- c(covertype_p_1_uniform_MMD_median,median_result)
  covertype_p_1_uniform_MMD_max <- c(covertype_p_1_uniform_MMD_max,maxresult)
  covertype_p_1_uniform_MMD_min <- c(covertype_p_1_uniform_MMD_min,minresult)
  covertype_p_1_uniform_MMD_average <- c(covertype_p_1_uniform_MMD_average,averageresult)
}

p1_mmd_covertype_data <- data.frame(mean=c(covertype_p_1_MMD_average,
                                           covertype_p_1_lpcoreset_MMD_average,
                                           covertype_p_1_uniform_MMD_average),
                                    median=c(covertype_p_1_MMD_median,
                                             covertype_p_1_lpcoreset_MMD_median,
                                             covertype_p_1_uniform_MMD_median),
                                    max=c(covertype_p_1_MMD_max,
                                          covertype_p_1_lpcoreset_MMD_max,
                                          covertype_p_1_uniform_MMD_max),
                                    min=c(covertype_p_1_MMD_min,
                                          covertype_p_1_lpcoreset_MMD_min,
                                          covertype_p_1_uniform_MMD_min),
                                    label=c(rep("1-probit",length(frac)),rep("root l2",length(frac)),rep("uniform",length(frac))),
                                    frac=c(rep(frac,3)))
p1_mmd_covertype_data
ggplot(p1_meandiff_covertype_data, aes(frac, median, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior mean", colour="method", fill="method") +
  theme_bw()
ggplot(p1_covdiff_covertype_data, aes(frac, median, fill=factor(label), colour=factor(label))) +
  geom_ribbon(aes(ymin=min, max=max), alpha=0.2, colour=NA) +
  geom_line() +
  labs(x="fraction", y="norm difference of posterior covariance", colour="method", fill="method") +
  theme_bw()
ggplot(p1_mmd_covertype_data, aes(frac, median, fill=factor(label), colour=factor(label))) +
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
    assign(paste("covertype_p2_nofix_frac",frac[i],sep = ""), Model)
  
}

for (i in 1:length(frac)) {
  
  print(paste("Model_Frac:",i,"chain:",j,sep = ""))
  sample_index<- Compute_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=1)
  Model <- multi_chain(N_sim=N_sim,burn_in=burn_in,
                       X=X[sample_index,], y=y[sample_index],initial_theda=probit_coefficients,
                       true_theta=probit_coefficients,Lp=3,M_iter=100,range=c(0.1,5),step=0.02,times = 1)
  assign(paste("covertype_p1_nofix_frac",frac[i],sep = ""), Model)
  
}
#####LP sampling#########
for (i in 1:length(frac)) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Lp_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=2)
    Model <- multi_chain(N_sim=N_sim,burn_in=burn_in,
                         X=X[sample_index,], y=y[sample_index],initial_theda=probit_coefficients,
                         true_theta=probit_coefficients,Lp=3,M_iter=100,range=c(0.1,5),step=0.02,times = 3)
    assign(paste("covertype_p_lpcoreset_nofix_frac",frac[i],sep = ""), Model)
  
}
#########uniform sampling
for (i in 1:length(frac)) {
    print(paste("Model_Frac:",i,"chain:",j,sep = ""))
    sample_index<- Uniform_coreset(sketch_size = sketch_size,coreset_size = N*frac[i],X=X,Lp=2)
    Model <- multi_chain(N_sim=N_sim,burn_in=burn_in,
                         X=X[sample_index,], y=y[sample_index],initial_theda=probit_coefficients,
                         true_theta=probit_coefficients,Lp=3,M_iter=100,range=c(0.1,5),step=0.02,times = 3)
    assign(paste("covertype_p_uniform_nofix_frac",frac[i],sep = ""), Model)
}

covertype_p_summary1 <- c()
for (i in 1:length(frac)){
  
  b <- get(paste("covertype_p2_nofix_frac",frac[i],sep = ""))$Lp.mean
  
  covertype_p_summary1<- rbind(covertype_p_summary1,b)
  
}
plot(frac,covertype_p_summary1,type = "l")

covertype_p_summary2 <- c()
for (i in 1:length(frac)){
  
  b <- get(paste("covertype_p_lpcoreset_nofix_frac",frac[i],sep = ""))$Lp.mean
  
  covertype_p_summary2<- rbind(covertype_p_summary2,b)
  
}
plot(frac,covertype_p_summary2,type = "l")

covertype_p_summary3 <- c()
for (i in 1:length(frac)){
  
  b <- get(paste("covertype_p_uniform_nofix_frac",frac[i],sep = ""))$Lp.mean
  
  covertype_p_summary3<- rbind(covertype_p_summary3,b)
  
}
plot(frac,covertype_p_summary2,type = "l")
lines(frac,covertype_p_summary3,col=2)

plot(frac,covertype_p_summary1,type = "l",ylim = c(0,3),xlab = "sample fraction", ylab = "estimated P")
lines(frac,covertype_p_summary2,col=2)
lines(frac,covertype_p_summary3,col=3)
legend("topright",legend=c("2-probit","root l2","uniform"), col=c(1,2,3), 
       bty = 'n', lwd = 2, inset = c(0, 0), lty = 1, cex = 0.73)
