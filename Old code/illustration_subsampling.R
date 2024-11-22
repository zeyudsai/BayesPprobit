set.seed(2)
rm(list = ls())
options(scipen = 200)


setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("Functions.R")

N <- 100000
#D should always be odd number since we use (D-1)/2 to generate the variables
D <- 3
N_sim <- 1000
burn_in<- 800
#gengeratedata and true parameterye

#X <- generate.data.cauthy(N,D)

X1 <- generate.data1(N,D)
X0 <- generate.data0(N,D)
X <- rbind(X1,X0)
y1 <- rep(1,N/2)
y0 <- rep(0,N/2)
label <- sample(1:N,N,replace = FALSE)
y <- c(y1,y0)
#Randomly setting some misclassification points
error_sample <- sample(N,size = 200,replace = FALSE)
y[error_sample][y[error_sample] == 1] <- 2
y[error_sample][y[error_sample] == 0] <- 1
y[error_sample][y[error_sample] == 2] <- 0

X <- X[label,]
y <- y[label]
data <- data.frame(cbind(X,y))
#Randomly setting some errors
noisy <- matrix(c(runif(10,40,50),runif(10,-45,-35),runif(20,-10,0),runif(10,40,50)
                  ,runif(10,0,10),runif(10,35,45),runif(20,-55,-45),runif(10,-35,-25), 
                  runif(50*(D-3),-10,10)),
                ncol=D-1)
noise <- data.frame(rep(1,50),noisy,c(rep(1,20),rep(0,30)))
names(noise) <- names(data)
X <- rbind(X,cbind(rep(1,50),noisy))
data <- rbind(data,noise)
#index1 <- which(data$V2>5&data$V3>5)
#index1
#index2 <- which(data$V2<(-5) & data$V3 < (-5) & data$V4 < (-5) & data$V5 < (-5))
#index2
#error_sample1 <- sample(index1,size = length(index1))
#error_sample2 <- sample(index2,size = length(index2)/2)
#data$y[error_sample1] <- 0
#error_sample1
#data$y[error_sample2] <-1
#error_sample2
data$label <- 0
data$label[data$y==0]=0
data$label[data$y==1]=1
#Compute coreset
sample_index<- Compute_coreset(sketch_size = 200,coreset_size = 200,X=X,Lp=2)



#Using MLE estimator to otain the maximum likelihood result
#Assign a as a symbol states NxxKDxx,which represents the sample size and number of variables.

MLE_model <- glm(data$y~X[,-1],family = binomial(link = probit))
MLE_logit <- glm(data$y~X[,-1],family = binomial(link = logit))

MLE_coefficients <- MLE_logit$coefficients
#XB <- X%*%MLE_coefficients
#plot(XB,y)
P_probit <- multi_chain(N_sim=N_sim,burn_in=burn_in,
                        X=X, y=data$y,initial_theda=MLE_model$coefficients,
                        true_theta=MLE_coefficients,Lp=5,M_iter=200,
                        range=c(0.1,5),step=0.05,times = 1)
P_probit
P_original_L1 <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                                X=X, y=data$y,initial_theda=MLE_model$coefficients,
                                true_theta=MLE_coefficients,Lp=1)
P_original_L1
P_original_L2 <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                                X=X, y=data$y,initial_theda=MLE_model$coefficients,
                                true_theta=MLE_coefficients,Lp=2)
P_original_L2

P_probit$Beta.mean
P_probit
P_original_L1
P_original_L2
MLE_logit
MLE_model
MLE_coefficients
plot(1:N_sim,P_probit$Lp_chain[[1]],type = "l")

frac <-seq(0.0001,0.005,0.0001)
frac[10]
coreset_size <- N*frac[10]
coreset_size
sketch_size <- D^2
par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=TRUE)





sample_index <- c()
for (i in 1:10) {
  index1<- Compute_coreset(sketch_size = sketch_size,coreset_size = N*frac[10],X=X,Lp=1)
  sample_index <- union(sample_index,index1)
}

sample_index1 <- c()
for (i in 1:10) {
  index1<- Compute_coreset(sketch_size = sketch_size,coreset_size = N*frac[10],X=X,Lp=2)
  sample_index1 <- union(sample_index1,index1)
}

sample_index2 <- c()
for (i in 1:10) {
  index1<- Uniform_coreset(sketch_size = sketch_size,coreset_size = N*frac[10],X=X,Lp=2)
  sample_index2 <- union(sample_index2,index1)
}

N <- nrow(X)
sample_index4 <- c()
for (i in 1:10) {
  scores<- one_shot_coreset(sketch_size = sketch_size, coreset_size = N*frac[10], X=X, Lp_max = 2.5)
  index1 <- sample(1:N,size = N*frac[10],
                         replace = FALSE, prob = scores)
  sample_index4 <- union(sample_index4,index1)
}
sample_index5 <- c()
for (i in 1:10) {
  index1<- Lp_coreset(sketch_size = sketch_size,coreset_size = N*frac[10],X=X,Lp=2)
  sample_index5 <- union(sample_index5,index1)
}


par(mfrow = c(1, 1))


m1 <- -P_probit$Beta.mean[2]/P_probit$Beta.mean[3]
b1 <- -P_probit$Beta.mean[1]/P_probit$Beta.mean[3] #3?
a <- seq(40,-40,length.out=1000)
yd1 <- m1*a+b1

MLE_model$coefficients
P_probit$Beta.mean

m2 <- -P_original_L1$post.mean[2]/P_original_L1$post.mean[3]
b2 <- -P_original_L1$post.mean[1]/P_original_L1$post.mean[3] #3?

yd2 <- m2*a+b2

m3 <- -P_original_L2$post.mean[2]/P_original_L2$post.mean[3]
b3 <- -P_original_L2$post.mean[1]/P_original_L2$post.mean[3] #3?
yd3 <- m3*a+b3
plot(data$V2,data$V3,col=data$color,
     xlab = "X2",ylab = "X3",main = "original data distribution")
lines(a,yd1,type="l")
lines(a,yd2,type="l",col=2)
lines(a,yd3,type="l",col=3)



linedata <- data.frame(a=a,yd=yd,yd1=yd1,yd2=yd2)

linedata_new <- data.frame(a=rep(a,3),yd=c(yd,yd1,yd2),label=c(rep("coreset(L1)",length(a))
                                                               ,rep("original(L1)",length(a)),
                                                               rep("original(p=1.21)",length(a))))

jpeg("original data illustration", units="in", width=8, height=5, res=300)
figure0 <- ggplot(data) +
  geom_point(aes(x=V2, y=V3, color=factor(label)))+
  guides(color=guide_legend(title="label"))+
  xlab("x1")+ylab("x2")+ggtitle("Original simulation data")

  new_scale_color() +
  geom_line(data=linedata_new,aes(a, yd, color=label))+
  guides(color=guide_legend(title="coreset"))+
  scale_color_manual(values = c("original(L1)" = "black",
                                "original(p=1.21)" = "blue",
                                "coreset(L1)" = "red"))+

dev.off()
par(mfrow = c(1, 1))

####################1##########
model1 <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                          X=X[sample_index,], y=data$y[sample_index],initial_theda=MLE_model$coefficients,
                          true_theta=MLE_coefficients,Lp=1)
m <- -model1$post.mean[2]/model1$post.mean[3]
b <- -model1$post.mean[1]/model1$post.mean[3]
a <- seq(40,-40,length.out=length(sample_index))

yd <- m*a + b
yd1 <- m1*a+b1

m2 <- -P_original_L1$post.mean[2]/P_original_L1$post.mean[3]
b2 <- -P_original_L1$post.mean[1]/P_original_L1$post.mean[3]
a <- seq(40,-40,length.out=length(sample_index))
yd <- m*a + b
yd1 <- m1*a+b1
yd2 <- m2*a+b2

plot(data[sample_index,]$V2,data[sample_index,]$V3,col=data[sample_index,]$color,
     xlab = "X2",ylab = "X3",main = "compressed data from 1-probit coreset(fraction=0.01%)",
     ylim = c(min(c(yd,yd1,yd2,data$V3)),max(c(yd,yd1,yd2,data$V3))),
                                               xlim=c(min(data$V2),max(data$V2)))
                                             
lines(a,yd,type="l",col=3)
lines(a,yd1,type="l")
lines(a,yd2,type="l",col=4)
legend("bottomright",legend=c("original","coreset(l1)","original(l1)"), col=c(1,3,4), 
       bty = 'n', lwd = 1, lty = 1, cex = 0.8)

plot1 <- data[sample_index,]

linedata <- data.frame(a=a,yd=yd,yd1=yd1,yd2=yd2)

linedata_new <- data.frame(a=rep(a,3),yd=c(yd,yd1,yd2),label=c(rep("coreset(L1)",length(a))
                                                               ,rep("original(L1)",length(a)),
                                                               rep("original(p=1.21)",length(a))))

jpeg("L1 coreset illustration", units="in", width=8, height=5, res=300)
figure1 <- ggplot(plot1) +
  geom_point(aes(x=V2, y=V3, color=factor(label)))+
  guides(color=guide_legend(title="label"))+
  new_scale_color() +
  geom_line(data=linedata_new,aes(a, yd, color=label))+
  guides(color=guide_legend(title="coreset"))+
  scale_color_manual(values = c("original(L1)" = "black",
                                "original(p=1.21)" = "blue",
                                "coreset(L1)" = "red"))+
  xlab("x1")+ylab("x2")+ggtitle("L1 coreset")
dev.off()


Xb1 <- X%*%(model1$post.mean)
Xb2 <- X%*%(P_original_L1$post.mean)
plot(sort(Xb1),sort(pgnorm(Xb1,alpha = Phi_p(1),beta = 1)),type="l")
lines(sort(Xb2),sort(pgnorm(Xb2,alpha = Phi_p(1),beta = 1)),col=2)
plot(1:D,model1$post.mean,type = "l",ylim = c(-3,3),main = "comparison of coefficients",xlab = "Beta-index",
     ylab = "parameters")
lines(P_original_L1$post.mean,col=2)
#############################

model1 <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                          X=X[sample_index1,], y=data$y[sample_index1],initial_theda=MLE_model$coefficients,
                          true_theta=MLE_coefficients,Lp=2)
m <- -model1$post.mean[2]/model1$post.mean[3]
b <- -model1$post.mean[1]/model1$post.mean[3]
a <- seq(40,-40,length.out=length(sample_index1))
yd <- m*a + b
yd1 <- m1*a+b1

m2 <- -P_original_L2$post.mean[2]/P_original_L2$post.mean[3]
b2 <- -P_original_L2$post.mean[1]/P_original_L2$post.mean[3]
a <- seq(40,-40,length.out=length(sample_index1))
yd2 <- m2*a + b2
plot(data[sample_index1,]$V2,data[sample_index1,]$V3,col=data[sample_index1,]$color,
     xlab = "X2",ylab = "X3",main = "compressed data from 2-probit coreset(fraction=0.01%)",
     ylim = c(min(c(yd,yd1,yd2,data$V3)),max(c(yd,yd1,yd2,data$V3))),
     xlim=c(min(data$V2),max(data$V2)))
lines(a,yd,type="l",col=3)
lines(a,yd1,type="l")
lines(a,yd2,type="l",col=4)
legend("bottomright",legend=c("original","coreset(l2)","original(l2)"), col=c(1,3,4), 
       bty = 'n', lwd = 1, lty = 1, cex = 0.8)

plot2 <- data[sample_index1,]

linedata <- data.frame(a=a,yd=yd,yd1=yd1,yd2=yd2)

linedata_new <- data.frame(a=rep(a,3),yd=c(yd,yd1,yd2),label=c(rep("coreset(l2)",length(a))
                                                               ,rep("original(l2)",length(a)),
                                                               rep("original(p=1.21)",length(a))))

jpeg("L2 coreset illustration", units="in", width=8, height=5, res=300)

figure2 <- ggplot(plot2) +
  geom_point(aes(x=V2, y=V3, color=factor(label)))+
  guides(color=guide_legend(title="label"))+
  new_scale_color() +
  geom_line(data=linedata_new,aes(a, yd, color=label))+
  guides(color=guide_legend(title="coreset"))+
  scale_color_manual(values = c("original(l2)" = "black",
                                "original(p=1.21)" = "blue",
                                "coreset(l2)" = "red"))+
  xlab("x1")+ylab("x2")+ggtitle("L2 coreset")
dev.off()



Xb1 <- X%*%(model1$post.mean)
Xb2 <- X%*%(P_original_L2$post.mean)
plot(sort(Xb1),sort(pgnorm(Xb1,alpha = Phi_p(2),beta = 2)),type="l")
lines(sort(Xb2),sort(pgnorm(Xb2,alpha = Phi_p(2),beta = 2)),col=2)

plot(1:D,model1$post.mean,type = "l",ylim = c(-3,3),main = "comparison of coefficients",xlab = "Beta-index",
     ylab = "parameters")
lines(P_original_L2$post.mean,col=2)


#plot(data$V2,data$V3,col=data$color,
#     xlab = "X2",ylab = "X3",main = "original data distribution")
#m1 <- -MLE_coefficients[2]/MLE_coefficients[3]
#b1 <- -MLE_coefficients[1]/MLE_coefficients[2]
#a <- seq(12,-12,length.out=1000)
#yd1 <- m1*a+b1
#lines(a,yd1,type="l")

#########uniform########3
model1 <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                          X=X[sample_index2,], y=data$y[sample_index2],initial_theda=MLE_model$coefficients,
                          true_theta=MLE_coefficients,Lp=1)
m <- -model1$post.mean[2]/model1$post.mean[3]
b <- -model1$post.mean[1]/model1$post.mean[3]
m
m2
b2
a <- seq(40,-40,length.out=length(sample_index1))
yd <- m*a + b
yd1 <- m1*a+b1


model2 <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                          X=X[sample_index2,], y=data$y[sample_index2],initial_theda=MLE_model$coefficients,
                          true_theta=MLE_coefficients,Lp=2)
m2 <- -model2$post.mean[2]/model2$post.mean[3]
b2 <- -model2$post.mean[1]/model2$post.mean[3]
a <- seq(40,-40,length.out=length(sample_index1))
yd2 <- m2*a + b2
plot(data[sample_index2,]$V2,data[sample_index2,]$V3,col=data[sample_index2,]$color,
     xlab = "X2",ylab = "X3",main = "compressed data from uniform coreset(fraction=0.01%)",
     ylim = c(min(c(yd,yd1,yd2,data$V3)),max(c(yd,yd1,yd2,data$V3))),
     xlim=c(min(data$V2),max(data$V2)))

lines(a,yd,type="l",col=3)
lines(a,yd1,type="l")
lines(a,yd2,type="l",col=4)
legend("bottomright",legend=c("original","coreset(l1)","coreset(l2)"), col=c(1,3,4), 
       bty = 'n', lwd = 1, lty = 1, cex = 0.8)

plot3 <- data[sample_index2,]

linedata <- data.frame(a=a,yd=yd,yd1=yd1,yd2=yd2)

linedata_new <- data.frame(a=rep(a,3),yd=c(yd,yd1,yd2),label=c(rep("coreset(l1)",length(a))
                                                                   ,rep("original",length(a)),
                                                                   rep("coreset(l2)",length(a))))

jpeg("Uniform illustration", units="in", width=8, height=5, res=300)
figure3 <- ggplot(plot3) +
  geom_point(aes(x=V2, y=V3, color=factor(label)))+
  guides(color=guide_legend(title="label"))+
  new_scale_color() +
  geom_line(data=linedata_new,aes(a, yd, color=label))+
  guides(color=guide_legend(title="coreset"))+
  scale_color_manual(values = c("original" = "black",
                                "coreset(l1)" = "red",
                                "coreset(l2)" = "blue"))+
  xlab("x1")+ylab("x2")+ggtitle("Uniform coreset")
dev.off()
Xb1 <- X%*%(model1$post.mean)
Xb2 <- X%*%(model2$post.mean)
plot(sort(Xb1),sort(pgnorm(Xb1,alpha = Phi_p(1),beta = 1)),type="l")
lines(sort(Xb2),sort(pgnorm(Xb2,alpha = Phi_p(2),beta = 2)),col=2)

plot(1:D,model1$post.mean,type = "l",ylim = c(-3,3),main = "comparison of coefficients",xlab = "Beta-index",
     ylab = "parameters")
lines(model2$post.mean,col=2)

########one shot coreset#########
model1 <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                          X=X[sample_index4,], y=data$y[sample_index4],initial_theda=MLE_model$coefficients,
                          true_theta=MLE_coefficients,Lp=1)
m <- -model1$post.mean[2]/model1$post.mean[3]
b <- -model1$post.mean[1]/model1$post.mean[3]
m
m2
b2
a <- seq(40,-40,length.out=length(sample_index1))
yd <- m*a + b
yd1 <- m1*a+b1


model2 <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                          X=X[sample_index4,], y=data$y[sample_index4],initial_theda=MLE_model$coefficients,
                          true_theta=MLE_coefficients,Lp=2)
m2 <- -model2$post.mean[2]/model2$post.mean[3]
b2 <- -model2$post.mean[1]/model2$post.mean[3]
a <- seq(40,-40,length.out=length(sample_index1))
yd2 <- m2*a + b2

model3 <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                          X=X[sample_index4,], y=data$y[sample_index4],initial_theda=MLE_model$coefficients,
                          true_theta=MLE_coefficients,Lp=4)
m3 <- -model3$post.mean[2]/model3$post.mean[3]
b3 <- -model3$post.mean[1]/model3$post.mean[3]
a <- seq(40,-40,length.out=length(sample_index1))
yd3 <- m3*a + b3

plot(data[sample_index4,]$V2,data[sample_index4,]$V3,col=data[sample_index4,]$color,
     xlab = "X2",ylab = "X3",main = "compressed data from uniform coreset(fraction=0.01%)",
     ylim = c(min(c(yd,yd1,yd2,data$V3)),max(c(yd,yd1,yd2,data$V3))),
     xlim=c(min(data$V2),max(data$V2)))

lines(a,yd1,type="l")
lines(a,yd,type="l",col=3)
lines(a,yd2,type="l",col=4)
lines(a,yd3,type="l",col=6)

legend("bottomright",legend=c("original","coreset(l1)","coreset(l2)","coreset(l4)"), col=c(1,3,4,6), 
       bty = 'n', lwd = 1, lty = 1, cex = 0.8)

plot4 <- data[sample_index4,]

linedata <- data.frame(a=a,yd=yd,yd1=yd1,yd2=yd2,yd3=yd3)

linedata_new <- data.frame(a=rep(a,4),yd=c(yd,yd1,yd2,yd3),label=c(rep("coreset(l1)",length(a))
                                                                   ,rep("original",length(a)),
                                                                   rep("coreset(l2)",length(a)),
                                                                   rep("coreset(l4)",length(a))))


jpeg("one-shot illustration", units="in", width=8, height=5, res=300)
figure4 <- ggplot(plot4) +
  geom_point(aes(x=V2, y=V3, color=factor(label)))+
  guides(color=guide_legend(title="label"))+
  new_scale_color() +
  geom_line(data=linedata_new,aes(a, yd, color=label))+
  guides(color=guide_legend(title="coreset"))+
  scale_color_manual(values = c("original" = "black",
                                "coreset(l1)" = "red",
                                "coreset(l2)" = "green",
                                "coreset(l4)" = "orange"))+
  xlab("x1")+ylab("x2")+ggtitle("One-shot coreset")
figure4
dev.off()

####################rootl2##########
model1 <- Lp_gibbssampler(N_sim=N_sim,burn_in=burn_in,
                          X=X[sample_index5,], y=data$y[sample_index5],initial_theda=MLE_model$coefficients,
                          true_theta=MLE_coefficients,Lp=1)
m <- -model1$post.mean[2]/model1$post.mean[3]
b <- -model1$post.mean[1]/model1$post.mean[3]
a <- seq(40,-40,length.out=length(sample_index5))
yd <- m*a + b
yd1 <- m1*a+b1
P_original_L2

m2 <- -P_original_L1$post.mean[2]/P_original_L1$post.mean[3]
b2 <- -P_original_L1$post.mean[1]/P_original_L1$post.mean[3]

m2 <- -MLE_logit$coefficients[2]/MLE_logit$coefficients[3]
b2 <- -MLE_logit$coefficients[1]/MLE_logit$coefficients[3]
a <- seq(40,-40,length.out=length(sample_index5))
yd <- m*a + b
yd1 <- m1*a+b1
yd2 <- m2*a+b2

plot(data[sample_index5,]$V2,data[sample_index5,]$V3,col=data[sample_index5,]$color,
     xlab = "X2",ylab = "X3",main = "compressed data from 1-probit coreset(fraction=0.01%)",
     ylim = c(min(c(yd,yd1,yd2,data$V3)),max(c(yd,yd1,yd2,data$V3))),
     xlim=c(min(data$V2),max(data$V2)))

lines(a,yd,type="l",col=3)
lines(a,yd1,type="l")
lines(a,yd2,type="l",col=4)
legend("bottomright",legend=c("original(p=1.21)","coreset(rootl2)","original(1-probit)"), col=c(1,3,4), 
       bty = 'n', lwd = 1, lty = 1, cex = 0.8)


plot5 <- data[sample_index5,]
linedata <- data.frame(a=a,yd=yd,yd1=yd1)
linedata_new <- data.frame(a=rep(a,2),yd=c(yd,yd1),label=c(rep("coreset(rootl2)",length(a))
                                                               ,rep("original(1-probit)",length(a))))

jpeg("rootl2 coreset illustration", units="in", width=8, height=5, res=300)
figure5 <- ggplot(plot5) +
  geom_point(aes(x=V2, y=V3, color=factor(label)))+
  guides(color=guide_legend(title="label"))+
  new_scale_color() +
  geom_line(data=linedata_new,aes(a, yd, color=label))+
  guides(color=guide_legend(title="coreset"))+
  scale_color_manual(values = c("original(1-probit)" = "black",
                                "coreset(rootl2)" = "red"))+
  xlab("x1")+ylab("x2")+ggtitle("Rootl2 coreset")

dev.off()
jpeg(" illustration summary", units="in", width=10, height=5, res=300)


plot_grid(figure0, figure1, figure2, figure3,figure4,figure5,
          labels = c("A", "B", "C","D","E","F"),
          ncol = 3, nrow = 2)
dev.off()
  # scale_color_manual(values = c("survival1"="dodgerblue3",
  #                               "survival10"="dodgerblue2",  
  #                               "survival50"="steelblue", 
  #                               "survival90"="dodgerblue1", 
  #                               "survival99"="cornflowerblue", 
  #                               "lit.data"="orange",
  #                               "Lung"="darkslategrey",
  #                               "shock.data"="violetred4")) 

## add the new scale here


  ## then add a new color aesthetic, all else as per usual
# 
#   scale_color_manual(NULL, values = "grey", guide = guide_legend(order = 2))
# 
#   geom_line(data=linedata,aes(a, yd1),color = "blue")+
#   geom_line(data=linedata,aes(a, yd2),color = "red")+
#   geom_line(data=linedata,aes(a, yd3),color = "green")



Xb1 <- X%*%(model1$post.mean)
Xb2 <- X%*%(model2$post.mean)
plot(sort(Xb1),sort(pgnorm(Xb1,alpha = Phi_p(1),beta = 1)),type="l")
lines(sort(Xb2),sort(pgnorm(Xb2,alpha = Phi_p(2),beta = 2)),col=2)

plot(1:D,model1$post.mean,type = "l",ylim = c(-3,3),main = "comparison of coefficients",xlab = "Beta-index",
     ylab = "parameters")
lines(model2$post.mean,col=2)


