

###Plotting for illustration figures
plot(Xb, dgnorm(Xb, mu = 0, alpha = sqrt(8), beta = 8),
     xlab = "x", ylab = expression(p(x)))
Xb <- subset(Xb,Xb<10)
Xb <- subset(Xb,-10<Xb)

density0.5 <- data.frame(Xb=Xb,density=dgnorm(Xb,alpha = p_scale(0.5),beta = 0.5))
density1 <- data.frame(Xb=Xb,density=dgnorm(Xb,alpha = p_scale(1),beta = 1))
density1.5 <- data.frame(Xb=Xb,density=dgnorm(Xb,alpha = p_scale(1.5),beta = 1.5))
density2 <- data.frame(Xb=Xb,density=dgnorm(Xb,alpha = p_scale(2),beta = 2))
density2.5 <- data.frame(Xb=Xb,density=dgnorm(Xb,alpha = p_scale(2.5),beta = 2.5))
density3 <- data.frame(Xb=Xb,density=dgnorm(Xb,alpha = p_scale(3),beta = 3))
density3.5 <- data.frame(Xb=Xb,density=dgnorm(Xb,alpha = p_scale(3.5),beta = 3.5))
density4 <- data.frame(Xb=Xb,density=dgnorm(Xb,alpha = p_scale(4),beta = 4))
density5 <- data.frame(Xb=Xb,density=dgnorm(Xb,alpha = p_scale(5),beta = 5))
density8<- data.frame(Xb=Xb,density=dgnorm(Xb,alpha = p_scale(8),beta = 8))
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

prob0.5 <- data.frame(Xb=Xb,pi=pgnorm(Xb,alpha = p_scale(0.5),beta = 0.5))
prob1 <- data.frame(Xb=Xb,pi=pgnorm(Xb,alpha = p_scale(1),beta = 1))
prob1.5 <- data.frame(Xb=Xb,pi=pgnorm(Xb,alpha = p_scale(1.5),beta = 1.5))
prob2 <- data.frame(Xb=Xb,pi=pgnorm(Xb,alpha = p_scale(2),beta = 2))
prob2.5 <- data.frame(Xb=Xb,pi=pgnorm(Xb,alpha = p_scale(2.5),beta = 2.5))
prob3 <- data.frame(Xb=Xb,pi=pgnorm(Xb,alpha = p_scale(3),beta = 3))
prob3.5 <- data.frame(Xb=Xb,pi=pgnorm(Xb,alpha = p_scale(3.5),beta = 3.5))
prob4 <- data.frame(Xb=Xb,pi=pgnorm(Xb,alpha = p_scale(4),beta = 4))
prob5 <- data.frame(Xb=Xb,pi=pgnorm(Xb,alpha = p_scale(5),beta = 5))
prob8<- data.frame(Xb=Xb,pi=pgnorm(Xb,alpha = p_scale(8),beta = 8))
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


df <- data.frame(Xb=Xb,aaa=dgnorm(Xb,alpha = p_scale(2),beta = 2))
df2 <- data.frame(Xb=Xb,aaa=pgnorm(Xb,alpha = p_scale(2),beta = 2))


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

