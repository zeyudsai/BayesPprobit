# Set seed for reproducible results'
library("expint")
library("pgnorm")
library("gnorm")
library(ggplot2)



set.seed(1234)
par(mfrow = c(1,1))
dens <- function(X, a, b, mu){
  Y <- abs(X-mu)/a
  Y <- exp(-Y^b)
  Y <- Y*b/(2*a)
  Y <-Y/ exp(log(gamma(1/b)))
  return(Y)
}
cdf<- function(X, a, b, mu){
  Y <- abs((X-mu)/a)^b
  Y <- gammainc(1/b,Y)
  Y <- Y * sign(X-mu) / 2 + 0.5
  return(Y)
}
cdf(1,1,1,1)

w <- 1.5

X <- seq(-3, 3, 0.01)
plot(X,dens(X,1,0.5,0), type = "l",ylim = c(0,1),col=1)
lines(X,dens(X,1,1,0),col=2)
lines(X,dens(X,1,1.5,0),col=3)
lines(X,dens(X,1,2,0))
lines(X,dens(X,1,3,0))
lines(X,dens(X,1,8,0))
# Plot generalized normal/exponential power density
# that corresponds to the standard normal density
phi_p <- function(p){
  result <- (1/p)^(1/p)
  return(result)
}
Phi_p <- function(p){
  result <- p^(1/p)
  return(result)
}
phi_p(0.5)
phi_p(1)
phi_p(2)
Phi_p(2)
sqrt(2)
1/sqrt(2)

truncated_gnorm <- function(n, range, mu, sd,alpha,beta) {
  
  # range is a vector of two values
  
  F.a <- pgnorm(min(range), alpha=alpha,beta=beta, mu = mu)
  F.b <- pgnorm(max(range), alpha=alpha,beta=beta, mu = mu)
  
  u <- runif(n, min = F.a, max = F.b)
  
  qgnorm(u, mu = mu,alpha = alpha,beta=beta)
  
}
truncated_gnorm(10,c(-1,0),mu=0,sd=1,alpha = 2,beta = 2)
truncated_gnorm(10,c(0,1),mu=0,sd=1,alpha = 2,beta = 2)

xs <- seq(-3, 3, length.out = 1000)

plot(Xb, dgnorm(Xb, mu = 0, alpha = sqrt(2), beta = 2), type = "l",
     xlab = "x", ylab = expression(p(x)))
plot(xs, dgnorm(xs, mu = 0, alpha = (1/2)^(1/2), beta = 2), type = "l",
     xlab = "x", ylab = expression(p(x)))


n <- length(xs)
plotdata <- data.frame(xs=rep(xs,8),
                       probability=c(dgnorm(xs, mu = 0, alpha = Phi_p(0.5), beta = 0.5),
                                     dgnorm(xs, mu = 0, alpha = Phi_p(1), beta = 1),
                                     dgnorm(xs, mu = 0, alpha = Phi_p(2), beta = 2),
                                     dgnorm(xs, mu = 0, alpha = Phi_p(3), beta = 3),
                                     dgnorm(xs, mu = 0, alpha = Phi_p(5), beta = 5),
                                     dgnorm(xs, mu = 0, alpha = Phi_p(6), beta = 6),
                                     dgnorm(xs, mu = 0, alpha = Phi_p(8), beta = 8),
                                     dgnorm(xs, mu = 0, alpha = Phi_p(100), beta = 100)),
                       cdf=c(pgnorm(xs, mu = 0, alpha = Phi_p(0.5), beta = 0.5),
                             pgnorm(xs, mu = 0, alpha = Phi_p(1), beta = 1),
                             pgnorm(xs, mu = 0, alpha = Phi_p(2), beta = 2),
                             pgnorm(xs, mu = 0, alpha = Phi_p(3), beta = 3),
                             pgnorm(xs, mu = 0, alpha = Phi_p(5), beta = 5),
                             pgnorm(xs, mu = 0, alpha = Phi_p(6), beta = 6),
                             pgnorm(xs, mu = 0, alpha = Phi_p(8), beta = 8),
                             pgnorm(xs, mu = 0, alpha = Phi_p(100), beta = 100)),
                       
                       p=c(rep(0.5,n),
                           rep(1,n),
                           rep(2,n),
                           rep(3,n),
                           rep(5,n),
                           rep(6,n),
                           rep(8,n),
                           rep(100,n)))

ggplot(plotdata,aes(x=xs,y=probability,color=factor(p))) +
  geom_line()+
  guides(color=guide_legend(title="p"))+
  # ggtitle("Probability density function",
          # subtitle = "p-generalized normal distribution"
          # )+
  theme_bw()+
scale_color_manual(values = c("0.5" = brewer.pal(8, "Dark2")[1],
                              "1" = brewer.pal(8, "Dark2")[2],
                              "2"=brewer.pal(8, "Dark2")[3],
                              "3"=brewer.pal(8, "Dark2")[4],
                              "5"=brewer.pal(8, "Dark2")[6],
                              "6"=brewer.pal(8, "Dark2")[7],
                              "8"=brewer.pal(8, "Dark2")[8],
                              "100"=brewer.pal(8, "Set1")[2]))+
  xlab("x")
dev.off()




ggplot(plotdata,aes(x=xs,y=cdf,color=factor(p))) +
  geom_line()+
  guides(color=guide_legend(title="p"))+
  # ggtitle("Cumulative density function",
          # subtitle = "p-generalized normal distribution"
          # )+
  theme_bw()+
  scale_color_manual(values = c("0.5" = brewer.pal(8, "Dark2")[1],
                                "1" = brewer.pal(8, "Dark2")[2],
                                "2"=brewer.pal(8, "Dark2")[3],
                                "3"=brewer.pal(8, "Dark2")[4],
                                "5"=brewer.pal(8, "Dark2")[6],
                                "6"=brewer.pal(8, "Dark2")[7],
                                "8"=brewer.pal(8, "Dark2")[8],
                                "100"=brewer.pal(8, "Set1")[2]))+
  ylab("cumulative probability")+xlab("x")
dev.off()

  