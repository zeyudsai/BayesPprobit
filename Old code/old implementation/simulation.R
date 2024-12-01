library(BayesPprobit)

library(mvtnorm)
library(Matrix)
library(ggplot2)
# 基础模拟数据生成函数

generate.data <- function(N,D){
  #x is generated from a multivariate normal distribution
  b <- D-1
  sigma <- matrix(rep(1,b*b),nrow = b)
  for (i in 1:b) {
    for (j in 1:b) {
      sigma[i,j] <- 2*(0.5^abs(i-j))
    }
  }

  #a <- rmvnorm(n=N,mean=c(rep(-5,(D-1)/5),rep(5,(D-1)/5),
  #                       rep(0,(D-1)/5)),sigma=sigma)

  a <- mvtnorm::rmvnorm(n=N,
                        mean=c(rep(-2,(D-1)/5),rep(2,(D-1)/5),
                          rep(-3,(D-1)/5),rep(3,(D-1)/5),
                          rep(0,(D-1)/5)),sigma=sigma)
  # Create n x D design matrix
  #X <- matrix(c(rep(1, N), a,b), ncol = D)
  X <- matrix(c(rep(1, N), a), ncol = D)/10

  return(X)
}

# 运行模拟
N <- 50000
#D should always be odd number since we use (D-1)/2 to generate the variables
D <- 11
N_sim <- 10000
burn_in<- 2000
#gengerate data and true parameters
X <- generate.data(N,D)
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
# 拟合初始MLE模型
MLE_probit <- glm(y.probit~X[,-1],family = binomial(link = probit))
MLE_logit <- glm(y.logit~X[,-1],
                 family = binomial(link = logit))

summary(MLE_logit)
summary(MLE_probit)
true_theta
# 存储所有结果
results <- list()

# 处理特殊情况: logit和probit
special_cases <- list(
  logit = list(y = y.logit,
               init = model_logit$coefficients,
               p_range = c(0.5, 3)),
  probit = list(y = y.probit,
                init = model_probit$coefficients,
                p_range = c(0.5, 5))
)

# 运行特殊情况
for(case in names(special_cases)) {
  message("\nRunning ", case, " model...")
  results[[case]] <- multi_chain(
    n_sim = 1000,
    burn_in = 10,
    X = X,
    y = special_cases[[case]]$y,
    initial_theta = special_cases[[case]]$init,
    true_theta = true_theta,
    initial_p = 2,
    mh_iter = 200,
    p_range = special_cases[[case]]$p_range,
    step_size = 0.05,
    n_chains = 2
  )
  print(results[[case]])
}

results$logit
results$probit
# 使用示例
print(all_results)
dim(X)
dim(X_rest)
all_results$probit
true_theta
all_results$logit

summary(model_logit)


p_values <- c(0.5, 1, 1.5, 2, 3, 4, 5, 8)
# 4. 生成不同p值下的响应变量
y_list <- list()
pi_list <- list()
for(p in p_values) {
  Xb <- X %*% true_theta + gnorm::pgnorm(N, mu = 0, alpha = p_scale(p),
                                         beta = p)


  # 计算概率
  pi_p <- gnorm::pgnorm(Xb, mu = 0, alpha = p_scale(p),
                        beta = p)
  pi_list[[paste0("pi_", p)]] <- pi_p

  # 生成二项分布响应
  y_p <- rbinom(N, 1, pi_p)
  y_list[[paste0("y_p", p)]] <- y_p
}


for(p in p_values) {
  message("\nRunning p = ", p, " model...")
  results[[paste0("p", p)]] <- multi_chain(
    n_sim = N_sim,
    burn_in = burn_in,
    X = X,
    y = y_list[[paste0("y_p", p)]],
    initial_theta = model_probit$coefficients,
    true_theta = true_theta,
    initial_p = 2,
    mh_iter = 200,
    p_range = c(0.5, 5),
    step_size = 0.05,
    n_chains = 2
  )
  print(results[[paste0("p", p)]])
}


