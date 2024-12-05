## -----------------------------------------------------------------------------
set.seed(88) 
n <- 1000
sigma <- 1
u <- runif(n)  # 生成 n 个 [0, 1) 之间的均匀随机数
rayleigh_n <- sigma * sqrt(-2 * log(u))  # Rayleigh 变换
hist(rayleigh_n, prob = TRUE, ylim = c(0, 0.8), main = expression(paste("Rayleigh Distribution ",sigma==1)), xlab = "x",  col = "lightblue")
# 绘制原本的图像
y <- seq(0,4,0.01)
yy <- (y/sigma^2)*exp(-y^2/(2*sigma^2))
lines(y,yy)

## -----------------------------------------------------------------------------
set.seed(88)
n <- 1000
sigma <- 2
u <- runif(n)  # 生成 n 个 [0, 1) 之间的均匀随机数
rayleigh_n <- sigma * sqrt(-2 * log(u))  # Rayleigh 变换
hist(rayleigh_n, prob = TRUE, ylim = c(0, 0.4), main = expression(paste("Rayleigh Distribution ",sigma==2)), xlab = "x",  col = "lightblue")
# 绘制原本的图像
y <- seq(0,8,0.01)
yy <- (y/sigma^2)*exp(-y^2/(2*sigma^2))
lines(y,yy)

## -----------------------------------------------------------------------------
set.seed(88)
n <- 1000
sigma <- 3
u <- runif(n)  # 生成 n 个 [0, 1) 之间的均匀随机数
rayleigh_n <- sigma * sqrt(-2 * log(u))  # Rayleigh 变换
hist(rayleigh_n, prob = TRUE, ylim = c(0, 0.3), main = expression(paste("Rayleigh Distribution ",sigma==3)), xlab = "x",  col = "lightblue")
# 绘制原本的图像
y <- seq(0,12,0.01)
yy <- (y/sigma^2)*exp(-y^2/(2*sigma^2))
lines(y,yy)

## -----------------------------------------------------------------------------
set.seed(123)
n <- 1000
u <- runif(n)
k <- as.integer(u>0.25)
x1 <- rnorm(1000,0,1)
x2 <- rnorm(1000,3,1)
x <- k*x1+k*x2
hist(x,prob=TRUE,col = "lightblue")

## ----echo=FALSE---------------------------------------------------------------
par(mfrow = c(2, 3))

set.seed(123)
n <- 1000
u <- runif(n)
k10 <- as.integer(u>0)
k11 <- as.integer(u>0.1)
k12 <- as.integer(u>0.2)
k13 <- as.integer(u>0.3)
k14 <- as.integer(u>0.4)
k15 <- as.integer(u>0.5)
x1 <- rnorm(1000,0,1)
x2 <- rnorm(1000,3,1)
x10 <- k10*x1+k10*x2
x11 <- k11*x1+k11*x2
x12 <- k12*x1+k12*x2
x13 <- k13*x1+k13*x2
x14 <- k14*x1+k14*x2
x15 <- k15*x1+k15*x2
hist(x10,prob=TRUE,col = "lightblue", xlab='x', main = expression(paste("Distribution ",p1==0)))
hist(x11,prob=TRUE,col = "lightblue", xlab='x', main = expression(paste("Distribution ",p1==0.1)))
hist(x12,prob=TRUE,col = "lightblue", xlab='x', main = expression(paste("Distribution ",p1==0.2)))
hist(x13,prob=TRUE,col = "lightblue", xlab='x', main = expression(paste("Distribution ",p1==0.3)))
hist(x14,prob=TRUE,col = "lightblue", xlab='x', main = expression(paste("Distribution ",p1==0.4)))
hist(x15,prob=TRUE,col = "lightblue", xlab='x', main = expression(paste("Distribution ",p1==0.5)))

par(mfrow = c(1, 1))

## ----echo=FALSE---------------------------------------------------------------
par(mfrow = c(2, 3))

set.seed(123)
n <- 1000
u <- runif(n)
k10 <- as.integer(u>0)
k11 <- as.integer(u>0.02)
k12 <- as.integer(u>0.04)
k13 <- as.integer(u>0.06)
k14 <- as.integer(u>0.08)
k15 <- as.integer(u>0.05)
x1 <- rnorm(1000,0,1)
x2 <- rnorm(1000,3,1)
x10 <- k10*x1+k10*x2
x11 <- k11*x1+k11*x2
x12 <- k12*x1+k12*x2
x13 <- k13*x1+k13*x2
x14 <- k14*x1+k14*x2
x15 <- k15*x1+k15*x2
hist(x10,prob=TRUE,col = "lightblue", xlab='x', main = expression(paste("Distribution ",p1==0)))
hist(x11,prob=TRUE,col = "lightblue", xlab='x', main = expression(paste("Distribution ",p1==0.02)))
hist(x12,prob=TRUE,col = "lightblue", xlab='x', main = expression(paste("Distribution ",p1==0.04)))
hist(x15,prob=TRUE,col = "lightblue", xlab='x', main = expression(paste("Distribution ",p1==0.05)))
hist(x13,prob=TRUE,col = "lightblue", xlab='x', main = expression(paste("Distribution ",p1==0.06)))
hist(x14,prob=TRUE,col = "lightblue", xlab='x', main = expression(paste("Distribution ",p1==0.08)))

par(mfrow = c(1, 1))

## -----------------------------------------------------------------------------
library(ggplot2)
set.seed(123)

# 设置参数
lambda <- 1          # Poisson 过程的平均事件数
alpha <- 2           # Gamma 分布的形状参数
beta <- 1            # Gamma 分布的尺度参数
T <- 10              # 模拟时间长度

simulate_compound_poisson_gamma <- function(lambda, alpha, beta, t) {
  # 生成 Poisson 随机变量 N
  N <- rpois(1, lambda * T)  # 在时间区间 T 内生成的事件数量
  # 生成 Gamma 随机变量并求和
  gamma_samples <- rgamma(N, shape = alpha, scale = beta)  # 生成 N 个 Gamma 样本
  total_amount <- sum(gamma_samples)  # 计算总和
  return(total_amount)
}
  
# 进行多次模拟
num_simulations <- 10000
results <- replicate(num_simulations, simulate_compound_poisson_gamma(lambda, alpha, beta, t))

# 计算均值和方差
mean_X <- mean(results)
var_X <- var(results)

# 输出结果
cat("X(10) 的均值:", mean_X, "\n")
cat("X(10) 的方差:", var_X, "\n")
cat("均值理论值为：",lambda*T*alpha*beta,"真实值与理论值的差距为：",abs(lambda*T*alpha*beta-mean_X),"\n")
cat("方差理论值为：",lambda*T*(alpha*beta^2+(alpha*beta)^2),"真实值与理论值的差距为：",abs(lambda*T*(alpha*beta^2+(alpha*beta)^2)-var_X), "\n")


## -----------------------------------------------------------------------------
library(ggplot2)
set.seed(123)

# 设置参数
lambda <- 5          # Poisson 过程的平均事件数
alpha <- 2           # Gamma 分布的形状参数
beta <- 1            # Gamma 分布的尺度参数
T <- 10              # 模拟时间长度

simulate_compound_poisson_gamma <- function(lambda, alpha, beta, t) {
  # 生成 Poisson 随机变量 N
  N <- rpois(1, lambda * T)  # 在时间区间 T 内生成的事件数量
  # 生成 Gamma 随机变量并求和
  gamma_samples <- rgamma(N, shape = alpha, scale = beta)  # 生成 N 个 Gamma 样本
  total_amount <- sum(gamma_samples)  # 计算总和
  return(total_amount)
}
  
# 进行多次模拟
num_simulations <- 10000
results <- replicate(num_simulations, simulate_compound_poisson_gamma(lambda, alpha, beta, t))

# 计算均值和方差
mean_X <- mean(results)
var_X <- var(results)

# 输出结果
cat("X(10) 的均值:", mean_X, "\n")
cat("X(10) 的方差:", var_X, "\n")
cat("均值理论值为：",lambda*T*alpha*beta,"真实值与理论值的差距为：",abs(lambda*T*alpha*beta-mean_X), "\n")
cat("方差理论值为：",lambda*T*(alpha*beta^2+(alpha*beta)^2),"真实值与理论值的差距为：",abs(lambda*T*(alpha*beta^2+(alpha*beta)^2)-var_X), "\n")


## -----------------------------------------------------------------------------
library(ggplot2)
set.seed(123)

# 设置参数
lambda <- 10          # Poisson 过程的平均事件数
alpha <- 2           # Gamma 分布的形状参数
beta <- 1            # Gamma 分布的尺度参数
T <- 10              # 模拟时间长度

simulate_compound_poisson_gamma <- function(lambda, alpha, beta, t) {
  # 生成 Poisson 随机变量 N
  N <- rpois(1, lambda * T)  # 在时间区间 T 内生成的事件数量
  # 生成 Gamma 随机变量并求和
  gamma_samples <- rgamma(N, shape = alpha, scale = beta)  # 生成 N 个 Gamma 样本
  total_amount <- sum(gamma_samples)  # 计算总和
  return(total_amount)
}
  
# 进行多次模拟
num_simulations <- 10000
results <- replicate(num_simulations, simulate_compound_poisson_gamma(lambda, alpha, beta, t))

# 计算均值和方差
mean_X <- mean(results)
var_X <- var(results)

# 输出结果
cat("X(10) 的均值:", mean_X, "\n")
cat("X(10) 的方差:", var_X, "\n")
cat("均值理论值为：",lambda*T*alpha*beta,"真实值与理论值的差距为：",abs(lambda*T*alpha*beta-mean_X), "\n")
cat("方差理论值为：",lambda*T*(alpha*beta^2+(alpha*beta)^2),"真实值与理论值的差距为：",abs(lambda*T*(alpha*beta^2+(alpha*beta)^2)-var_X), "\n")


## -----------------------------------------------------------------------------
# 设置参数
set.seed(12345)  # 固定随机数种子
m <- 10000  # 蒙特卡罗模拟次数
y <- runif(m)

# Beta分布的概率密度函数
beta_cdf <- function(x){
  return (mean(30 * x^3 * y^2 * (1-x * y)^2))
}

## -----------------------------------------------------------------------------
x <- seq(0.1,1,length=10)
cdf_beta <- numeric(length(t))
for (i in 1:length(x)){
  # 估计CDF
  cdf_beta[i] <- beta_cdf(x[i])
}
# 标准
Phi <- pbeta(x,3,3)
round(rbind(x, cdf_beta, Phi), 3)

## -----------------------------------------------------------------------------
# 设置参数
set.seed(123)  # 固定随机数种子
m <- 10000  # 蒙特卡罗模拟次数

# 模拟函数 Rayleigh分布的概率密度函数
Rayleigh_cdf <- function(x, sigma, antithetic = TRUE) {
  u <- runif(m/2)
  if (!antithetic) v <- runif(m/2) else v<- (1 - u)
  u <- c(u, v)
  cdf_ray <- numeric(length(x))
  for (i in 1:length(x)) {
    cdf_ray[i] <- (mean((x[i]^2 * u / sigma^2) *exp(-(x[i]*u)^2/(2* sigma^2))))
  }
  cdf_ray
}

## -----------------------------------------------------------------------------
# 标准的瑞利分布
sigma <- 2
x <- seq(0.1,2,length=20)

Phi_ray <- (1-exp(-x^2 / (2*sigma^2)))
cdf_ray1<- Rayleigh_cdf(x,sigma,antithetic=FALSE)
cdf_ray2<- Rayleigh_cdf(x,sigma)

round(rbind(x, cdf_ray1, cdf_ray2, Phi_ray), 4)

## -----------------------------------------------------------------------------
m <- 1000
MC1 <- MC2 <- numeric(m)
x <- 1
for (i in 1:m) {
MC1[i] <- Rayleigh_cdf(x, sigma, anti = FALSE)
MC2[i] <- Rayleigh_cdf(x, sigma)
}
print(sd(MC1))
print(sd(MC2))
print((var(MC1) - var(MC2))/var(MC1))

## -----------------------------------------------------------------------------
m <- 10000 #模拟次数
set.seed(123)

theta.hat <- se <- numeric(2)

# 原函数
g <- function(x) {
  exp(-x^2/2 + log((x^2)/sqrt(2*pi))) * (x > 1)
}

# 重要性函数1
x1 <- rexp(m, 1) #生成指数分布的随机数
i <- c(which(x1 < 1))
x1[i] <- 0
fg1 <- g(x1) / exp(-x1)
theta.hat[1] <- mean(fg1)
se[1] <- sd(fg1)

# 重要性函数2
x2 <- rnorm(m,1,1) #生成正态分布的随机数
i <- c(which(x2 < 1))
x2[i] <- 0
fg2 <- g(x2) / dnorm(x2,1,1)
theta.hat[2] <- mean(fg2)
se[2] <- sd(fg2)

rbind(theta.hat, se)

## -----------------------------------------------------------------------------
k<-seq(0,5,0.01)
kk<-g(k)
plot(k,kk)
plot(k,exp(-k)*(k>1))
kkk<-1/sqrt(2*pi)*exp(-(k-1)^2/2)*(k>1)
plot(k,kkk)

## -----------------------------------------------------------------------------
n = c(1e4,2*1e4,4*1e4,6*1e4,8*1e4)
quick_sort<-function(x){
  num<-length(x)
  if(num==0||num==1){return(x)
  }else{
    a<-x[1]
    y<-x[-1]
    lower<-y[y<a]
    upper<-y[y>=a]
    return(c(quick_sort(lower),a,quick_sort(upper)))}
}
time = numeric(length(n))
for (i in 1:length(n)){
  test<-sample(1:(n[i]))
  sort_result = quick_sort(test)
  time[i]=system.time(quick_sort(test))[1]
}
# system.time(quick_sort(sample(1:(2*1e4))))[1]

## ----eval = FALSE-------------------------------------------------------------
#  # m是模拟次数
#  # time_vc存储每一次模拟的时间，是一个1*m的向量
#  # time_mean存储平均时间，是一个1*5的向量
#  m=100
#  time_vc = numeric(m)
#  time_mean = numeric(length(n))
#  for (i in 1:length(n)){
#    for (j in 1:m){
#      test<-sample(1:(n[i]))
#      time_vc[j] = system.time(quick_sort(test))[1]
#    }
#    time_mean[i] = mean(time_vc)
#  }
#  print(time_mean)

## ----eval = FALSE-------------------------------------------------------------
#  nn = seq(1,4*1e4,5*1e3)
#  time_mean_ = numeric(length(nn))
#  for (i in 1:length(nn)){
#    for (j in 1:m){
#      test<-sample(1:(nn[i]))
#      time_vc[j] = system.time(quick_sort(test))[1]
#    }
#    time_mean_[i] = mean(time_vc)
#  }
#  print(time_mean_)

## -----------------------------------------------------------------------------
# 运行时间太长了这里记录运行结果
nn = seq(1,4*1e4,5*1e3)
time_mean_ = c(0.0000,0.0085,0.0122,0.0260,0.0320,0.0404,0.0431,0.0467)
model <- lm(time_mean_ ~ nn) 
plot(nn,time_mean_,main = "散点图及回归线", xlab = "nlog(n)", ylab = "平均运行时间",pch = 19, col = "blue")
abline(model, col = "red", lwd = 2)  # 使用回归模型绘制回归线

## -----------------------------------------------------------------------------
get_data <- function(n,mu,sigma){rnorm(n,mu,sigma)}
calculate_data <- function(sample.x){
  mu <- mean(sample.x)
  sigma <- sd(sample.x)
  SK <- 1/n/sigma^3*sum((sample.x-mu)^3)
  return(SK)
}
show_data <- function(SK.li, q.li, m, n){
  hat.SK.q <- quantile(SK.li,q.li)
  sd.Sk.q <- sqrt(q.li*(1-q.li)/m/dnorm(hat.SK.q)^2)
  SK.q_large <- sqrt(6/n)*qnorm(q.li)
  return (data.frame('large.sample'=SK.q_large,'estimate'=hat.SK.q,'estimate.sd'=sd.Sk.q))
}
m <- 1000
n <- 50000
mu <- 0
sigma <- 1
hat.SK <- numeric(m)
q.li <- c(0.025, 0.05, 0.95, 0.975)
for (j in 1:m){
  sample.j <- get_data(n,mu,sigma)
  hat.SK[j] <- calculate_data(sample.j)
}
show_data(hat.SK,q.li,m,n)

## -----------------------------------------------------------------------------
get_data <- function(n){
  x <- rnorm(n)
  y <- rnorm(n)
  return(list(x=x,y=y))
}
calculate_data <- function(sample.x,sample.y){
  x <- sample.x
  y <- sample.y
  p1 <- cor.test(x,y,method='pearson')$p.value
  p2 <- cor.test(x,y,method='kendall')$p.value
  p3 <- cor.test(x,y,method='spearman')$p.value
  return(list(p1=p1,p2=p2,p3=p3))
}
show_data <- function(p1.li,p2.li,p3.li, alpha){
  c(mean(p1.li<=alpha),mean(p2.li<=alpha),mean(p3.li<=alpha))
}
m <- 5e4
n <- 10
alpha <- 0.05
p1 <- p2 <- p3 <- numeric(m)
for (j in 1:m){
  sample.j <- get_data(n)
  x <- sample.j$x
  y <- sample.j$y
  p.j <- calculate_data(x,y)
  p1[j] <- p.j$p1
  p2[j] <- p.j$p2
  p3[j] <- p.j$p3
}
show_data(p1,p2,p3,alpha)

## -----------------------------------------------------------------------------
get_data <- function(n){
    y <- rexp(n)
    x <- 1 + 1.5*y^{0.5} + rt(n,2)
    return(list(x=x,y=y))
 }
 m <- 5e3
 n <- 100
 alpha <- 0.05
 p1 <- p2 <- p3 <- numeric(m)
 for (j in 1:m){
    sample.j <- get_data(n)
    x <- sample.j$x
    y <- sample.j$y
    p.j <- calculate_data(x,y)
    p1[j] <- p.j$p1
    p2[j] <- p.j$p2
    p3[j] <- p.j$p3
 }
 show_data(p1,p2,p3,alpha)

## -----------------------------------------------------------------------------
n <- 10000
x.bar <- 0.651
y.bar <- 0.676
p.hat <- 0.6635
z <- (x.bar-y.bar)/sqrt(2*p.hat*(1-p.hat)/n)
p.value <- 2*min(pnorm(z),1-pnorm(z))
p.value

## -----------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
set.seed(123)  # 设置随机数种子，保证结果可重复
# 参数设置
N <- 1000       
# 总假设数
n_null <- 950   # 原假设数
n_alt <- 50     
m <- 10000      
alpha <- 0.1    
# 初始化结果矩阵
# 备择假设数
# 模拟次数
# 显著性水平
results <- matrix(0, nrow = 3, ncol = 2)
rownames(results) <- c("FWER", "FDR", "TPR")
colnames(results) <- c("Bonferroni校正", "B-H校正")
 
 # 加载必要的库
library(ggplot2)
 # 初始化存储变量
fwer_bonf <- numeric(m)
fdr_bonf <- numeric(m)
tpr_bonf <- numeric(m)
fwer_bh <- numeric(m)
fdr_bh <- numeric(m)
tpr_bh <- numeric(m)
 # 执行模拟
for (i in 1:m) {
  # 生成 p 值
  p_null <- runif(n_null)               # 原假设的 p 值
  p_alt <- rbeta(n_alt, 0.1, 1)         # 备择假设的 p 值
  p_values <- c(p_null, p_alt)          # 合并所有 p 值
  true_labels <- c(rep(0, n_null), rep(1, n_alt))  # 标签：0-原假设，1-备择假设
  
  # Bonferroni 校正
  p_bonf <- p.adjust(p_values, method = "bonferroni")
  rejected_bonf <- which(p_bonf < alpha)
  
  V_bonf <- sum(true_labels[rejected_bonf] == 0)  # 误拒的原假设数
  S_bonf <- sum(true_labels[rejected_bonf] == 1)  # 正确拒绝的备择假设数
  R_bonf <- length(rejected_bonf)                 # 总共拒绝的假设数
  
  fwer_bonf[i] <- as.numeric(V_bonf > 0)          # 至少一个假阳性
  fdr_bonf[i] <- ifelse(R_bonf > 0, V_bonf / R_bonf, 0)  # 错误发现率
  tpr_bonf[i] <- S_bonf / n_alt                   # 真阳性率
  
  # B-H 校正
  p_bh <- p.adjust(p_values, method = "BH")
  rejected_bh <- which(p_bh < alpha)
  
  V_bh <- sum(true_labels[rejected_bh] == 0)
  S_bh <- sum(true_labels[rejected_bh] == 1)
  R_bh <- length(rejected_bh)
  
  fwer_bh[i] <- as.numeric(V_bh > 0)
  fdr_bh[i] <- ifelse(R_bh > 0, V_bh / R_bh, 0)
  tpr_bh[i] <- S_bh / n_alt
  
  # 可视化 p 值分布（密度图）
  if (i == 1) { # 只绘制第一次模拟的 p 值分布
    # 创建一个数据框，用于 ggplot2 绘图
    p_data <- data.frame(
      p_value = p_values,
      hypothesis = factor(true_labels, labels = c("Null Hypothesis", "Alternative Hypothesis"))
    )
    
    # 绘制密度图
    p_plot <- ggplot(p_data, aes(x = p_value, fill = hypothesis, color = hypothesis)) +
      geom_density(alpha = 0.4) +  # 绘制密度曲线
      geom_vline(xintercept = alpha, linetype = "dashed", color = "red") +  # 显著性水平线
      labs(title = "P-value Density for Null and Alternative Hypotheses",
           x = "P-value",
           y = "Density") +
      theme_minimal()
    
    # 打印图形
    print(p_plot)
  }
  
  # 可视化 p 值分布（散点图）
  if (i == m) { # 只绘制第m次模拟的 p 值分布
    # 创建一个数据框，用于 ggplot2 绘图
    p_data <- data.frame(
      p_value = p_values,
      hypothesis = factor(true_labels, labels = c("Null Hypothesis", "Alternative Hypothesis"))
    )
  
    # 绘制散点图
    p_plot <- ggplot(p_data, aes(x = hypothesis, y = p_value, color = hypothesis)) +
      geom_jitter(width = 0.2, height = 0) +  # 使用散点抖动效果
      geom_hline(yintercept = alpha, linetype = "dashed", color = "red") +  # 显著性水平线
      labs(title = "P-value Distribution for Null and Alternative Hypotheses",
         x = "Hypothesis Type",
         y = "P-value") +
      theme_minimal()
  
    # 打印图形
    print(p_plot)
  }
}


## -----------------------------------------------------------------------------
# 计算平均值
results["FWER", "Bonferroni校正"] <- mean(fwer_bonf)
results["FDR", "Bonferroni校正"] <- mean(fdr_bonf)
results["TPR", "Bonferroni校正"] <- mean(tpr_bonf)
results["FWER", "B-H校正"] <- mean(fwer_bh)
results["FDR", "B-H校正"] <- mean(fdr_bh)
results["TPR", "B-H校正"] <- mean(tpr_bh)
print(results)

## -----------------------------------------------------------------------------
# 加载 boot 包和数据
library(boot)
data("aircondit")
x <- aircondit$hours
 # 样本大小
n <- length(x)
# 计算 λ 的 MLE
lambda_hat <- n / sum(x)
# 定义用于 Bootstrap 的统计量函数
boot_lambda <- function(data, indices) {
  # 从数据中抽样
  sample_data <- data[indices]
  # 计算 λ 的 MLE
  n_boot <- length(sample_data)
  sum_boot <- sum(sample_data)
  lambda_boot <- n_boot / sum_boot
  return(lambda_boot)
}
# 设置 Bootstrap 重复次数
R <- 10^5
# 进行 Bootstrap
set.seed(123)
boot_results <- boot(data = x, statistic = boot_lambda, R = R)

# 估计偏差
bias_lambda <- mean(boot_results$t) - lambda_hat
cat(sprintf("偏差估计为：%.6f\n", bias_lambda))
se_lambda <- sd(boot_results$t)
cat(sprintf("标准差估计为：%.6f\n", se_lambda))

## -----------------------------------------------------------------------------
# 将 λ 转换为 θ = 1/λ
boot_theta <- boot_results
boot_theta$t <- 1 / boot_theta$t
theta_hat <- 1 / lambda_hat

# 定义一个用于 Bootstrap 方法的统计量函数，用于计算平均故障间隔时间（θ = 1/λ）
mean_time_boot <- function(data, indices) {
  # 从原始数据中抽取 Bootstrap 样本
  # 'data' 是原始数据集，'indices' 是抽样索引（包含有放回的抽样位置）
  sample_data <- data[indices]
  # 计算抽样样本的大小（观测值数量）
  n_boot <- length(sample_data)
  # 计算抽样样本中所有故障时间的总和
  sum_boot <- sum(sample_data)
  # 计算 Bootstrap 样本的最大似然估计（MLE）λ的估计值
  # λ 的 MLE 公式为 λ_hat = n / Σx_i，其中 n 是样本大小，Σx_i 是故障时间总和
  lambda_hat <- n_boot / sum_boot
  # 计算平均故障间隔时间 θ 的估计值，θ = 1/λ
  theta_boot <- 1 / lambda_hat
  # 返回 θ 的估计值作为 Bootstrap 样本的统计量
  return(theta_boot)
}
set.seed(133)  # 设置随机数种子，保证结果可重复
# 数据：空调设备的故障时间（小时）
failure_times <- c(3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487)
 # 使用bootstrap方法计算
bootstrap_mean_time_results <- boot(data = failure_times, statistic = mean_time_boot, R = 10000)

# 标准正态法置信区间
norm_ci <- boot.ci(bootstrap_mean_time_results, type = "norm")
 # 基本法置信区间
basic_ci <- boot.ci(bootstrap_mean_time_results, type = "basic")
 # 百分位法置信区间
perc_ci <- boot.ci(bootstrap_mean_time_results, type = "perc")
 # BCa 法置信区间
bca_ci <- boot.ci(bootstrap_mean_time_results, type = "bca")

ci_methods <- c("标准正态法", "基本法", "百分位法", "BCa 法")
lower_bounds <- c(norm_ci$normal[2], basic_ci$basic[4], perc_ci$percent[4], bca_ci$bca[4])
upper_bounds <- c(norm_ci$normal[3], basic_ci$basic[5], perc_ci$percent[5], bca_ci$bca[5])
data.frame(ci_methods,lower_bounds,upper_bounds)

## -----------------------------------------------------------------------------
library(bootstrap)
data = data(scor)

## -----------------------------------------------------------------------------
# 主成分分析
original_pca <- prcomp(scor, center = TRUE, scale. = TRUE)
# 提取前五个主成分的特征值
explained_variances <- (original_pca$sdev)^2
# 计算原始数据的PC1在前五个主成分的比例
original_ratio <- explained_variances[1] / sum(explained_variances[1:5])

## -----------------------------------------------------------------------------
# 假设scor是一个数据框，含有数据集的样本
set.seed(123)  # 设置随机种子
n <- nrow(scor)  # 获取数据集的行数
# 重复次数
N = 1000
boot_ratios <- numeric(N)

for (i in 1:N){
  # 有放回的抽样
  bootstrap_sample <- scor[sample(1:n, n, replace = TRUE),]
  # 执行主成分分析
  boot_pca <- prcomp(bootstrap_sample, center = TRUE, scale. = TRUE)
  # 提取前五个主成分的特征值
  explained_variances <- (boot_pca$sdev)^2
  # 计算原始数据的PC1在前五个主成分的比例
  boot_ratios[i] <- explained_variances[1] / sum(explained_variances[1:5])
}
# 计算偏差
bias <- mean(boot_ratios) - original_ratio
# 计算标准误
standard_error <- sd(boot_ratios)

cat("原始数据的样本估计值为", original_ratio,'\n')
cat("bootstrap方法的偏差为", bias,'\n')
cat("bootstrap方法的标准误为",standard_error)

## -----------------------------------------------------------------------------
library(DAAG)
attach(ironslag)

## -----------------------------------------------------------------------------
a <- seq(10, 40, .1)
par(mfrow = c(2, 2))

L1 <- lm(magnetic ~ chemical)
plot(chemical, magnetic, main="Linear", pch=16)
yhat1 <- L1$coef[1] + L1$coef[2] * a
lines(a, yhat1, lwd=2)
L2 <- lm(magnetic ~ chemical + I(chemical^2))
plot(chemical, magnetic, main="Quadratic", pch=16)
yhat2 <- L2$coef[1] + L2$coef[2] * a + L2$coef[3] * a^2
lines(a, yhat2, lwd=2)
L3 <- lm(log(magnetic) ~ chemical)
plot(chemical, magnetic, main="Exponential", pch=16)
logyhat3 <- L3$coef[1] + L3$coef[2] * a
yhat3 <- exp(logyhat3)
lines(a, yhat3, lwd=2)
L4 <- lm(magnetic ~ chemical + I(chemical^2)+ I(chemical^3))
plot(chemical, magnetic, main="cubic", pch=16)
yhat4 <- L4$coef[1] + L4$coef[2] * a + L4$coef[3] * a^2 + L4$coef[4] * a^3
lines(a, yhat4, lwd=2)

## -----------------------------------------------------------------------------
n <- length(magnetic)
e1 <- e2 <- e3 <- e4 <- numeric(n)
for (k in 1:n) {
  y <- magnetic[-k]
  x <- chemical[-k]
  
  J1 <- lm(y ~ x)
  yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
  e1[k] <- magnetic[k] - yhat1
  
  J2 <- lm(y ~ x + I(x^2))
  yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] + J2$coef[3] * chemical[k]^2
  e2[k] <- magnetic[k] - yhat2
  
  J3 <- lm(log(y) ~ x)
  logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
  yhat3 <- exp(logyhat3)
  e3[k] <- magnetic[k] - yhat3
  
  J4 <- lm(y ~ x + I(x^2) + I(x^3))
  yhat4 <- J4$coef[1] + J4$coef[2] * chemical[k] + J4$coef[3] * chemical[k]^2 + J4$coef[4] * chemical[k]^3
  e4[k] <- magnetic[k] - yhat4
}
# 根据误差比较
cat("均值为", c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2)),'\n')
cat("R^2为", c(summary(J1)$r.squared,summary(J2)$r.squared,summary(J3)$r.squared,summary(J4)$r.squared))

## -----------------------------------------------------------------------------
data("chickwts")
# 选择两个饲料组作为样本
feed1 <- chickwts$weight[chickwts$feed == "horsebean"]
feed2 <- chickwts$weight[chickwts$feed == "linseed"]

# Cramér-von Mises统计量计算函数
cramer_von_mises <- function(x, y) {
  n1 <- length(x)
  n2 <- length(y)
  # 计算经验分布函数
  ecdf_x <- ecdf(x)
  ecdf_y <- ecdf(y)
  # 合并并排序样本
  combined <- sort(c(x, y))
  # 计算Cramér-von Mises统计量
  stat <- sum((ecdf_x(combined) - ecdf_y(combined))^2) * (n1 * n2 / (n1 + n2))
  return(stat)
}
# 置换检验函数
cramer_von_mises_permutation_test <- function(x, y, n_permutations = 1000) {
  # 计算原始统计量
  observed_stat <- cramer_von_mises(x, y)
  # 合并样本
  combined <- c(x, y)
  # 初始化置换统计量向量
  permuted_stats <- numeric(n_permutations)
  # 进行置换
  for (i in 1:n_permutations) {
    permuted_sample <- sample(combined)
    new_x <- permuted_sample[1:length(x)]
    new_y <- permuted_sample[(length(x) + 1):length(combined)]
    permuted_stats[i] <- cramer_von_mises(new_x, new_y)
  }
  # 计算p值
  p_value <- mean(permuted_stats >= observed_stat)
  return(list(observed_stat = observed_stat, p_value = p_value))
}
# 进行置换检验
set.seed(123)
result <- cramer_von_mises_permutation_test(feed1, feed2, n_permutations = 1000)

# 输出结果
cat("Cramer-von Mises检验统计量结果:", result$observed_stat, "\n")
cat("P-value:", result$p_value, "\n")

## -----------------------------------------------------------------------------
# 生成示例数据
set.seed(123)
x <- rnorm(50)
y <- rnorm(50)

# 定义置换检验函数
spearman_permutation_test <- function(x, y, n_permutations = 1000) {
  # 计算原始样本的Spearman相关系数
  observed_stat <- cor(x, y, method = "spearman")
  # 初始化置换统计量向量
  permuted_stats <- numeric(n_permutations)
  # 进行置换
  for (i in 1:n_permutations) {
    permuted_y <- sample(y)  # 随机打乱y
    permuted_stats[i] <- cor(x, permuted_y, method = "spearman")
  }
  # 计算置换p值
  p_value <- mean(abs(permuted_stats) >= abs(observed_stat))
  return(list(observed_stat = observed_stat, p_value = p_value))
}

# 运行置换检验
set.seed(123)
result_permutation <- spearman_permutation_test(x, y, n_permutations = 1000)

# 计算Spearman相关系数和p值
result_cor_test <- cor.test(x, y, method = "spearman")

# 输出结果
cat("Spearman相关系数:", result_permutation$observed_stat, "\n")
cat("Permutation test p值:", result_permutation$p_value, "\n")
cat("cor.test p值:", result_cor_test$p.value, "\n")
cat("p值的差距为：", abs(result_permutation$p_value-result_cor_test$p.value))

## -----------------------------------------------------------------------------
# 设置目标分布的密度函数
target_density <- function(x) {
  return (1/(pi*(1+x^2)))
}

## -----------------------------------------------------------------------------
# iter迭代次数，start是起始点
metropolis_hastings <- function(iter, start, proposal_sd) {
  k = 0
  samples <- numeric(iter) # 保存样本
  samples[1] <- start      # 初始值 
  for (i in 2:iter) {
    # 从提议分布中采样，生成新状态
    proposal <- rnorm(1, mean = samples[i - 1], sd = proposal_sd)
    # 计算接受率
    acceptance_ratio <- target_density(proposal) / target_density(samples[i - 1])
    # 以接受率的最小值为界进行判断
    if (runif(1) < acceptance_ratio) {
      samples[i] <- proposal # 接受提议样本
    } else {
      samples[i] <- samples[i - 1] # 拒绝提议，保留当前样本
      k <- k+1
    }
  }
  return(samples)
}

## -----------------------------------------------------------------------------
iterations <- 10000      # 迭代次数
start_value <- 0         # 起始点
proposal_sd <- 1         # 提议分布的标准差
set.seed(123)
samples <- metropolis_hastings(iterations, start_value, proposal_sd)

index <- 1000:10000
y1 <- samples[index]
plot(index, y1, type="l", main="", ylab="x")

## -----------------------------------------------------------------------------
sample_deciles <- quantile(samples[1000:10000], probs = seq(0, 1, by = 0.1))

# 设置位置参数和尺度参数
location <- 0  # 默认位置参数
scale <- 1     # 默认尺度参数

# 计算十分位数
decile_probs <- seq(0, 1, by = 0.1)  # 十分位数对应的概率
deciles <- qcauchy(decile_probs, location, scale)

# 输出结果
results <- data.frame(
  Sample_Quantiles = sample_deciles,
  Theoretical_Quantiles = deciles,
  abs_diff = abs(sample_deciles-deciles)
)
print(results)

## -----------------------------------------------------------------------------
N = 10000
X = matrix(0,N,2)
a = 2
b = 3
n = 10
burn = 1000

# 链的初始值
x0 <- rbinom(1, n, 0.5)  # 随机初始化 x
y0 <- rbeta(1,a,b)  # 随机初始化 y
X[1, ] <- c(x0, y0)

# 进行Gibbs采样
for (i in 2:N) {
  y <- X[i-1, 2]
  X[i, 1] <- rbinom(1,n,y)
  x <- X[i, 1]
  X[i, 2] <- rbeta(1,x+a,n-x+b)
}

# 丢弃前1000个数据
b <- burn + 1
x <- X[b:N, ]

# 可视化Gibbs采样
plot(X, main="", cex=.5, xlab='x',
  ylab='y', ylim=range(x[,2]))

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi)
  #row means
  B <- n * var(psi.means)
  #between variance est.
  psi.w <- apply(psi, 1, "var") #within variances
  W <- mean(psi.w)
  #within est.
  v.hat <- W*(n-1)/n + (B/n)
  #upper variance est.
  r.hat <- v.hat / W
  #G-R statistic
  return(r.hat)
}

set.seed(123)
sigma <- 1 #parameter of proposal distribution
k <- 4
#number of chains to generate
n <- 15000
#length of chains
b <- 1000
#burn-in length
#choose overdispersed initial values
x0 <- c(-1, 0, 1, 2)
#generate the chains
X <- matrix(0, nrow=k, ncol=n)

for (i in 1:k)
  X[i, ] <- metropolis_hastings(n, x0[i], sigma)

#compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
# 进行矩阵变化，psi是一个矩阵
for (i in 1:nrow(psi))
  psi[i,] <- psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi))

#plot psi for the four chains
par(mfrow=c(2,2))
for (i in 1:k)
  plot(psi[i, (b+1):n], type="l",
    xlab=i, ylab=bquote(psi))
par(mfrow=c(1,1)) #restore default

#plot the sequence of R-hat statistics
rhat <- rep(0, n)
for (j in (b+1):n)
  rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
abline(h=1.44, lty=2)

## -----------------------------------------------------------------------------
k = 100
a = c(1,2,3)
d = 10

# 计算第k项的函数
numk = function(k, a){
  l1 = (-1)^k/(factorial(k)*2^k)
  l2 = exp(log(norm(a, type = "2")^(2*k+2))-log(2*k+1)-log(2*k+2))
  l3 = exp(lgamma((d+1)/2)+lgamma(k+3/2)-lgamma(k+d/2+1))
  result = l1*l2*l3
  return(result)
}

## -----------------------------------------------------------------------------
num = 100
for (num in 1:20){
  k = 0:num
  odd = 2*(0:num)+1
  even = 2*(0:num)+2
  i = rep(c(1,-1),length = num+1)
  a = c(1,2)
  d = 10
  l3 = exp(lgamma((d+1)/2)+lgamma(k+3/2)-lgamma(k+d/2+1))
  l1 = i/factorial(k)/2^k
  l2 = norm(a, type = "2")^even/odd/even
  print(sum(l1*l2*l3))
}

## -----------------------------------------------------------------------------
c_k = function(a,k){
  result = sqrt(a^2*k/(k+1-a^2))
  return(result)
}

integral = function(u,k){
  (1+u^2/(k-1))^(-k/2)
}

expr = function(a,k){
  integral1 = function(u){
    integral(u,k)
  }
  c = c_k(a,k-1)
  l1 = 2*exp(lgamma(k/2)-lgamma((k-1)/2))/sqrt(pi*(k-1))
  res = l1 * integrate(integral1,lower = 0, upper = c)$value
  return(res)
}

# 进行迭代
result1 = function(a){
  expr(a,k) - expr(a,k+1)
}

resl = function(k){
  r = uniroot(result1, c(0.01,sqrt(k)-0.01))$root
  return(r)
}

for (k in c(4:15)){
  print(resl(k))
}

## -----------------------------------------------------------------------------
# 初始化观测数据
Y <- c(0.54, 0.48, 0.33, 0.43, 1.00, 1.00, 0.91, 1.00, 0.21, 0.85)
tau <- 1  # 截断值

# E-M算法的实现
EM_algorithm <- function(Y, tau, max_iter = 100, tol = 1e-6) {
  n <- length(Y)
  
  # 初始化lambda
  lambda <- 1 
  
  # 进行E-M迭代
  for (iter in 1:max_iter) {
    # E步：计算潜在变量T的期望
    T_expect <- rep(NA, n)
    for (i in 1:n) {
      if (Y[i] < tau) {
        T_expect[i] <- Y[i]
      } else {
        T_expect[i] <- tau + 1 / lambda
      }
    }
    
    # M步：最大化对数似然函数，更新lambda
    # 计算新的lambda
    lambda_new <- n / sum(T_expect)
    
    # 检查是否收敛
    if (abs(lambda_new - lambda) < tol) {
      break
    }
    
    # 更新lambda
    lambda <- lambda_new
  }
  
  return(lambda)
}

# 运行E-M算法来估计lambda
lambda_EM <- EM_algorithm(Y, tau)
cat("E-M估计的lambda值:", lambda_EM, "\n")

# MLE的估计
MLE_lambda <- 2*sum(Y[Y < tau]) / length(Y[Y < tau])
cat("MLE估计的lambda值:", MLE_lambda, "\n")

## -----------------------------------------------------------------------------
library(boot)
A1 = rbind(c(2, 1, 1), c(1, -1, 3))
b1 <- c(2, 3)
a <- c(4, 2, 9)
simplex(a = a, A1 = A1, b1 = b1, maxi = FALSE)

## -----------------------------------------------------------------------------
data(mtcars)
formulas <- list(
  mpg ~ disp, # 探索线性关系
  mpg ~ I(1 / disp), # 探索倒数关系
  mpg ~ disp + wt, # mpg与disp和wt的线性关系
  mpg ~ I(1 / disp) + wt # mpg 与 disp 的倒数和 wt 的线性关系
)
# 使用 lapply() 循环应用每个公式
models_ex3 <- lapply(formulas, function(formula) {
  lm(formula, data = mtcars)  # 拟合线性回归模型
})
# 打印模型摘要
models_ex3

## -----------------------------------------------------------------------------
# 创建一个空列表来存储模型
models_loop_ex3 <- list()

# 使用for循环拟合线性模型
for (i in seq_along(formulas)) {
  models_loop_ex3[[i]] <- lm(formulas[[i]], data = mtcars)
}

print(models_loop_ex3)

## -----------------------------------------------------------------------------
set.seed(123)
# 定义 bootstrap 列表
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), replace = TRUE)
  mtcars[rows, ]
})

# 使用线性模型拟合
fit_model <- function(data) {
  lm(mpg ~ disp, data = data)
}

# 使用 lapply 将 fit_model 函数应用于每个 bootstrap 复制
models_ex4 <- lapply(bootstraps, fit_model)

# 创建一个空列表来存储模型
models_loop_ex4 <- list()

j = 0
for (i in bootstraps){
  j = j+1
  models_loop_ex4[[j]] = lm(i$mpg ~ i$disp,data = i)
}

print(models_ex4)
print(models_loop_ex4)


## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared
# 对于ex3的模型
rex3 = lapply(models_ex3,rsq)
rex3_loop = lapply(models_loop_ex3,rsq)
rex4 = lapply(models_ex4,rsq)
rex4_loop = lapply(models_loop_ex4,rsq)
df1 = data.frame(
  Name = c("mpg ~ disp","mpg ~ I(1 / disp)","mpg ~ disp + wt","mpg ~ I(1 / disp) + wt"),
  lapply_R2 = as.double(unlist(rex3)),
  loop_R2 = as.double(unlist(rex3_loop))
)
df2 = data.frame(
  Name = c("bootstraps1","bootstraps2","bootstraps3","bootstraps4","bootstraps5","bootstraps6",
           "bootstraps7","bootstraps8","bootstraps9","bootstraps10"),
  lapply_R2 = as.double(unlist(rex4)),
  loop_R2 = as.double(unlist(rex4_loop))
)
print(df1)
print(df2)

## -----------------------------------------------------------------------------
trials <- replicate(
  100,
  t.test(rpois(10, 10), rpois(7, 10)),
  simplify = FALSE
)

## -----------------------------------------------------------------------------
sapply(trials, function(ttest) ttest$p.value)

## -----------------------------------------------------------------------------
my_list <- list(1, 2, 3, 4, 5)

# 使用 Map() 和 vapply() 实现 lapply() 的功能
result <- Map(function(x) vapply(x, function(y) y^2, numeric(1)), my_list)

data.frame(
  before = as.double(unlist(my_list)),
  after = as.double(unlist(result))
)

## -----------------------------------------------------------------------------
x = c(24,12,18)
y = c(31,71,45)
data = data.frame(x,y)
# 计算行列和
datasum = sum(data)
# 按列求和
sr <- rowSums(data)
# 按行求和
sc <- colSums(data)
# 计算对应的期望值
E = outer(sc,sr)/datasum
# 计算实际观测值与期望值的偏差和
chisq = sum((t(as.matrix(data))-E)^2/E)
chisq

## -----------------------------------------------------------------------------
Fruit = c("Apple", "Banana", "Apple", "Orange", "Apple", "Banana", "Orange","Apple", "Banana", "Apple", "Orange", "Apple",
          "Apple", "Orange", "Apple")
Gender = c("Male", "Female", "Female", "Male", "Female", "Male", "Female","Female", "Male", "Female","Female", "Male",
           "Female","Male","Male")

data = data.frame(Gender, Fruit)
# 确定出现的元素
fruituni = unique(Fruit)
genderuni = unique(Gender)

# 计算出现的数量
count_occurrences <- function(x, B) {
  sum(x == B)
}

Fru_counts_male <- sapply(fruituni, count_occurrences, B = subset(data, Gender == "Male"))
Fru_counts_female <- sapply(fruituni, count_occurrences, B = subset(data, Gender == "Female"))

data.frame(Fru_counts_female,Fru_counts_male)

## -----------------------------------------------------------------------------
library(Rcpp)
library(SA24204125)

# 设置参数
N <- 10000
a <- 2
b <- 3
n <- 10
burn <- 1000

# 运行 Gibbs 采样
result <- gibbsCpp(N, a, b, n, burn)

# 查看结果
plot(result, main="", cex=.5, xlab='x',ylab='y', ylim=range(result[,2]))

## -----------------------------------------------------------------------------
N = 10000
a = 2
b = 3
n = 10
burn = 1000
set.seed(12345)

gibbs_R = function(N,a,b,n,burn){
  X = matrix(0,N,2)
  # 链的初始值
  x0 <- rbinom(1, n, 0.5)  # 随机初始化 x
  y0 <- rbeta(1,a,b)  # 随机初始化 y
  X[1, ] <- c(x0, y0)
  
  # 进行Gibbs采样
  for (i in 2:N) {
    y <- X[i-1, 2]
    X[i, 1] <- rbinom(1,n,y)
    x <- X[i, 1]
    X[i, 2] <- rbeta(1,x+a,n-x+b)
  }
  
  # 丢弃前1000个数据
  b <- burn + 1
  x <- X[b:N, ]
}

X = gibbs_R(N,a,b,n,burn)
# 可视化Gibbs采样
plot(X, main="", cex=.5, xlab='x',ylab='y', ylim=range(X[,2]))

## -----------------------------------------------------------------------------
par(mfrow = c(1, 2))
qqplot(X[,1], result[,1], main = "QQ Plot of X", xlab = "X from R", ylab = "X from Cpp")
abline(0, 1, col = "red", lwd = 2)  # 添加 y = x 的参考线
qqplot(X[,2], result[,2], main = "QQ Plot of Y", xlab = "Y from R", ylab = "Y from Cpp")
abline(0, 1, col = "red", lwd = 2)  # 添加 y = x 的参考线

## -----------------------------------------------------------------------------
library(microbenchmark)
library(SA24204125)
ts<-microbenchmark(Rgibbs = gibbs_R(N,a,b,n,burn), Cppgibbs = gibbsCpp(N, a, b, n, burn))
summary(ts)[,c(1,3,5,6)]

