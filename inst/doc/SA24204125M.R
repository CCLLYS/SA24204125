## ----eval = FALSE, include=FALSE----------------------------------------------
#  # 安装包
#  devtools::install()

## ----warning = FALSE, include=FALSE-------------------------------------------
library(SA24204125)
library(Rcpp)
library(microbenchmark)
library(deSolve)
trapz = pracma::trapz

## ----eval=FALSE---------------------------------------------------------------
#  RK4ODE <- function(fun, a, b, n, y0){
#    x <- numeric()
#    y <- numeric()
#    h <- (b-a)/n
#    x[1] <- a
#    y[1] <- y0
#    for(k in 1:n){
#      K1 <- fun(x[k],y[k])
#      K2 <- fun(x[k]+h/2,y[k]+h/2*K1)
#      K3 <- fun(x[k]+h/2,y[k]+h/2*K2)
#      K4 <- fun(x[k]+h,y[k]+h*K3)
#      y[k+1] <- y[k]+h/6*(K1+2*K2+2*K3+K4)
#      x[k+1] <- x[k]+h
#    }
#    result <- data.frame(x=x,y=y)
#    return(result)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  DataFrame RK4ODECpp(Function fun, double a, double b, int n, double y0) {
#    NumericVector x(n + 1);
#    NumericVector y(n + 1);
#    double h = (b - a) / n;
#  
#    x[0] = a;
#    y[0] = y0;
#  
#    for (int k = 0; k < n; ++k) {
#      double K1 = as<double>(fun(x[k], y[k]));
#      double K2 = as<double>(fun(x[k] + h / 2, y[k] + h / 2 * K1));
#      double K3 = as<double>(fun(x[k] + h / 2, y[k] + h / 2 * K2));
#      double K4 = as<double>(fun(x[k] + h, y[k] + h * K3));
#  
#      y[k + 1] = y[k] + h / 6 * (K1 + 2 * K2 + 2 * K3 + K4);
#      x[k + 1] = x[k] + h;
#    }
#    return DataFrame::create(Named("x") = x, Named("y") = y);
#  }

## -----------------------------------------------------------------------------
fun <- function(x,y){x-y+1}
a=0
b=1
n=100
y0=1
result <- RK4ODE(fun, a, b, n, y0)
plot(result$x,result$y,type='l',xlab='x',ylab='y',col='skyblue')
curve(x+exp(-x),add=T,col='red')
print(head(result))

## ----eval=FALSE---------------------------------------------------------------
#  Romberg <- function(f, a, b, TOL) {
#    # 初始化参数
#    n <- 1
#    h <- b - a
#    delt <- 1  # 误差初值
#    k <- 0  # 迭代次数
#    R <- matrix(0, nrow = 10, ncol = 10)  # 创建一个4x4的矩阵
#    R[1, 1] <- h / 2 * (f(a) + f(b))
#  
#    while (delt > TOL && k<10) {
#      # 如果两次计算的误差大于给定误差则进入循环
#      k <- k + 1
#      h <- h / 2
#      s <- 0
#      for (j in 1:n) {
#        x <- a + h * (2 * j - 1)
#        s <- s + f(x)
#      }
#      R[k + 1, 1] <- R[k, 1] / 2 + h * s
#      n <- 2 * n
#      for (i in 1:k) {
#        R[k + 1, i + 1] <- ((4^i) * R[k + 1, i] - R[k, i]) / (4^i - 1)
#      }
#      delt <- abs(R[k + 1, k] - R[k + 1, k + 1])  # 计算前后两次的值
#    }
#    s <- R[k + 1, k + 1]
#    return(list(s = s, R = R, k = k))
#  }

## ----eval=FALSE---------------------------------------------------------------
#  List RombergCpp(Function f, double a, double b, double TOL) {
#    int n = 1;                 // 初始区间数
#    double h = b - a;          // 初始步长
#    double delt = 1.0;         // 误差初值
#    int k = 0;                 // 迭代次数
#    int max_iter = 10;         // 最大迭代次数
#    NumericMatrix R(max_iter, max_iter); // 初始化R矩阵，大小为max_iter x max_iter
#  
#    // 初始化 R[0, 0]
#    R(0, 0) = h / 2 * (as<double>(f(a)) + as<double>(f(b)));
#  
#    while (delt > TOL && k < max_iter - 1) {
#      k++;
#      h /= 2;
#      double sum = 0.0;
#      // 计算新的点
#      for (int j = 1; j <= n; j++) {
#        double x = a + h * (2 * j - 1);
#        sum += as<double>(f(x));
#      }
#      // 更新 R[k, 0]
#      R(k, 0) = R(k - 1, 0) / 2 + h * sum;
#      n *= 2;
#      // 计算更高阶结果
#      for (int i = 1; i <= k; i++) {
#        R(k, i) = (pow(4, i) * R(k, i - 1) - R(k - 1, i - 1)) / (pow(4, i) - 1);
#      }
#      // 更新误差
#      delt = fabs(R(k, k - 1) - R(k, k));
#    }
#    // 最终结果
#    double s = R(k, k);
#    return List::create(
#      Named("result") = s,
#      Named("iterations") = k,
#      Named("R_matrix") = R
#    );
#  }

## -----------------------------------------------------------------------------
f <- function(x) exp(x)*cos(x)+x^2+x^3*sin(x)  # 积分函数s
a <- 0  # 积分下限
b <- 1  # 积分上限
TOL <- 1e-16  # 误差容限
result <- Romberg(f, a, b, TOL)
cat("积分结果为", result$s)

## -----------------------------------------------------------------------------
fun1 <- function(x,y,parms){
  list(x-y+1)
}
ts<-microbenchmark(RK4ODE = RK4ODE(fun, a, b, n, y0), 
                   RK4ODECpp = RK4ODECpp(fun, a, b, n, y0),
                   odeR = ode(y =y0, times = seq(0,1,0.01), func = fun1, parms = NULL),
                   odeRRK4 = rk4(y =y0, times = seq(0,1,0.01), func = fun1, parms = NULL)
                   )
summary(ts)[,c(1,3,5,6)]

## -----------------------------------------------------------------------------
f <- function(x) exp(x)*cos(x)+x^2+x^3*sin(x)  # 积分函数s
a <- 0  # 积分下限
b <- 1  # 积分上限
TOL <- 1e-16  # 误差容限
print(c(Romberg(f, a, b, TOL)$s,RombergCpp(f, a, b, TOL)$result,integrate(f,lower = a, upper = b)$value,trapz(seq(a, b, length.out = 1000), f(seq(a, b, length.out = 1000)))))

## -----------------------------------------------------------------------------
ts<-microbenchmark(Romberg = Romberg(f, a, b, TOL), 
                   RombergCpp = RombergCpp(f, a, b, TOL),
                   resultR = integrate(f,lower = a, upper = b),
                   trapz = trapz(seq(a, b, length.out = 1000), f(seq(a, b, length.out = 1000)))
                   )
summary(ts)[,c(1,3,5,6)]

