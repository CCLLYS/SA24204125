#' @import microbenchmark
#' @import deSolve
#' @import ggplot2
#' @import boot
#' @import bootstrap
#' @import DAAG
#' @importFrom pracma trapz
#' @importFrom Rcpp evalCpp
#' @importFrom stats rnorm rgamma
#' @useDynLib SA24204125
NULL

#' @title The fourth-order Runger-Kutta method
#' @description This model is used to calculate the initial values of ordinary differential equations
#' @param fun the ordinary differential equations
#' @param a the upper bound of the interval
#' @param b the lower bound of the interval
#' @param n the length of the interval
#' @param y0 initial value
#' @return The solution of the equation on the interval [a,b], the data type is dataframe.
#' @examples
#' \dontrun{
#' fun <- function(x,y){
#' x-y+1
#' }
#' a=0
#' b=1
#' n=10
#' y0=1
#' result <- RK4ODE(fun, a, b, n, y0)
#' }
#' @export
RK4ODE <- function(fun, a, b, n, y0){
  x <- numeric()
  y <- numeric()
  h <- (b-a)/n
  x[1] <- a
  y[1] <- y0
  for(k in 1:n){
    K1 <- fun(x[k],y[k])
    K2 <- fun(x[k]+h/2,y[k]+h/2*K1)
    K3 <- fun(x[k]+h/2,y[k]+h/2*K2)
    K4 <- fun(x[k]+h,y[k]+h*K3)
    y[k+1] <- y[k]+h/6*(K1+2*K2+2*K3+K4)
    x[k+1] <- x[k]+h
  }
  result <- data.frame(x=x,y=y)
  return(result)
}

#' @title The Romberg quadrature method
#' @description This method is used to solve numerical integrals
#' @param f the function
#' @param a the upper bound of the interval
#' @param b the lower bound of the interval
#' @param TOL integration error limit
#' @return A list is returned, with the first value s being the integral result, 
#' the second value R being the iterative process, and the third value k being the number of iterations
#' @examples
#' \dontrun{
#' f <- function(x) exp(x)*x^2+exp(x)*sin(x)+exp(x)*cos(x)+2**x  # 积分函数
#' a <- 0  # 积分下限
#' b <- 1  # 积分上限
#' TOL <- 1e-16  # 误差容限
#' result <- Romberg(f, a, b, TOL)
#' }
#' @export
Romberg <- function(f, a, b, TOL) {
  # 初始化参数
  n <- 1
  h <- b - a
  delt <- 1  # 误差初值
  k <- 0  # 迭代次数
  R <- matrix(0, nrow = 10, ncol = 10)  # 创建一个4x4的矩阵
  R[1, 1] <- h / 2 * (f(a) + f(b))
  
  while (delt > TOL && k<10) {
    # 如果两次计算的误差大于给定误差则进入循环
    k <- k + 1
    h <- h / 2
    s <- 0
    for (j in 1:n) {
      x <- a + h * (2 * j - 1)
      s <- s + f(x)
    }
    R[k + 1, 1] <- R[k, 1] / 2 + h * s
    n <- 2 * n
    for (i in 1:k) {
      R[k + 1, i + 1] <- ((4^i) * R[k + 1, i] - R[k, i]) / (4^i - 1)
    }
    delt <- abs(R[k + 1, k] - R[k + 1, k + 1])  # 计算前后两次的值
  }
  s <- R[k + 1, k + 1]
  return(list(s = s, R = R, k = k))
}



