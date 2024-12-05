#include <Rcpp.h>
using namespace Rcpp;

//' @title The fourth-order Runger-Kutta method using Cpp
//' @description This model is used to calculate the initial values of ordinary differential equations
//' @param fun the ordinary differential equations
//' @param a the upper bound of the interval
//' @param b the lower bound of the interval
//' @param n the length of the interval
//' @param y0 initial value
//' @return The solution of the equation on the interval [a,b], the data type is dataframe.
//' @examples
//' \dontrun{
//' fun <- function(x,y){
//' x-y+1
//' }
//' a=0
//' b=1
//' n=10
//' y0=1
//' library(Rcpp)
//' sourceCpp(paste0('D:/R/Rwork/Statistics_computing/RK.cpp'))
//' RK4ODECpp = RK4ODECpp(fun, a, b, n, y0)
//' }
//' @export
// [[Rcpp::export]]
DataFrame RK4ODECpp(Function fun, double a, double b, int n, double y0) {
  NumericVector x(n + 1);
  NumericVector y(n + 1);
  double h = (b - a) / n;
  
  x[0] = a;
  y[0] = y0;
  
  for (int k = 0; k < n; ++k) {
    double K1 = as<double>(fun(x[k], y[k]));
    double K2 = as<double>(fun(x[k] + h / 2, y[k] + h / 2 * K1));
    double K3 = as<double>(fun(x[k] + h / 2, y[k] + h / 2 * K2));
    double K4 = as<double>(fun(x[k] + h, y[k] + h * K3));
    
    y[k + 1] = y[k] + h / 6 * (K1 + 2 * K2 + 2 * K3 + K4);
    x[k + 1] = x[k] + h;
  }
  return DataFrame::create(Named("x") = x, Named("y") = y);
}

//' @title The Romberg quadrature method using Cpp
//' @description This method is used to solve numerical integrals
//' @param f the function
//' @param a the upper bound of the interval
//' @param b the lower bound of the interval
//' @param TOL integration error limit
//' @return A list is returned, with the first value s being the integral result, 
//' the second value R being the iterative process, and the third value k being the number of iterations
//' @examples
//' \dontrun{
//' f <- function(x) exp(x)*x^2+exp(x)*sin(x)+exp(x)*cos(x)+2**x  # 积分函数
//' a <- 0  # 积分下限
//' b <- 1  # 积分上限
//' TOL <- 1e-16  # 误差容限
//' resultCpp <- RombergCpp(f, a, b, TOL)
//' }
//' @export
// [[Rcpp::export]]
List RombergCpp(Function f, double a, double b, double TOL) {
  int n = 1;                 // 初始区间数
  double h = b - a;          // 初始步长
  double delt = 1.0;         // 误差初值
  int k = 0;                 // 迭代次数
  int max_iter = 10;         // 最大迭代次数
  NumericMatrix R(max_iter, max_iter); // 初始化R矩阵，大小为max_iter x max_iter
  
  // 初始化 R[0, 0]
  R(0, 0) = h / 2 * (as<double>(f(a)) + as<double>(f(b)));
  
  while (delt > TOL && k < max_iter - 1) {
    k++;
    h /= 2;
    double sum = 0.0;
    // 计算新的点
    for (int j = 1; j <= n; j++) {
      double x = a + h * (2 * j - 1);
      sum += as<double>(f(x));
    }
    // 更新 R[k, 0]
    R(k, 0) = R(k - 1, 0) / 2 + h * sum;
    n *= 2;
    // 计算更高阶结果
    for (int i = 1; i <= k; i++) {
      R(k, i) = (pow(4, i) * R(k, i - 1) - R(k - 1, i - 1)) / (pow(4, i) - 1);
    }
    // 更新误差
    delt = fabs(R(k, k - 1) - R(k, k));
  }
  // 最终结果
  double s = R(k, k);
  return List::create(
    Named("result") = s,
    Named("iterations") = k,
    Named("R_matrix") = R
  );
}

 
 