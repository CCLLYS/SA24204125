#include <Rcpp.h>
using namespace Rcpp;

//' @title The gibbs function
//' @description A function is used in homework
//' @param N size
//' @param a the param of beta
//' @param b the param of beta
//' @param n the size of x
//' @param burn start from burn.
//' @return gibbs result
//' @examples
//' \dontrun{
//' N <- 10000
//' a <- 2
//' b <- 3
//' n <- 10
//' burn <- 1000
//' result <- gibbsCpp(N, a, b, n, burn)
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix gibbsCpp(int N, int a, int b, int n, int burn) {
  NumericMatrix X(N, 2); // 初始化存储结果的矩阵
  
  // 随机初始化 x 和 y
  X(0, 0) = R::rbinom(n, 0.5);  // 初始化 x
  X(0, 1) = R::rbeta(a, b);     // 初始化 y
  
  for (int i = 1; i < N; ++i) {
    double y = X(i - 1, 1);                   // 获取上一轮 y 值
    X(i, 0) = R::rbinom(n, y);               // 采样 x
    int x = X(i, 0);                         // 获取新的 x 值
    X(i, 1) = R::rbeta(x + a, n - x + b);    // 采样 y
  }
  
  // 丢弃前 burn 个样本
  NumericMatrix result(N - burn, 2);
  for (int i = burn; i < N; ++i) {
    result(i - burn, _) = X(i, _);
  }
  
  return result;
}