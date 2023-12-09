#include <Rcpp.h>
#include <vector>
#include <numeric>
#include <algorithm>

using namespace Rcpp;


std::vector<int> argsort(std::vector<double> xs) {
  std::vector<int> result(xs.size());
  std::iota(result.begin(), result.end(), 0);
  std::sort(result.begin(), result.end(), [&] (int i, int j) {
    return xs[i] < xs[j];
  });
  
  return result;
}


//' @title A function assign ranks to data
//' @description A function assign ranks to data
//' @param xs data to be assigned with rank
//' @return Rank vector of input data 
//' @export
// [[Rcpp::export]]
NumericVector rankdata(NumericVector xs) {
  int n = xs.size();
  NumericVector result(n);
  
  std::vector<int> idx = argsort(Rcpp::as<std::vector<double>>(xs));
  
  for (int i = 0; i < n; i++) {
    result[idx[i]] = i;
  }
  
  return result;
}
