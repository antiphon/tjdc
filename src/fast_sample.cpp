#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int sample_fast_c(NumericVector prob) {
  double u = Rcpp::as<double>(Rcpp::runif(1));
  double s = 0;
  for( int i = 0; i < prob.size(); i++){
    s  += prob(i);
    if(u < s){ return i+1;}
  }
  return prob.size();
}
