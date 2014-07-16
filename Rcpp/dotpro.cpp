#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector dotpro(NumericVector u, NumericVector v) {
	return u*v;
}
