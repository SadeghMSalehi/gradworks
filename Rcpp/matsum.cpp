#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
double matsum(NumericMatrix m) {
	int nrow = m.nrow();
	int ncol = m.ncol();
	
	double total = 0;
	for (int i = 0; i < nrow; i++) {
		for (int j = 0; j < ncol; j++) {
			total += m(i,j);
		}
	}
	return total;
}
