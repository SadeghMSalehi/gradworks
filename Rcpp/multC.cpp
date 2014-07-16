#include <Rcpp.h>
#include <iostream>


using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix multC(NumericMatrix m, NumericVector v) {
	int ncols = m.cols();
	int nrows = m.rows();

	NumericMatrix out(nrows,1);
	for (int j = 0; j < nrows; j++) {
		double row_total = 0;
		for (int k = 0; k < ncols; k++) {
			row_total += (m(j,k) * v(k));
		}
		out(j,0) = row_total;
	}
	return out;
}
