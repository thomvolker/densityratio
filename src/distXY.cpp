#include <Rcpp.h>
using namespace Rcpp;

//' Create a Gram matrix with squared Euclidean distances between
//' observations in the input matrix x and the input matrix y
//' @param x A numeric input matrix
//' @param y A numeric input matrix with the same variables as x
//' @param nrx Number of rows of the input matrix x
//' @param nry Number of rows of the input matrix y
//' @param nc Number of columns of the input matrices x and y
//' @export
// [[Rcpp::export]]

NumericMatrix distXY(NumericMatrix x, NumericMatrix y, int nrx, int nry, int nc) {
  double dev, dist;
  NumericMatrix out(nrx, nry);

  for(int i = 0; i < nrx; i++) {
    for(int j = 0; j < nry; j++) {
      dist = 0;
      for (int c = 0; c < nc; c++) {
        dev = x(i,c) - y(j,c);
        dist += dev*dev;
      }
      out(i,j) = dist;
    }
  }
  return out;
}
