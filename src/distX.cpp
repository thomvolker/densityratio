#include <Rcpp.h>
using namespace Rcpp;

//' Create a Gram matrix with squared Euclidean distances between all
//' observations in the input matrix
//'
//' @param x A numeric input matrix
//' @param nr Number of rows of the input matrix
//' @param nc Number of columns of the input matrix

 // [[Rcpp::export]]
NumericVector distX(NumericMatrix x, int nr, int nc) {
  double dev, dist;
  NumericVector out(nr*(nr-1)/2);
  int count = 0;

  for(int i = 0; i < nr; i++) {
    for(int j = i+1; j < nr; j++) {
      dist = 0;
      for (int c = 0; c < nc; c++) {
        dev = x(i,c) - x(j,c);
        dist += dev*dev;
        }
      out[count] = dist;
      count++;
      }
    }
  return out;
 }
