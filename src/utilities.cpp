#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace std;

//const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat mvrnormArma(const int n,
                      const arma::vec mu,
                      const arma::mat Sigma) {
  //returns random multivariate normal vectors with mean mu and covariance Sigma
  //input : integer n for the number of vectors you'd like to draw
  //      : vector mu for the mean
  //      : matrix Sigma for the covariance - needs to be psd
  int ncols = Sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(Sigma);
}