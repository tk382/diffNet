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



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat solve_tmp_c(const arma::vec x){
  arma::mat xx = arma::zeros<arma::mat>(2,2);
  xx(0,0) = x.size();
  xx(0,1) = sum(x);
  xx(1,0) = sum(x);
  xx(1,1) = sum(x%x);
  return inv(xx);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double get_score_c(const arma::vec x,
                  const arma::vec y1,
                  const arma::vec y2) {
  //returns score statistics between vectors y1 and y2 wrt x
  int n = x.size();
  if(n != y1.size() | n != y2.size()){
    stop("x, y1, y2 should have the same length");
  }
  arma::vec w = y1+y2;
  arma::mat xx = solve_tmp_c(x);
  arma::vec intercepts = arma::ones<arma::vec>(n);
  arma::mat xt = join_cols(intercepts, x);
  double hatsigma = sum(w%w)/(n-1);
  arma::mat lhs = arma::zeros<arma::mat>(1,2);
  lhs(0,0) = sum(w%w/hatsigma - 1);
  lhs(0,1) = sum((w%w/hatsigma - 1) % x);
  arma::mat S = lhs * xx * lhs.t()/2;
  return S(0,0);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double get_degree_c(const arma::vec x,
                    const arma::vec y,
                    const arma::mat Y){
  //degree statistic for i'th gene based on matrix Y
  int K = Y.n_cols;
  double d = 0;
  for (int i=0; i < K; ++i){
    d = d + get_score_c(x, y, Y.col(i));
  }
  return d;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double get_eta_c(const double rho12,
                 const double rho23,
                 const double rho13){
  double eta = 6;
  eta += 8*rho12 + 8*rho13 + 4*rho23;
  eta += 2*rho12*rho12 + 2*rho13*rho13 + 2*rho23*rho23;;
  eta += 8*rho12*rho13 + 4*rho12*rho23 + 4*rho13*rho23;
  double denom = (2 * (2*rho12+2) * (2*rho13+2));
  eta = eta/denom - .5;
  return eta;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double dgenrayleigh_c(const double t,
                      const double alpha,
                      const double beta,
                      const bool logg = false){
  double out = 0;
  out += log(2) + log(alpha) + log(beta);
  out -= beta*t*t;
  out += (alpha - 1) * log(1-exp(-beta*t*t));
  if (logg){
    return out;
  }else{
    return exp(out);
  }
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double cgenrayleigh_c(const double t,
                      const double alpha,
                      const double beta,
                      const bool logg = false){
  double out = alpha * log(1-exp(-beta * t * t));
  if(logg){
    return out;
  }else{
    return exp(out);
  }
}
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat get_H_c(const arma::mat Sigma){
  int K = Sigma.n_rows;
  arma::mat est_H = arma::zeros<arma::mat>(K-1, K-1);
  double eta;
  for (int i=1; i < K-1; ++i){
    for (int j=(i+1); j < K; ++j){
      eta = get_eta_c(Sigma(0,i), Sigma(i,j), Sigma(j,0));
      est_H(i-1, j-1) = eta;
      est_H(j-1, i-1) = eta;
    }
  }
  est_H.diag().ones();
  return est_H;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double get_A1_c(int k,
                int p,
                arma::mat TT,
                arma::mat RR,
                int n){
  double first = 24 * (k-1) * (p-1);
  double second = -24 * n * sum(TT.diag() % RR.diag());
  double third = 0;
  for (int i=0; i < n; ++i){
    for (int j=0; j < n; ++j){
      double tmp = RR(i,i) * RR(j,j);
      double tmp2 = RR(i,j) * RR(i,j) * 2;
      third += TT(i,j) * (tmp + tmp2);
    }
  }
  third *= 6 * n;
  return first + second + third;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double get_A2_c(int k,
                int p,
                arma::mat TT,
                arma::mat RR,
                int n){
  double first = -24 * (p-1) * (p-1);
  double second = 36 * n * accu(TT.diag() % TT.diag());
  double third = -48 * accu(TT%TT);
  double fourth = 0;
  for (int i=0; i < (n-1); ++i){
    for (int j=(i+1); j < n; ++j){
      fourth += TT(i,j) * (TT(j,j) * RR(i,i) + TT(i,i) * RR(j,j));
    }
  }
  for (int i=0; i < n; ++i){
    fourth += TT(i,i) * TT(i,i) * RR(i,i);
  }
  fourth *= (-24) * n;
  // cout << first << "\n" << second << "\n" << third << "\n" << fourth << "\n";
  return (first + second + third + fourth);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double get_A3_c(int k,
                int p,
                arma::mat TT,
                arma::mat RR,
                int n){
  double first = 0;
  for (int i = 0; i < n; ++i){
    for (int j = 0; j < n; ++j){
      first += TT(i,i) * TT(j,j) * TT(i,j);
    }
  }
  first *= 24 * n;
  double second = 16 * n * accu(TT%TT%TT);
  return first + second;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List get_TT_RR_c(arma::vec A){
  int n = A.size();
  arma::mat X(n,2);
  X.col(0) = arma::ones<arma::vec>(n);
  X.col(1) = A;
  double barA = mean(A);
  double H = sum((A - barA)%(A-barA));
  arma::mat TT(n,n);
  for (int i = 0; i < (n-1); ++i){
    for (int j = (i+1); j < n; ++j){
      TT(i,j) = (A(i) - barA) * (A(j) - barA) / H;
      TT(j,i)= TT(i,j);
    }
  }
  for (int i=0; i < n; ++i){
    TT(i,i) = (A(i)- barA) * (A(i) - barA) / H;
  }
  arma::mat RR = X * (X.t() * X).i() * X.t();
  return Rcpp::List::create(Rcpp::Named("TT") = TT,
                            Rcpp::Named("RR") = RR);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double final_correction_c(arma::vec A,
                        double C,
                        int k,
                        int p){
  Rcpp::List TTRR = get_TT_RR_c(A);
  arma::mat TT = as<arma::mat>(TTRR["TT"]);
  arma::mat RR = as<arma::mat>(TTRR["RR"]);
  int n = A.size();
  double A1 = get_A1_c(k,p,TT,RR,n);
  double A2 = get_A2_c(k,p,TT,RR,n);
  double A3 = get_A3_c(k,p,TT,RR,n);

  double first = A3 * C / (p-1) / (p+1) / (p+3) * (C*C + (p+3)*C + (p+1)*(p+3));
  double second = C * (C+p+1) * (A2-3*A3) / (p-1) / (p+1);
  double third = C * (3 * A3 - 2 * A2 + A1) / (p-1);
  return ((first + second + third) / 12 / n);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double cubic_correction_c(arma::vec A,
                          double C,
                          int k,
                          int p){
  Rcpp::List TTRR = get_TT_RR_c(A);
  arma::mat TT = as<arma::mat>(TTRR["TT"]);
  arma::mat RR = as<arma::mat>(TTRR["RR"]);
  int n = A.size();
  double A1 = get_A1_c(k,p,TT,RR,n);
  double A2 = get_A2_c(k,p,TT,RR,n);
  double A3 = get_A3_c(k,p,TT,RR,n);
  double coef1 = A3/(12*n*(p-1)*(p+1)*(p+3));
  double coef2 = (A2-2*A3)/(12*n*(p-1)*(p+1));
  double coef3 = (A3-A2+A1) / (12*n*(p-1))  + 1;
  return coef1*C*C*C + coef2*C*C + coef3 * C;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec cubic_coeff_c(arma::vec A,
                          double C,
                          int k,
                          int p){
  Rcpp::List TTRR = get_TT_RR_c(A);
  arma::mat TT = as<arma::mat>(TTRR["TT"]);
  arma::mat RR = as<arma::mat>(TTRR["RR"]);
  int n = A.size();
  double A1 = get_A1_c(k,p,TT,RR,n);
  double A2 = get_A2_c(k,p,TT,RR,n);
  double A3 = get_A3_c(k,p,TT,RR,n);
  arma::vec coef = arma::zeros<arma::vec>(3);
  coef(0) = (A3-A2+A1) / (12*n*(p-1))  + 1;
  coef(1) = (A2-2*A3)/(12*n*(p-1)*(p+1));
  coef(2) = A3/(12*n*(p-1)*(p+1)*(p+3));
  return coef;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List simulate_c(arma::vec A,
                      arma::mat Sigma,
                      bool use_Hhat = true) {
  const int K = Sigma.n_rows;
  const int n = A.size();
  arma::mat Y = mvrnormArma(n, arma::zeros<arma::vec>(K), Sigma);
  arma::mat H;
  if (use_Hhat){
    arma::mat Sigmahat = cov(Y);
    H = get_H_c(Sigmahat);
  }else{
    H = get_H_c(Sigma);
  }
  arma::mat X = arma::ones<arma::mat>(n,2);
  X.col(1) = A;
  arma::vec scores(K-1);
  for (int k = 1; k < K; ++k){
    //arma::vec tmpw = Y.col(0) + Y.col(k);
    double score = get_score_c(A, Y.col(0), Y.col(k));
    scores(k-1) = score;
    // arma::vec coef = arma::solve(X, tmpw);
    // arma::vec w = tmpw - X*coef;
    // double sigmahat = sum(w%w) / n;
    // arma::vec tmp1 = w%w / sigmahat - 1;
    // arma::vec tmp2 = A % tmp1;
    // double tmp3 = sum(tmp2);
    // scores(k-1) = tmp3*tmp3 / (2*n);
  }
  return Rcpp::List::create(Rcpp::Named("Y") = Y,
                            Rcpp::Named("scores") = scores,
                            Rcpp::Named("Hhat") = H);
}
























