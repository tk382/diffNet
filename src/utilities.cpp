#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace std;

//const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double get_grad_c(const double sigw,
                      const double sigv,
                      const arma::vec uw,
                      const arma::vec uv,
                      const arma::vec x){
  arma::vec tmp1 = uw%uw / (sigw*  sigw);
  arma::vec tmp2 = uv%uv / (sigv * sigv);
  arma::vec tmp3 = arma::zeros<arma::vec>(uw.size());
  tmp3.fill(1/sigw);
  arma::vec tmp4 = arma::zeros<arma::vec>(uw.size());
  tmp4.fill(1/sigv);
  return sum(x % (tmp3 - tmp1 - tmp4 + tmp2));
}

// [[Rcpp::export]]
double const_c(const double sigw,
        const double sigv){
  double tmp = 2 * (1/(sigw*sigw) + 1/(sigv*sigv));
  return 1/tmp;
}

//' get_score_c
//' The order of y1 and y2 does not matter.
//'
//' @param y1 observation of the first variable
//' @param y2 observation of the second variable
//' @param x covariate vector
//' @param coef coefficients for the cubic function. see cubic_coeff_c
//'
//' @export
// [[Rcpp::export]]
double get_score_c(const arma::vec x,
            const arma::vec y1,
            const arma::vec y2){
  arma::vec w = y1 + y2;
  arma::vec v = y1 - y2;
  int n = w.size();
  double sigw = sum(w%w)/n;
  double sigv = sum(v%v)/n;
  double grad = get_grad_c(sigw, sigv, w, v, x);
  double fisherinf = const_c(sigw, sigv);
  double q = grad * grad / n * fisherinf;
  return q;
}


//' mvrnormArma
//'
//' This function returns n samples of multivariate normal distribution with mean mu and variance Sigma.
//'
//' @param n Sample size
//' @param mu Mean vector
//' @param Sigma Variance-covariance matrix
//' @export
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


// [[Rcpp::export]]
arma::mat solve_tmp_c(const arma::vec x){
  arma::mat xx = arma::zeros<arma::mat>(2,2);
  xx(0,0) = x.size();
  xx(0,1) = sum(x);
  xx(1,0) = sum(x);
  xx(1,1) = sum(x%x);
  return inv(xx);
}

//' get_score_w_c
//'
//' This function takes the sum of two variables, w, the difference of them, v, and one dimensional covariate x
//' and returns the score test statistic
//'
//' @param x Covariate vector
//' @param w sum of variable1 and variable2
//' @export
// [[Rcpp::export]]
double get_score_w_c(const arma::vec x,
                   const arma::vec w,
                   const arma::vec v) {
  int n = w.size();
  double sigw = sum(w%w)/n;
  double sigv = sum(v%v)/n;
  double grad = get_grad_c(sigw, sigv, w, v, x);
  double fisherinf = const_c(sigw, sigv);
  double q = grad * grad / n * fisherinf;
  return q;
}

//' get_score_W_c
//'
//' This function returns a vector of score statistics from one dimensional covariate x and coexpression matrix W
//' W should be in the form of (y1+y2, y1+y3, ... y1+yK)
//' The function returns K-1 scores for each pair
//'
//' @param x Covariate vector
//' @param W Sums for each pairs of variables
//' @param V Differences for each pairs of variables
//' @export
// [[Rcpp::export]]
arma::vec get_score_W_c(const arma::vec x,
                     const arma::mat W,
                     const arma::mat V) {
  int n = x.size();
  int K = W.n_cols;
  arma::vec out = arma::zeros<arma::vec>(K);
  for (int i=0; i < K; ++i){
    out(i) = get_score_w_c(x, W.col(i), V.col(i));
  }
  return out;
}



//' get_degree_c
//'
//' This function returns a sum statistic d for vector y tested with all other variables in matrix Y against covariate x
//'
//' @param x Covariate vector
//' @param y Variable of interest
//' @param Y All other variables to test with y
//' @export
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

//' get_degree_c
//'
//' This function returns a sum statistic d for vector y tested with all other variables in matrix Y against covariate x
//'
//' @param x Covariate vector
//' @param W Matrix of [y1+y2, y1+y3, .., y1+yK]
//' @param V Matrix of [y1-y2, y1-y3, ..., y1-yK]
//' @export
// [[Rcpp::export]]
double get_degree_w_c(const arma::vec x,
                      const arma::mat W,
                      const arma::mat V){
  //degree statistic for i'th gene based on matrix Y
  int K = W.n_cols;
  double d = 0;
  for (int i=0; i < K; ++i){
    d = d + get_score_w_c(x, W.col(i), V.col(i));
  }
  return d;
}

//' @export
// [[Rcpp::export]]
double get_eta_c(const double rho12,
                 const double rho23,
                 const double rho13){
  double num = (rho23 + 2 * rho12 * rho23)* (rho12*rho12+1) * (rho13*rho13 + 1);
  num = num + rho12 * rho13* (6 + 2 * rho12 + 2 * rho13 + 2 * rho23);
  num = num - rho12 * (rho13*rho13 + 1) * (3*rho13 + rho13 + 2 * rho12*rho23);
  num = num - rho13 * (rho12*rho12 + 1) * (3*rho12 + rho12 + 2 * rho13 * rho23);
  double denom = (1-rho12*rho12)*(1-rho13*rho13) * sqrt(1+rho12*rho12)*sqrt(1+rho13*rho13);
  return num/denom;
}

//' @export
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

//' @export
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


//' @export
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


//' @export
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


//' @export
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


//' @export
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


//' @export
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



//' small_sample_correction_c
//'
//' This function returns a sum statistic d for vector y tested with all other variables in matrix Y against covariate x
//'
//' @param x Covariate vector
//' @param C chi-squared threshold for desired FDR rate
//' @param k number of covariates + 1
//' @param p number of covariates + 1
//' @export
// [[Rcpp::export]]
double small_sample_correction_c(arma::vec x,
                        double C,
                        int k,
                        int p=2){
  Rcpp::List TTRR = get_TT_RR_c(x);
  arma::mat TT = as<arma::mat>(TTRR["TT"]);
  arma::mat RR = as<arma::mat>(TTRR["RR"]);
  int n = x.size();
  double A1 = get_A1_c(k,p,TT,RR,n);
  double A2 = get_A2_c(k,p,TT,RR,n);
  double A3 = get_A3_c(k,p,TT,RR,n);

  double first = A3 * C / (p-1) / (p+1) / (p+3) * (C*C + (p+3)*C + (p+1)*(p+3));
  double second = C * (C+p+1) * (A2-3*A3) / (p-1) / (p+1);
  double third = C * (3 * A3 - 2 * A2 + A1) / (p-1);
  return ((first + second + third) / 12 / n);
}



//' @export
// [[Rcpp::export]]
double cubic_correction_c(arma::vec x,
                          double C,
                          int k,
                          int p){
  Rcpp::List TTRR = get_TT_RR_c(x);
  arma::mat TT = as<arma::mat>(TTRR["TT"]);
  arma::mat RR = as<arma::mat>(TTRR["RR"]);
  int n = x.size();
  double A1 = get_A1_c(k,p,TT,RR,n);
  double A2 = get_A2_c(k,p,TT,RR,n);
  double A3 = get_A3_c(k,p,TT,RR,n);
  double coef1 = A3/(12*n*(p-1)*(p+1)*(p+3));
  double coef2 = (A2-2*A3)/(12*n*(p-1)*(p+1));
  double coef3 = (A3-A2+A1) / (12*n*(p-1))  + 1;
  return coef1*C*C*C + coef2*C*C + coef3*C;
}


//' @export
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

//' simulate_c
//'
//' Given fixed covariate and covariance matrix, this function simulates multivariate normal Y
//' and returns the score vector for each variable
//'
//' @param x Covariate vector
//' @param Sigma true covariance matrix of matrix Y
//' @param use_Hhat Boolean argument whether to use true Sigma or estimate Sigma
//' @export
// [[Rcpp::export]]
Rcpp::List simulate_c(arma::vec x,
                      arma::mat Sigma,
                      bool use_Hhat = true) {
  const int K = Sigma.n_rows;
  const int n = x.size();
  arma::mat Y = mvrnormArma(n, arma::zeros<arma::vec>(K), Sigma);
  arma::mat H;
  if (use_Hhat){
    arma::mat Sigmahat = cov(Y);
    H = get_H_c(Sigmahat);
  }else{
    H = get_H_c(Sigma);
  }
  arma::mat X = arma::ones<arma::mat>(n,2);
  X.col(1) = x;
  arma::vec scores(K-1);
  for (int k = 1; k < K; ++k){
    double score = get_score_c(x, Y.col(0), Y.col(k));
    scores(k-1) = score;
  }
  return Rcpp::List::create(Rcpp::Named("Y") = Y,
                            Rcpp::Named("scores") = scores,
                            Rcpp::Named("Hhat") = H);
}



//' @export
// [[Rcpp::export]]
arma::mat shuffle_x_c(arma::vec x, int B){
  arma::mat X(x.size(), B);
  for (int b=0; b < B; ++b){
    X.col(b) = shuffle(x);
  }
  return X;
}

//' @export
// [[Rcpp::export]]
arma::mat store_W_c(const arma::vec y,
                    const arma::mat smallY){
  int n = smallY.n_rows;
  int K = smallY.n_cols;
  arma::mat W(n, K);
  for (int k = 0; k < K; ++k){
    W.col(k) = y + smallY.col(k);
  }
  return W;
}

//' @export
// [[Rcpp::export]]
arma::mat store_V_c(const arma::vec y,
                    const arma::mat smallY){
  int n = smallY.n_rows;
  int K = smallY.n_cols;
  arma::mat V(n, K);
  for (int k = 0; k < K; ++k){
    V.col(k) = y - smallY.col(k);
  }
  return V;
}


//' bootstrap_c
//'
//' This function performs permutation test for given covariate x B times
//' Also takes as input W instead of Y matrix
//'
//' @param x Covariate vector
//' @param B number of permutations
//' @param W transformed variable matrix
//' @export
// [[Rcpp::export]]
arma::mat bootstrap_c(const arma::vec x,
                      const int B,
                      const arma::mat W,
                      const arma::mat V){
  arma::mat Xb = shuffle_x_c(x, B);
  const int K = W.n_cols+1;
  arma::mat out(B, K-1);
  for (int b = 0; b < B; ++b){
    for (int k = 0; k < (K-1); ++k){
      out(b, k) = get_score_w_c(Xb.col(b), W.col(k), V.col(k));
    }
  }
  return out;
}










