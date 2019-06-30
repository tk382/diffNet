get_grad = function(sigw, sigv, uw, uv, x){
  return(sum(x * (1/sigw - (uw^2)/(sigw^2) - 1/sigv + (uv^2)/(sigv^2))))
  # var1 = 4 + sigw - sigv # 4 + 4rho
  # var2 = 4 - sigw + sigv # 4 - 4rho
  # part1 = x/var1 * (1-uw^2/var1)
  # part2 = x/var2 * (1-uv^2/var2)
  # return(sum(part1-part2))
}

const = function(sigw, sigv){
  tmp = 2 * (1/sigw^2 + 1/sigv^2)
  return(1/tmp)
  # var1 = (4 + sigw - sigv)^2
  # var2 = (4 - sigw + sigv)^2
  # denom = 1/var1 + 1/var2
  # return(2 * 1/denom)
}

post_score = function(score, coef){
  roots = polyroot(c(-score, coef))
  return(Re(roots)[abs(Im(roots)) < 1e-6][1])
}

get_q = function(y1, y2, x, coef, correction = TRUE){
  w = y1 + y2
  v = y1 - y2
  n = length(w)
  sigw = sum(w^2)/n
  sigv = sum(v^2)/n
  grad = get_grad(sigw, sigv, w, v, x)
  fisherinf = const(sigw, sigv)
  if(correction){
    q = post_score(grad^2 / length(x) * fisherinf, coef)
  }else{
    q = grad^2 / length(x) * fisherinf
  }

  return(q)
}

get_cor_r = function(rho12, rho23, rho13){
  num = (rho23 + 2 * rho12 * rho23)* (rho12^2+1) * (rho13^2 + 1)
  num = num + rho12 * rho13* (6 + 2 * rho12 + 2 * rho13 + 2 * rho23)
  num = num - rho12 * (rho13^2 + 1) * (3*rho13 + rho13 + 2 * rho12*rho23)
  num = num - rho13 * (rho12^2 + 1) * (3*rho12 + rho12 + 2 * rho13 * rho23)
  # num = rho12^3*rho13^3 - 3*rho12^2*rho13^2*rho23
  # num = num + 2*rho12*rho13*rho23^2 - rho12^2*rho23
  # num = num - rho13^2*rho23 - rho12 * rho13 + rho23
  # num = num + rho12^3 * rho13 + rho13^3 * rho12
  denom = (1-rho12^2)*(1-rho13^2) * sqrt(1+rho12^2)*sqrt(1+rho13^2)
  return(num / denom)
}

get_est_H = function(Sigma){
  K = nrow(Sigma)
  est_H = matrix(0, K-1, K-1)
  diag(est_H) = 1
  for (i in 1:(K-1)){
    for (j in (i+1):K){
      eta = get_cor_r(Sigma[1, i], Sigma[i,j], Sigma[j,1])
      est_H[i-1, j-1] = est_H[j-1, i-1] = eta
    }
  }
  return(est_H)
}



#three h functions to define $\sigma^2$.
#$\sigma^2 = 2 + 2\rho$, ranges from 2 to 4 (assumes rho is positive)


#' Returns rho
#'
#' @param X scaled covariate
#' @param alpha length 2 vector for intercept + slope for X
#'
#' @export
h1_fisher = function(X, alpha){
  tmp = X %*% alpha
  return((exp(tmp)-1)/(exp(tmp)+1))
}

#' Returns rho
#'
#' @param X scaled covariate
#' @param alpha length 2 vector for intercept + slope for X
#'
#' @export
h2_sqrt = function(X, alpha){
  tmp = (X %*% alpha) / sqrt(1 + (X %*% alpha)^2)
  return(tmp*2-1)
}


#' Returns rho
#'
#' @param X scaled covariate
#' @param alpha length 2 vector for intercept + slope for X
#'
#' @export
h3_cdf = function(X, alpha){
  tmp = pnorm(X %*%alpha, 0, 10)*2 - 1
  return(tmp)
}

#' Returns rho
#'
#' @param X scaled covariate
#' @param alpha length 2 vector for intercept + slope for X
#'
#' @export
h4_sin = function(X, alpha){
  tmp = sin(X %*% alpha * 2)
  return(tmp)
}

#' Returns rho
#'
#' @param X scaled covariate
#' @param alpha length 2 vector for intercept + slope for X
#'
#' @export
h5_gumbel = function(X, alpha){
  tmp = pgumbel(X %*% alpha,  loc = 1, scale=2)
  return(tmp)
}


#' Returns rho
#'
#' @param X scaled covariate
#' @param alpha length 2 vector for intercept + slope for X
#'
#' @export
h6_quadratic = function(X, alpha){
  sigma1 = (X %*% alpha-0.1)^2 - 0.99
  # sigma1 = pmin(sigma1, 0.9)
  # sigma1 = pmax(sigma1, -0.9)
  return(sigma1)
}


# powercheck = function(out, B, n){
#   size = as.character(n)
#   df = data.frame(n = c("LM", "LA", "Fisher", "Normal CDF", "sin"),
#                   size = c(sum(out$lm < 0.05)/B,
#                            sum(out$la < 0.05)/B,
#                            sum(out$fisher < 0.05)/B,
#                            sum(out$cdf < 0.05)/B,
#                            sum(out$sin < 0.05)/B))
#   colnames(df)[2] = size
#   return(df)
# }
powercheck = function(out, B, n){
  size = as.character(n)
  df = data.frame(n = c("LM", "LA", "Fisher"),
                  size = c(sum(out$lm < 0.05)/B,
                           sum(out$la < 0.05)/B,
                           sum(out$fisher < 0.05)/B))
  colnames(df)[2] = size
  return(df)
}
