
#' Returns the likelihood assuming fisher
#'
#' @param alpha length 2 vector for intercept + slope for X
#' @param opt covariate vector X, sum of two vectors Y
#'
#' @export
likelihood_fisher = function(alpha, opt = list(X, U)){
  sigma2 = h1_fisher(X, alpha)
  l = sum(dnorm(U, 0, sqrt(sigma2), log = TRUE))
  return(l)
}

#' Returns the likelihood assuming linear
#'
#'
#' @param alpha length 2 vector for intercept + slope for X
#' @param opt covariate vector X, sum of two vectors Y
#'
#' @export
likelihood_sqrt = function(alpha, opt = list(X, U)){
  sigma2 = h2_sqrt(X, alpha)
  l = sum(dnorm(U, 0, sqrt(sigma2), log = TRUE))
  return(l)
}

#' Returns the likelihood assuming gaussian cdf function
#'
#' @param alpha length 2 vector for intercept + slope for X
#' @param opt covariate vector X, sum of two vectors Y
#'
#' @export
likelihood_cdf = function(alpha, opt = list(X, U)){
  sigma2 = h3_cdf(X, alpha)
  l = sum(dnorm(U, 0, sqrt(sigma2), log=TRUE))
  return(l)
}

#' Returns the likelihood assuming sin function
#'
#' @param alpha length 2 vector for intercept + slope for X
#' @param opt covariate vector X, sum of two vectors Y
#'
#' @export
likelihood_sin = function(alpha, opt = list(X, U)){
  sigma2 = h4_sin(X, alpha)
  l = sum(dnorm(U, 0, sqrt(sigma2), log=TRUE))
  return(l)
}

#' Returns the likelihood assuming cdf of exp function
#'
#' @param alpha length 2 vector for intercept + slope for X
#' @param opt covariate vector X, sum of two vectors Y
#'
#' @export
likelihood_gumbel = function(alpha, opt = list(X, U)){
  sigma2 = h5_gumbel(X, alpha)
  l = sum(dnorm(U, 0, sqrt(sigma2), log=TRUE))
  return(l)
}

#' Returns the likelihood assuming cdf of exp function
#'
#' @param alpha length 2 vector for intercept + slope forl X
#' @param opt covariate vector X, sum of two vectors Y
#'
#' @export
likelihood_quadratic = function(alpha, opt = list(X, U)){
  sigma2 = h6_quadratic(X, alpha)
  l = sum(dnorm(U, 0, sqrt(sigma2), log=TRUE))
  return(l)
}

#' Returns the likelihood assuming null
#'
#' @param alpha length 2 vector for intercept + slope for X
#' @param opt covariate vector X, sum of two vectors Y
#'
#' @export
likelihood_null = function(rho, opt = list(X, U)){
  l = sum(dnorm(U, 0, sqrt(2+2*rho), log=TRUE))
  return(l)
}



#' Returns the likelihood assuming null
#'
#' @param n sample size
#' @param B number of simulations
#'
#' @export
null_simulation = function(n, B){

  likelihood_fisher = function(alpha, opt = list(X, U)){
    sigma2 = h1_fisher(X, alpha)
    l = sum(dnorm(U, 0, sqrt(sigma2), log = TRUE))
    return(l)
  }

  likelihood_sqrt = function(alpha, opt = list(X, U)){
    sigma2 = h2_sqrt(X, alpha)
    l = sum(dnorm(U, 0, sqrt(sigma2), log = TRUE))
    return(l)
  }

  likelihood_cdf = function(alpha, opt = list(X, U)){
    sigma2 = h3_cdf(X, alpha)
    l = sum(dnorm(U, 0, sqrt(sigma2), log=TRUE))
    return(l)
  }

  likelihood_sin = function(alpha, opt = list(X, U)){
    sigma2 = h4_sin(X, alpha)
    l = sum(dnorm(U, 0, sqrt(sigma2), log=TRUE))
    return(l)
  }

  likelihood_gumbel = function(alpha, opt = list(X, U)){
    sigma2 = h5_gumbel(X, alpha)
    l = sum(dnorm(U, 0, sqrt(sigma2), log=TRUE))
    return(l)
  }

  likelihood_qudratic = function(alpha, opt = list(X,U)){
    sigma2 = h6_quadratic(X, alpha)
    l = sum(dnorm(U, 0, sqrt(sigma2), log=TRUE))
    return(l)
  }

  likelihood_null = function(sigma2, opt = list(X, U)){
    l = sum(dnorm(U, 0, sqrt(sigma2), log=TRUE))
    return(l)
  }


  likelihood_quadratic = function(alpha, opt = list(X, U)){
    sigma2 = h6_quadratic(X, alpha)
    l = sum(dnorm(U, 0, sqrt(sigma2), log=TRUE))
    return(l)
  }


  X = scale(runif(n)) * sqrt(n) / sqrt(n-1)
  X = cbind(rep(1, n), X)
  coef = cubic_coeff_c(X[,2], pchisq(0.95, 1), 2, 2)

  p1 = p2 = p3 = p4 = p5 = p6 = bp = rep(NA, B)
  for (i in 1:B){
    U = mvrnorm(n, rep(0,2), matrix(c(1, 0.5, 0.5, 1), nrow=2))
    U = scale(U) * sqrt(n) / sqrt(n-1)
    U = rowSums(U)
    mle0 = find.mle(likelihood_null,
                    1,
                    method = "optim",
                    opt = list(X = X, U = U),
                    lower=1e-6, upper=100)$par
    mle1 = find.mle(likelihood_fisher,
                    c(0,0),
                    method="optim",
                    opt = list(X = X, U = U),
                    lower=-5, upper=5)$par
    mle2 = find.mle(likelihood_cdf,
                    c(0,0),
                    method="optim",
                    opt = list(X = X, U = U),
                    lower=-5,upper=5)$par
    mle3 = find.mle(likelihood_gumbel,
                    c(0,0),
                    method="optim",
                    opt = list(X = X, U = U),
                    lower=-5, upper=5)$par
    mle4 = find.mle(likelihood_sin,
                    c(0,0),
                    method="optim",
                    opt=list(X=X, U=U),
                    lower=-5, upper=5)$par
    mle5 = find.mle(likelihood_sqrt,
                    c(0,0),
                    method="optim",
                    opt=list(X=X, U=U),
                    lower=-5, upper=5)$par
    mle6 = find.mle(likelihood_quadratic,
                    c(1,1),
                    method = "optim",
                    opt = list(X=X, U=U),
                    lower = -3.9, upper=3.9)$par

    l0 = likelihood_null(mle0, opt=list(X, U))
    l1 = likelihood_fisher(mle1, opt=list(X=X, U=U))
    l2 = likelihood_cdf(mle2, opt=list(X=X, U=U))
    l3 = likelihood_gumbel(mle3, opt=list(X=X, U=U))
    l4 = likelihood_sin(mle4, opt=list(X=X, U=U))
    l5 = likelihood_sqrt(mle5, opt=list(X=X, U=U))
    l6 = likelihood_quadratic(mle6, opt = list(X=X, U=U))

    llr1 = l0 - l1
    llr2 = l0 - l2
    llr3 = l0 - l3
    llr4 = l0 - l4
    llr5 = l0 - l5
    llr6 = l0 - l6

    p1[i] = 1-pchisq(-2 * llr1, 1)
    p2[i] = 1-pchisq(-2 * llr2, 1)
    p3[i] = 1-pchisq(-2 * llr3, 1)
    p4[i] = 1-pchisq(-2 * llr4, 1)
    p5[i] = 1-pchisq(-2 * llr5, 1)
    p6[i] = 1-pchisq(-2 * llr6, 1)

    s = get_score_w_c(x = X[,2], w = U)
    s = post_score(s, coef)
    bp[i] = 1-pchisq(s, 1)
  }
  return(list(fisher = p1, cdf = p2, gumbel = p3, sin = p4, sqrt = p5, quadratic = p6, bp = bp))
}

#' Returns the likelihood assuming null
#'
#' @param n sample size
#' @param B number of simulations
#' @param alpha effect size
#' @param type data generating method
#'
#' @export
alt_simulation = function(n, B, alphas = seq(-1.2, 1.2, length=7), type = c("fisher", "sqrt", "cdf", "gumbel", "sin")){
  likelihood_fisher = function(alpha, opt = list(X, U)){
    sigma2 = h1_fisher(X, alpha)
    l = sum(dnorm(U, 0, sqrt(sigma2), log = TRUE))
    return(l)
  }

  likelihood_sqrt = function(alpha, opt = list(X, U)){
    sigma2 = h2_sqrt(X, alpha)
    l = sum(dnorm(U, 0, sqrt(sigma2), log = TRUE))
    return(l)
  }

  likelihood_cdf = function(alpha, opt = list(X, U)){
    sigma2 = h3_cdf(X, alpha)
    l = sum(dnorm(U, 0, sqrt(sigma2), log=TRUE))
    return(l)
  }

  likelihood_sin = function(alpha, opt = list(X, U)){
    sigma2 = h4_sin(X, alpha)
    l = sum(dnorm(U, 0, sqrt(sigma2), log=TRUE))
    return(l)
  }

  likelihood_gumbel = function(alpha, opt = list(X, U)){
    sigma2 = h5_gumbel(X, alpha)
    l = sum(dnorm(U, 0, sqrt(sigma2), log=TRUE))
    return(l)
  }

  likelihood_null = function(sigma2, opt = list(X, U)){
    l = sum(dnorm(U, 0, sqrt(sigma2), log=TRUE))
    return(l)
  }

  likelihood_quadratic = function(alpha, opt = list(X, U)){
    sigma2 = h6_quadratic(X, alpha)
    l = sum(dnorm(U, 0, sqrt(sigma2), log=TRUE))
    return(l)
  }

  p1 = p2 = p3 = p4 = p5 = p6 = bp = matrix(NA, B, length(alphas))
  X = scale(runif(n)) * sqrt(n) / sqrt(n-1)
  X = cbind(rep(1, n), X)
  for (i in 1:length(alphas)){
    a = alphas[i]
    coef = cubic_coeff_c(X[,2], pchisq(0.95, 1), 2, 2)
    if (type=="fisher"){ sigma1 = h1_fisher(X, c(0.1, a)) }
    if (type=="sqrt"){ sigma1 = h2_sqrt(X, c(0.1, a))}
    if (type=="cdf"){ sigma1 = h3_cdf(X, c(0.1, a))}
    if (type==("sin")){ sigma1 = h4_sin(X, c(0.1, a))}
    if (type==("gumbel")){ sigma1 = h5_gumbel(X, c(0.1, a))}
    if (type=="quadratic"){sigma1 = h6_quadratic(X, c(0.1, a))}
    rho = (sigma1-2)/2
    for (b in 1:B){
      #genearte Y
      Y = matrix(NA, n, 2)
      for (c in 1:n){
        Sigma = matrix(rho[c], 2, 2); diag(Sigma) = 1
        Y[c,] = mvrnormArma(1, rep(0,2), Sigma)
      }
      U = rowSums(Y)

      mle0 = find.mle(likelihood_null,
                      var(U),
                      method="optim",
                      opt = list(X = X, U = U))$par
      mle1 = find.mle(likelihood_fisher,
                      c(var(U),0),
                      method="optim",
                      opt = list(X = X, U = U),
                      lower=-5, upper=5)$par
      mle2 = find.mle(likelihood_cdf,
                      c(var(U),0),
                      method="optim",
                      opt = list(X = X, U = U),
                      lower=-5, upper=5)$par
      mle3 = find.mle(likelihood_gumbel,
                      c(var(U),0),
                      method="optim",
                      opt = list(X = X, U = U),
                      lower=-5,
                      upper=5)$par
      mle4 = find.mle(likelihood_sin,
                      c(var(U), 0),
                      method = "optim",
                      opt = list(X=X, U=U),
                      lower=-5, upper=5)$par
      mle5 = find.mle(likelihood_sqrt,
                      c(var(U), 0),
                      method="optim",
                      opt=list(X=X, U=U),
                      lower=-5, upper=5)$par
      mle6 = find.mle(likelihood_quadratic,
                      c(var(U), 0),
                      method = "optim",
                      opt=list(X=X, U=U),
                      lower = -5, upper=5)$par


      l0 = likelihood_null(mle0)
      l1 = likelihood_fisher(mle1)
      l2 = likelihood_cdf(mle2)
      l3 = likelihood_gumbel(mle3)
      l4 = likelihood_sin(mle4)
      l5 = likelihood_sqrt(mle5)
      l6 = likelihood_sqrt(mle6)

      llr1 = l0 - l1
      llr2 = l0 - l2
      llr3 = l0 - l3
      llr4 = l0 - l4
      llr5 = l0 - l5
      llr6 = l0 - l6

      p1[b,i] = 1-pchisq(-2 * llr1, 1)
      p2[b,i] = 1-pchisq(-2 * llr2, 1)
      p3[b,i] = 1-pchisq(-2 * llr3, 1)
      p4[b,i] = 1-pchisq(-2 * llr4, 1)
      p5[b,i] = 1-pchisq(-2 * llr5, 1)
      p6[b,i] = 1-pchisq(-2 * llr6, 1)

      s = get_score_w_c(x = X[,2], w = U)
      s = post_score(s, coef)
      bp[b,i] = 1-pchisq(s, 1)
    }
  }

  return(list(fisher = p1, cdf = p2, gumbel = p3, sin = p4, sqrt = p5, qudratic = p6, bp = bp))
}


power = function(out, alphas){
  B = nrow(out$fisher)
  err_fisher = err_cdf = err_gumbel = err_sin = err_sqrt = err_bp = rep(NA, length(alphas))
  for (i in 1:length(alphas)){
    err_fisher[i] = sum(out$fisher[1:B,i] < 0.05) / B
    err_cdf[i] = sum(out$cdf[1:B,i] < 0.05) / B
    err_gumbel[i] = sum(out$gumbel[1:B,i] < 0.05) / B
    err_sin[i] = sum(out$sin[1:B,i] < 0.05) / B
    err_sqrt[i] = sum(out$sqrt[1:B,i] < 0.05) / B
    err_bp[i] = sum(out$bp[1:B,i] < 0.05) / B
  }

  power_df = data.frame(alpha = alphas,
                        score_test = err_bp,
                        fisher = err_fisher,
                        normal.cdf = err_cdf,
                        gumbel.cdf = err_gumbel,
                        sin.cdf = err_sin,
                        sqrt = err_sqrt)
  return(power_df)
}


