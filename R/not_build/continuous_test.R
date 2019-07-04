continuous_test = function(x, w){
  n = length(x)
  mubeta    = continuous_mubeta(x,w)$x
  mu1       = continuous_mu(x,w)
  mu2       = mubeta[1]; beta2 = mubeta[2]
  thetahat1 = exp(mu1);   rhohat1 = (thetahat1-1)/(thetahat1+1)
  thetahat2 = exp(mu2 + beta2 * x); rhohat2 = (thetahat2-1)/(thetahat2+1)
  lik1   = sum(dnorm(w, 0, sqrt(rhohat1+1), log=TRUE))
  lik2   = 0
  for (i in 1:n){
    lik2 = lik2 + dnorm(w[i], 0, sqrt(rhohat2[i]+1), log=TRUE)
  }
  chistat   = 2*(lik2-lik1)
  p         = pchisq(chistat, 1, lower.tail = FALSE)
  return(list(chi = chistat, p = p, mu1 = mu1, mu2 = mu2, beta2 = beta2))
}
