pairwise_test = function(L1, L2, y1, y2){
  #returns the chi square statistics and p-value
  #for the effects of L1 and L2 on the correlation between y1 and y2
  if(length(L1)!=length(y1) || length(L2)!=length(y2)){
    stop('local ancestry and expression level should have equal sample size')
    return(list(L1=c(NA, NA), L2=c(NA, NA)))
  }
  #First, effects of L1 on cor(y1,y2)
  ind0 = which(L1==0); ind1 = which(L1==1); ind2 = which(L1==2)
  dat = as.data.frame(cbind(y1,y2,L1,L2))
  dat0 = dat[ind0, ]; n0 = nrow(dat0)
  dat1 = dat[ind1, ]; n1 = nrow(dat1)
  dat2 = dat[ind2, ]; n2 = nrow(dat2)
  if(n0<2 || n1<2 || n2<2){
    print('not enough sample sizes for each L1, assigning NA')
  }
  rho0 = cor(dat0$y1, dat0$y2); z0 = fishertrans(rho0); sig0 = 1/(n0-3);
  rho1 = cor(dat1$y1, dat1$y2); z1 = fishertrans(rho1); sig1 = 1/(n1-3);
  rho2 = cor(dat2$y1, dat2$y2); z2 = fishertrans(rho2); sig2 = 1/(n2-3);
  mu_nu    = rep(mu_null(z0,z1,z2,sig0,sig1,sig2), 3)
  mubeta   = mubeta_alt(z0,z1,z2,sig0,sig1,sig2)
  mu_alt1   = mubeta$mu
  beta_alt1 = mubeta$beta
  mubeta   = c(mu_alt1, mu_alt1+beta_alt1, mu_alt1+2*beta_alt1)
  Sigma    = diag(c(sig0, sig1, sig2))
  l_alt    = dmvnorm(c(z0,z1,z2), mubeta, Sigma, log=TRUE)
  l_null   = dmvnorm(c(z0,z1,z2), mu_nu, Sigma, log=TRUE)
  chistat1 = 2*(l_alt-l_null)
  pval1    = pchisq(2*(l_alt - l_null), 1, lower.tail = FALSE)

  #Second, effects of L2 on cor(y1,y2)
  ind0 = which(L2==0); ind1 = which(L2==1); ind2 = which(L2==2)
  dat0 = dat[ind0, ]; n0 = nrow(dat0)
  dat1 = dat[ind1, ]; n1 = nrow(dat1)
  dat2 = dat[ind2, ]; n2 = nrow(dat2)
  if(n0<2 || n1<2 || n2<2){
    print('not enough sample sizes for each L2, assigning NA')
  }
  rho0 = cor(dat0$y1, dat0$y2); z0 = fishertrans(rho0); sig0 = 1/(n0-3);
  rho1 = cor(dat1$y1, dat1$y2); z1 = fishertrans(rho1); sig1 = 1/(n1-3);
  rho2 = cor(dat2$y1, dat2$y2); z2 = fishertrans(rho2); sig2 = 1/(n2-3);

  mu_nu    = rep(mu_null(z0,z1,z2,sig0,sig1,sig2), 3)
  mubeta   = mubeta_alt(z0,z1,z2,sig0,sig1,sig2)
  mu_alt2   = mubeta$mu
  beta_alt2 = mubeta$beta
  mubeta   = c(mu_alt2, mu_alt2+beta_alt2, mu_alt2+2*beta_alt2)
  Sigma    = diag(c(sig0, sig1, sig2))
  l_alt    = dmvnorm(c(z0,z1,z2), mubeta, Sigma, log=TRUE)
  l_null   = dmvnorm(c(z0,z1,z2), mu_nu, Sigma, log=TRUE)

  chistat2 = 2*(l_alt-l_null)
  pval2    = pchisq(2*(l_alt - l_null), 1, lower.tail = FALSE)

  return(list(L1=c(chistat1,pval1,mu_alt1,beta_alt1),
              L2=c(chistat2,pval2,mu_alt2,beta_alt2)))
}
