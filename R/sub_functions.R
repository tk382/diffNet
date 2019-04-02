fishertrans = function(rho, n){
  log((1+rho)/(1-rho))/2
}

mu_null = function(z0,z1,z2,sig0,sig1,sig2){
  Phi = sig0*sig1+sig1*sig2+sig2*sig0
  Lambda = sig1*sig2*z0 + sig0*sig2*z1 + sig0*sig1*z2
  return(Lambda/Phi)
}

mubeta_alt = function(z0,z1,z2,sig0,sig1,sig2){
  #sig0, sig1, sig2 are actually sig0^2, sig1^2, sig^2
  Phi = sig0*sig1+sig1*sig2+sig2*sig0
  Lambda = sig1*sig2*z0 + sig0*sig2*z1 + sig0*sig1*z2
  beta = (sig2+2*sig1)*Lambda - (sig2*z1+2*sig1*z2)*Phi
  beta = beta / (sig2+2*sig1) / (sig0*sig2+2*sig0*sig1 - 2*Phi)
  mu = sig1*sig2*z0+sig2*sig0*(z1-beta) + sig0*sig1*(z2-2*beta)
  mu = mu/Phi
  return(list(mu=mu, beta=beta))
}
