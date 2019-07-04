pairwise_test_AAEA = function(A, y1, y2){

  ind0 = which(A==0); ind1 = which(A==1);
  dat = as.data.frame(cbind(y1,y2))
  dat0 = dat[ind0, ]; n0 = nrow(dat0)
  dat1 = dat[ind1, ]; n1 = nrow(dat1)

  rho0 = cor(dat0$y1, dat0$y2); z0 = fishertrans(rho0); sig0 = 1/(n0-3);
  rho1 = cor(dat1$y1, dat1$y2); z1 = fishertrans(rho1); sig1 = 1/(n1-3);

  z_score = (z1-z0) / sqrt(1/(n1-3) + 1/(n0-3))
  p = 2*pnorm(-abs(z_score))
  return(c(z_score, p, rho0, rho1))
}
