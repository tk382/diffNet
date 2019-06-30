continuous_mubeta = function(x, w){
  n = length(x)
  x_start = c(0,0)
  return(nleqslv(x_start, both_foc, control=list(ftol=1e-10)))
}

continuous_mu = function(x,w){
  return(log(sum(w^2)/sum(2-w^2)))
}

first_foc = function(mu, beta, x, w){
  theta = exp(mu+beta*x)
  w2    = w^2
  return(sum(2/(theta + 1) - w2/theta))
}

second_foc = function(mu, beta, x, w){
  theta = exp(mu+beta*x)
  w2    = w^2
  return(sum(2*x/(theta+1) - w2*x/theta))
}

both_foc = function(mubeta){
  mu = mubeta[1]
  beta = mubeta[2]
  return(c(first_foc(mu,beta,x,w), second_foc(mu,beta,x,w)))
}

visualize_first_foc = function(murange, betarange, x, w){
  n = length(murange); m = length(betarange)
  mat = matrix(0, n*m, 3)
  mat[,1] = rep(murange, each=m)
  mat[,2] = rep(betarange, n)
  for (i in 1:nrow(mat)){
    mat[i,3] = first_foc(mu=mat[i,1], beta=mat[i,2], x, w)
  }
  #scatterplot3d(x=mat[,2], y=mat[,1], mat[,3], ylab='mu', xlab='beta')
  return(mat)
}

visualize_second_foc = function(murange, betarange, x, w){
  n = length(murange); m = length(betarange)
  mat = matrix(0, n*m, 3)
  mat[,1] = rep(murange, each=m)
  mat[,2] = rep(betarange, n)
  for (i in 1:nrow(mat)){
    mat[i,3] = second_foc(mu=mat[i,1], beta=mat[i,2], x, w)
  }
  scatterplot3d(x=mat[,2], y=mat[,1], z=mat[,3], ylab='mu', xlab='beta')
  return(mat)
}


