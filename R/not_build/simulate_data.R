simulate_data = function(mu,beta,n){
  #set.seed(1)
  #x = runif(n)
  x = seq(0,1,length=n)
  theta = exp(mu + beta * x)
  rho = (theta-1)/(theta+1)
  y = matrix(0,n,2)
  for (i in 1:n){
    y[i,] = mvrnorm(1, mu=c(0,0), matrix(c(1/2, rho[i]/2, rho[i]/2, 1/2), nrow=2))
  }
  ind = sort(x, index.return=TRUE)$ix
  x = x[ind]; rho = rho[ind]; y = y[ind, ]
  return(list(y=y,x=x,mu=mu,beta=beta,rho=rho))
}
