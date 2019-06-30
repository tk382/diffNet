
solve_tmp = function(x){
  xx = matrix(0,2,2);
  xx[1,1] = length(x);
  xx[1,2] = sum(x);
  xx[2,1] = sum(x);
  xx[2,2] = sum(x^2);
  xx = solve(xx)
  return(xx)
}

get_score = function(x, y1, y2){
  n = length(x)
  if(length(x)!=length(y1) | length(x)!= length(y2) | length(y1)!=length(y2)){
    stop("x, y1, y2 should have the same length")
  }
  w = y1+y2
  xx = solve_tmp(x)
  xt = cbind(rep(1,n), x)
  hatsigma = sum(w^2)/n
  lhs = matrix(c(sum(w^2/hatsigma-1), sum((w^2/hatsigma-1)*x)), ncol=2)
  S = lhs %*% xx %*% t(lhs)
  S = S/2
  return(S)
}


