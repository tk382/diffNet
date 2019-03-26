#### honda correction ####
get_A1 = function(k,p,TT,RR, n){
  first = 24 * (k-1) * (p-1)
  second = -24 * n * sum(diag(TT) * diag(RR))
  third = 0
  for (i in 1:n){
    for (j in 1:n){
      third = third + TT[i,j] * (RR[i,i] * RR[j,j] + 2 * RR[i,j]^2)
    }
  }
  third = third * 6 * n
  return(first + second + third)
}


get_A2 = function(k, p, TT, RR, n){
  first = -24 * (p-1)^2
  second = 36 * n * sum(diag(TT)^2)
  third = -48 * sum(TT^2)
  fourth = 0
  for (i in 1:n){
    for (j in 1:n){
      fourth = fourth + TT[i,j] * TT[j,j] * RR[i,i]
    }
  }
  fourth = fourth * (-24)* n
  return(first + second + third + fourth)
}

get_A3 = function(k,p,TT,RR,n){
  first = 0
  for (i in 1:n){
    for (j in 1:n){
      first = first + TT[i,i] * TT[j,j] * TT[i,j]
    }
  }
  first = first * 24 * n
  second = 16 * n * sum(TT^3)
  return(first + second)
}

final_correction = function(A, C){
  #A = as.numeric(scale(rnorm(100)) * sqrt(100) / sqrt(99))
  n = length(A)
  X = cbind(rep(1,n), A)
  barA = mean(A) #always 0
  H = sum((A-mean(A))^2) #always n
  TT = matrix(NA, n, n)
  diag(TT) = (A-mean(A))^2 / H
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      TT[i,j] = TT[j,i] = (A[i] - mean(A))*(A[j] - mean(A)) / H
    }
  }
  RR = X %*% solve(crossprod(X)) %*% t(X)

  k = p = ncol(X)
  A1 = get_A1(k,p,TT,RR,n)
  A2 = get_A2(k,p,TT,RR,n)
  A3 = get_A3(k,p,TT,RR,n)

  first = A3 * C / (p-1) / (p+1) / (p+3) * (C^2 + (p+3)*C + (p+1)*(p+3))
  second = C * (C+p+1) * (A2-3*A3) / (p-1) / (p+1)
  third = C * (3 * A3 - 2 * A2 + A1) / (p-1)
  return((first + second + third) / 12 / n)
}

