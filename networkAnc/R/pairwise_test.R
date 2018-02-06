pairwise_test = function(L1,L2,y1,y2,n1=NA,n2=NA){
  if(is.na(n1)){n1 = length(L1)}
  if(is.na(n2)){n2 = length(L2)}
  if(length(L1)!=length(y1) || length(L2)!=length(y2)){
    stop('local ancestry and expression level should have equal sample size')
  }


}
