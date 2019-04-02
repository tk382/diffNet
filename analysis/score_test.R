library(R.utils)
sourceDirectory("R", modifiedOnly=FALSE)
Rcpp::sourceCpp("src/utilities.cpp")

# read data from the university HPC
library(data.table)
Y = setDF(fread("/Volumes/nicolae-lab-1/users/tae/Muscle_Skeletal_selected_exp.txt",
                header=T)) #selected expression of A effect pvalue < 0.1
pairs = read.table("/Volumes/nicolae-lab-1/users/tae/Muscle_Skeletal_selected_pairs.txt",
                   header=T) #among Y, marginal correlation > 0.7
A = read.table("/Volumes/im-lab/nas40t2/tae/differentialeQTLs/data/summary/finalA.txt", header=T,
               stringsAsFactors = FALSE)


colnames(Y) = gsub("[.]", "-", colnames(Y))
int = intersect(colnames(Y), A$subject)
Y = Y[, match(int, colnames(Y))]
Y = t(scale(t(Y)))
A = A[match(int, A$subject), ]

A = scale(A$A)

#### genes with ancestry effect AND high marginal correlation ####
s1 = rep(0,nrow(pairs))
for (i in 1:nrow(pairs)){
  y1 = as.numeric(Y[pairs[i,1], ])
  y2 = as.numeric(Y[pairs[i,2], ])
  y1 = residuals(lm(y1~A))
  y2 = residuals(lm(y2~A))
  s1[i] = get_score(y1, y2, A)
}
hist(s1, 100)
qqplot(rchisq(1000,1), s1); abline(0,1,col = 'red')



#### random pairs among genes with ancestry effect and without marcor ####
s2 = rep(0,1000)
for (i in 1:1000){
  inds = sample(nrow(Y), 2)
  y1 = as.numeric(Y[inds[1], ])
  y2 = as.numeric(Y[inds[2], ])
  y1 = resid(lm(y1~A))
  y2 = resid(lm(y2~A))
  s2[i] = get_score(y1, y2, A)
}
hist(s2,100)
qqplot(rchisq(1000,1), s2); abline(0,1,col='red')



#### entire expression ####
load("/Volumes/nicolae-lab/users/tae/expression/Muscle_Skeletal_exp.Rdata")
Y = obj
ind = which(colSums(is.na(Y))==nrow(Y))
Y = Y[,-ind]
ind = which(rowSums(is.na(Y))==ncol(Y))
Y = Y[-ind, ]
dim(Y)
colnames(Y) = gsub("[.]", "-", colnames(Y))
int = intersect(colnames(Y), A$subject)
Y = Y[, match(int, colnames(Y))]
Y = t(scale(t(Y)))
A = A[match(int, A$subject), ]
s3 = rep(0,1000)
for (i in 1:1000){
  inds = sample(nrow(Y), 2)
  y1 = as.numeric(Y[inds[1], ])
  y2 = as.numeric(Y[inds[2], ])
  y1 = resid(lm(y1~A$A))
  y2 = resid(lm(y2~A$A))
  s3[i] = get_score(y1, y2, A$A)
}
hist(s3,100)
qqplot(rchisq(1000,1), s3); abline(0,1,col='red')


#### perfect null simulated ####
library(MASS)
s = rep(0,10000)
for (i in 1:10000){
  y = mvrnormArma(56, c(0,0), matrix(c(1,.5,.5,1),2,2))
  y1 = y[,1]
  y2 = y[,2]
  x = rnorm(56)
  s[i]= get_score(x,y1,y2)
}
qqplot(rchisq(10000,1), s); abline(0,1,col='red')




#### genes with ancestry effect AND high marginal correlation ####
#### compute all pairs 734 * 734 and investigate row Sums ####
id = unique(c(pairs$row, pairs$col))
label = 1:length(id)
names(label) = id
s1 = rep(0,nrow(pairs))
Q = matrix(NA, length(id), length(id))
qq = qr(A)
YR = t(qr.resid(qq, t(Y)))

for (i1 in 1:(length(id)-1)){
  for (i2 in (i1):length(id)){
    ind1 = id[i1]
    ind2 = id[i2]
    y1 = as.numeric(YR[ind1,])
    y2 = as.numeric(YR[ind2,])
    Q[i1,i2] = get_score(y1, y2, A)
  }
}
qqplot(as.numeric(Q[upper.tri(Q)]))

Q[is.na(Q)] = 0
Q = Q + t(Q)

degree_stat = rowSums(Q)/(ncol(Q)-1)




