#select genes from results of muscle skeletal tissue
select_pairs = function(n, p_cutoff, marcor_cutoff, tissue){
  Aresult = read.table(paste0('/gpfs/data/im-lab/nas40t2/tae/differentialeQTLs/codes/Leffect/marginal/result_global/',
                              tissue,
                              '_p.txt'), header=TRUE, stringsAsFactors=FALSE)
  genes = Aresult$genenames[which(Aresult$A_p < p_cutoff)]
  load(paste0('/gpfs/data/nicolae-lab/users/tae/expression/',
                          tissue,
                          '_exp.Rdata'))
  expr = obj
  rm(obj)
  ind = which(colSums(is.na(expr))==nrow(expr))
  expr = expr[, -ind]
  ind = which(rowSums(is.na(expr))==ncol(expr))
  expr = expr[-ind, ]
  exp2 = expr[genes, ]

  subjects = read.table('/gpfs/data/im-lab/nas40t2/tae/differentialeQTLs/data/summary/subjects_full436.txt',
                        stringsAsFactors = FALSE)$V1
  exp = exp2; rm(exp2); gc()

  sequ = list()
  step = round(nrow(exp)/5)
  for (i in 1:4){
    sequ[[i]] = (step * (i-1) + 1) : (step * i)
  }
  sequ[[5]] = (step*4+1) : nrow(exp)
  k = 1
  out = list()
  for (i in 1:5){
    for (j in 1:5){
      col1 = sequ[[i]]
      col2 = sequ[[j]]
      out[[k]] = matrix(0, length(col1), length(col2))
      out[[k]] = (exp[col1,]) %*% t(exp[col2, ]) / ncol(exp)
      k = k+1
      print(k)
    }
  }
  bigmat = matrix(0, nrow(exp), nrow(exp))
  for (i in 1:5){
    for (j in 1:5){
      bigmat[sequ[[i]], sequ[[j]]] = out[[5*(i-1)+j]]
    }
  }

  where = which(abs(bigmat)>marcor_cutoff, arr.ind = TRUE)
  where2 = where[where[,1] != where[,2], ]

  write.table(exp, paste0(tissue,'_selected_exp.txt'), col.names=TRUE, row.names=FALSE, quote = FALSE)
  write.table(where2, paste0(tissue,'_selected_pairs.txt'), col.names=TRUE, row.names=FALSE, quote = FALSE)
  # write.table(genes, paste0(tissue,'_selected_names.txt'), col.names=FALSE, row.names=FALSE, quote = FALSE)
}
