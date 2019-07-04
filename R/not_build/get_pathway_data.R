library(KEGGREST)
browseVignettes("KEGGREST")
out = keggGet("H00079")
out = out[[1]]
genes = sapply(strsplit(out$GENE, " [(]"), function(x) x[1])

library(biomaRt)
ensembl = useMart("ENSEMBL_MART_ENSEMBL",
                  dataset = "hsapiens_gene_ensembl")
results = getBM(attributes = c('hgnc_symbol','ensembl_gene_id'),
                filters = "hgnc_symbol",
                values = genes,
                mart = ensembl)


results$ensembl_gene_id


