library(biomaRt)

args = commandArgs(trailingOnly=T)

chr = args[1] #"19"
start = args[2] #"q13.32"
end = args[3] #"q13.32"
out.file = args[4]

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
results <- getBM(attributes = c("external_gene_id"),
                 filters = c("chromosome_name", "band_start", "band_end"), values = list(chromosome_name=chr, band_start=start, band_end=end),
                 mart = mart)
write.table(results, out.file, quote=F, row.names=F, col.names=F)

