library(dplyr)
gencode33 <- read.table("gencode.v33.annotation.gtf", sep="\t", header=FALSE,  check.names = FALSE, skip = 5)

gencode33 <- gencode33 %>%
  filter(gencode33$V3=="gene")

  #write.table(gencode33, "gencode.v33.txt") ## make some modification in Notepad++
#gencode33 <- read.table("gencode.v33.txt", sep=" ", header=FALSE,  check.names = FALSE)
gene_id <- sapply(strsplit(gencode33$V9,split = " "),'[',2)
gene_id <- gsub(";", "", gene_id)
gene_type <- sapply(strsplit(gencode33$V9,split = " "),'[',4)
gene_type <- gsub(";", "",gene_type)
gene_name <- sapply(strsplit(gencode33$V9,split = " "),'[',6)
gene_name <- gsub(";", "", gene_name)
hgnc_id <- sapply(strsplit(gencode33$V9,split = " "),'[',10)
hgnc_id <- gsub(";", "", hgnc_id)
hgnc_id <- gsub("HGNC:", "", hgnc_id)

gencode33 <- cbind(gencode33, gene_id, gene_name, gene_type)
gencode33$V9 <- NULL
gencode33$V6 <- NULL; gencode33$V8 <-NULL; gencode33$V3 <- NULL
colnames(gencode33) <- c("chr", "annotation", "start", "end", "strand", "gene_id","gene_name", "gene_type")
rm(gene_id, gene_name, gene_type, hgnc_id, symbol)

save.image("gencode.v33.Rdata")

#####################################################
#####################################################
gtf <- rtracklayer::import('gencode.v33.annotation.gtf')
#gtf_df=as.data.frame(gencode.v33.annotation.gtf)
gtf <- as.data.frame(gtf)