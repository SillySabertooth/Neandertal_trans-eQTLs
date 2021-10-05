# https://vbaliga.github.io/verify-that-r-packages-are-installed-and-loaded/
## load or install&load all
packages = c("data.table",
             "dplyr", "tidyr")
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)


# setwd("~/")

### process Human TFTG database

tfgt <- fread("tftg_tf.tsv", header = F, fill = T)
# some double-checked assignments according to ncbi
tfgt$V1[tfgt$V1 == "1-Oct"] <- "POU2F1"
tfgt$V1[tfgt$V1 == "2-Oct"] <- "POU2F2"
tfgt$V1[tfgt$V1 == "4-Oct"] <- "POU5F1"
tfgt$V1[tfgt$V1 == "2-Dec"] <- "BHLHE41"
unique(tfgt$V1)
colnames(tfgt) <- c("TF","gene_numb")

# mart_export.txt is a table from www.ensembl.org/biomart/martview/ with NCBI Gene ID and Gene name/NCBI gene (formerly Entrezgene) accession
mart <- fread("mart_export.txt")
colnames(mart) <- c("gene_numb","regulated_gene")
mart <- mart[!is.na(mart$gene_numb),]
mart$gene_numb %>% unique() %>% length()
mart$regulated_gene %>% unique() %>% length()

#
mart_tgft <- merge(tfgt, mart, all.x = T, allow.cartesian=TRUE)
mart_tgft <- mart_tgft %>% unique()
mart_tgft <- mart_tgft[!is.na(mart_tgft$TF)]
mart_tgft <- mart_tgft[!is.na(mart_tgft$gene_numb)]
mart_tgft$regulated_gene[is.na(mart_tgft$regulated_gene)] <- mart_tgft$gene_numb[is.na(mart_tgft$regulated_gene)] 

sum(is.na(mart_tgft$regulated_gene))

mart_tgft <- mart_tgft[,2:3]
mart_tgft$Human_TFTG_DB <- "Human_TFTG_DB"
#mart_tgft[duplicated(mart_tgft[,c(1,2)])]
write.table(mart_tgft,"Human_TFTG_DB.tsv", row.names = F, quote = F, sep = "\t")

