# https://vbaliga.github.io/verify-that-r-packages-are-installed-and-loaded/
## load or install&load all
packages = c("data.table","dplyr")
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)


setwd("~/Neandertal_trans-eQTLs/aSNPs_proportion_enrichment/")
library(data.table)
library(dplyr)


dt_all <- fread("all.trans_eQTL.tested.tsv.gz")
removing <- fread("remove_10206.txt", header = F)$V1
dt_all_removed <- dt_all

dt_all_removed <- dt_all[!dt_all$chr_pos %in% removing,]
table(dt_all_removed$AssessedAllele_AF_rounded)

# 
dt_m <- fread("arch.trans_eQTL.tested.tsv")
dt_all_removed <- dt_all_removed[!dt_all_removed$chr_pos %in% dt_m$chr_pos,]

removing_small <- fread("remove_77_0.5.txt", header = F)$V1

dt_m_removed <- dt_m

dt_m_removed <- dt_m[!dt_m$chr_pos %in% removing_small,]


t_subset <- data.frame(table(dt_m_removed$AssessedAllele_AF_rounded))
dt <- fread("list_of_sign_SNPs.eQTLGen.tsv", header = F)$V1



# making the test
same_af <- dt_all_removed[dt_all_removed$AssessedAllele_AF_rounded %in% dt_m_removed$AssessedAllele_AF_rounded,]
table(same_af$AssessedAllele_AF_rounded)

#n=2
#might be 16
get_snps_bins <- function(n){
  temp <- same_af[same_af$AssessedAllele_AF_rounded == t_subset$Var1[n],16]$chr_pos
  return(sample(temp,t_subset$Freq[n]))
}

res <- c()
for (i in 1:100){
r1 <- unlist(lapply(1:nrow(t_subset), function(x) get_snps_bins(x)))
res <- c(res,sum(r1 %in% dt)) #simple associations
#res <- c(res,sum(r1 %in% height$chr_pos)) #simple height test

}


res
num <- sum(dt_m_removed$chr_pos %in% dt)
res_sum <- data.frame(
           median_5_95=paste0(median(res)," (",paste(quantile(res, probs = c(0.05, 0.95)), collapse = " - "),")"),
           d1=sum(res >= num),
           d2=sum(res <= num),
           pval=(min(c(sum(res >= num),sum(res <= num),
                       50))*2)/100,
           name="0.5_ilsntnfree_15aSNPs")



