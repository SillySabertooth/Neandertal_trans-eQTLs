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

# the idea is to test whether we see significant trans-eQTLs among variants of modern human ancestry 
# as often as in archaic ones taking into account frequency-matched


##### pre processing
# for the first, take all analyzed variants from eQTLGen (we did the collapsing beyond of this analysis)
# https://www.eqtlgen.org/trans-eqtls.html
dt_all <- fread("all.trans_eQTL.tested.tsv.gz")
removing <- fread("remove_10206.txt", header = F)$V1 # read the preformed vector by which we can leave only one SNP per linked by R2 SNPs; made with collapsing_SNPs_R2
dt_all_removed <- dt_all

dt_all_removed <- dt_all[!dt_all$chr_pos %in% removing,]
table(dt_all_removed$AssessedAllele_AF_rounded) # check the number of SNPs per bin of AF

# do the same with archaic variants
dt_m <- fread("arch.trans_eQTL.tested.tsv")
dt_all_removed <- dt_all_removed[!dt_all_removed$chr_pos %in% dt_m$chr_pos,] # remove archaic variants from the vector of all SNPs

removing_small <- fread("remove_77_0.5.txt", header = F)$V1
dt_m_removed <- dt_m
dt_m_removed <- dt_m[!dt_m$chr_pos %in% removing_small,]


t_subset <- data.frame(table(dt_m_removed$AssessedAllele_AF_rounded)) # create a table of archaic SNPs AF bins 
dt <- fread("list_of_sign_SNPs.eQTLGen.tsv", header = F)$V1 # read the list of significant trans-eQTL SNPs



##### the test
same_af <- dt_all_removed[dt_all_removed$AssessedAllele_AF_rounded %in% dt_m_removed$AssessedAllele_AF_rounded,] # leave in the table only AF matched variants
table(same_af$AssessedAllele_AF_rounded) #check the result

#n=2 #debugging
get_snps_bins <- function(n){
  # here we take all modern SNPs that fall into the given AF bin, and do a subset to match the number of archaic SNPs in the given AF bin
  temp <- same_af[same_af$AssessedAllele_AF_rounded == t_subset$Var1[n],16]$chr_pos # 16 is the last column with chr_pos
  return(sample(temp,t_subset$Freq[n]))
}

res <- c() # a empty vector to safe the result, we will test 100 random AF and size matches sets of modern human variants to be significant trans-eQTLs
for (i in 1:100){
r1 <- unlist(lapply(1:nrow(t_subset), function(x) get_snps_bins(x))) # iterate through all bins to create the similar AF distribution as archaic SNPs have
res <- c(res,sum(r1 %in% dt)) #simple associations # check how many of them are significant
#res <- c(res,sum(r1 %in% height$chr_pos)) #simple height test

} # take ~10 seconds


res # ! since we don't put here any seed, due to random, each launch the result of the test will slightly differ
num <- sum(dt_m_removed$chr_pos %in% dt) # how many of archaic variants are significant trans-eQTLs

# res_sum includes number of significant archaic trans-eQTLs, 
# median and 5/95CI of significant trans-eQTLs of modern origin in 100 random sets
# number of sets that have more or equal significant variants and less or equal
# the pval for two sided test
res_sum <- data.frame(
           aSNP_sign_numb=num,
           median_5_95_CI=paste0(median(res)," (",paste(quantile(res, probs = c(0.05, 0.95)), collapse = " - "),")"),
           d1=sum(res >= num),
           d2=sum(res <= num),
           pval=(min(c(sum(res >= num),sum(res <= num),
                       50))*2)/100,
           name="0.5R2_collapsed")
