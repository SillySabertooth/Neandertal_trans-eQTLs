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


setwd("~/Neandertal_trans-eQTLs/collapsing_SNPs_R2/")

library(data.table)
library(dplyr)
library(tidyr)


dt <- fread("all.txt.geno.ld") # a dataset where the R2 counted across aSNPs that are linked to TFs

# preprocess the dataset - leave SNPs above R2 threshold, mutate the table to be symmetrical 
colnames(dt)[5] <- "R2"

dt <- dt[dt$R2 >= 0.5,]
#dt <- dt[dt$R2 >= 0.8,]

dt$POS1 <- paste(dt$CHR,dt$POS1, sep = "_")
dt$POS2 <- paste(dt$CHR,dt$POS2, sep = "_")

dt <- rbind(dt, dt[,c(1,3,2,4,5)] %>% setnames(.,c("POS1","POS2"),c("POS2","POS1")))
dt %>% unique()



#chr=11 #debugging

tmp_list <- list() # the list for saving stuff

for (chr in unique(dt$CHR)){
  tmp <- dt[dt$CHR == chr,] # tmp table
  
  while (nrow(tmp) != 0){
    # we iterate through a given chr till all pairs of snps that have link >= threshold will be united in the one haplo, 
    # it works well for archaic haplos since they inherit in blocks, while for modern, the more full linkage could be quiried (all snps in the haplo have  r2 > 0.5),
    # on the edge states, it could unite positions when 1 <> 2 <> 3 but 1 <!> 3
    # and also the algorithm is not meant to distinguish very close haplotypes, it's more to collapse obvious related SNPs in an automatic way
  
    random_pos <- sample(tmp$POS1,1)
    tmp_list[random_pos] <- list(random_pos)
  
    while(nrow(tmp) != 0 & any(unlist(tmp_list[random_pos]) %in% tmp$POS1)){
      
      # attach to the list position that mathces position in the list
      tmp_list[random_pos] <- list(unique(c(unlist(tmp_list[random_pos]),
                                     tmp$POS2[tmp$POS1 %in% unlist(tmp_list[random_pos])])))
      
      # remove them from the table
      #tmp <- tmp[!tmp$POS1 %in% unlist(tmp_list[random_pos]),]
      tmp <- tmp[!tmp$POS2 %in% unlist(tmp_list[random_pos]),]
  
  }

  }
}


remove_list <- tmp_list # creates and fills a list of SNPs that can be removed from the original one in order to have a list of unrelated SNPs
#i <- "169099483" #debugging
for (i in names(remove_list)){
  remove_list[i] <- list(unlist(remove_list[i])[unlist(remove_list[i]) != i])
}

save_to <- as.character(unlist(remove_list)) # %>% length()
write.table(save_to,"remove_background_0.5.txt",col.names = F, row.names = F, quote = F, sep = "\t")

