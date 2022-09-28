# https://vbaliga.github.io/verify-that-r-packages-are-installed-and-loaded/
## load or install&load all
packages = c("data.table", "ggplot2",
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

# setwd("~/") # 
# args <- c("chr12_56740682.for_r2.tmp.geno.ld","12","56740682") # for the debugging

args <- commandArgs(trailingOnly = TRUE)
data <- fread(args[1])
colnames(data)[5] <- "R2"


# A. plot the whole initial region
# https://stackoverflow.com/questions/15624656/label-points-in-geom-point

# make false string to plot target aSNP
fake_str <- tibble(CHR=args[2],
                   POS1=as.numeric(args[3]),
                   POS2=as.numeric(args[3]),
                   N_INDV=2504,
                   R2=5)
plotty <- rbind(data,fake_str)


# add the fiest heatmap to the plotting list  
plot_list = list()
i <- 1
plot_list[[i]] = ggplot(data = plotty, aes(x=as.character(POS1), y=as.character(POS2), fill=R2)) + 
  geom_tile()+
  ggtitle(paste0("provided_region.",gsub(".tmp.geno.ld","",args[1])))+
  xlab("POS1")+
  ylab("POS2")+
  scale_fill_gradient('R2', limits=c(0, 1), breaks = c(0, 0.2, 0.5, 0.8, 1),  low = "yellow", high = "darkblue", na.value = "red")
i=i+1




# B. defining a haplotype

# prehotfix N1: seprate table and sort it for both direction from the target aSNP
snps_up <- data[data$POS1 == args[3],]
snps_up <- snps_up[order(snps_up$POS2),]

snps_down <- data[data$POS2 == args[3],]
snps_down <- snps_down[order(-snps_down$POS1),]

# iteration function
iterate <- function(vect, tresh){
  var = 1
  num = 1
  while (var == 1) {
    val = vect[num+1]
    if (length(vect) == num){
      num=num+1
      var = 0
    } else if (vect[num] >= tresh | val >= tresh){
      num=num+1
      next
    } else {var = 0}
  }
  return(num-1)
}

# tr=0.5 # for the debugging
# here we infer haplotype for two R2 tresholds: 0,5 and 0,8

for (tr in c(0.5,0.8)){
  # check whether there is SNPs up, and if yes, write down SNP position
  if (iterate(snps_up$R2, tr) == 0){
    pos_up <- as.numeric(args[3])
  } else {
    pos_up <- snps_up$POS2[iterate(snps_up$R2, tr)]
  }
  # getting some statistic about the upper part of the inferred haplotype
  # SNPs that have R2 lower than treshold in of the inferred haplotype
  pos_up_bad <- snps_up[snps_up$POS2 <= pos_up & snps_up$R2 < tr,c(1,3,5)] 
  colnames(pos_up_bad)[2] <- "POS"
  # check how much from 5 SNPs out of the haplotype bounds in a row have a lower threshold than needed 
  # should be more than 1, bigger is better
  pos_up_tail <- sum(snps_up$R2[(iterate(snps_up$R2, tr)+1):(iterate(snps_up$R2, tr)+5)] < tr) 
  
  
  # the same idea for the downstream part of the region
  if (iterate(snps_down$R2, tr) == 0){
    pos_down <- as.numeric(args[3])
  } else {
    pos_down <- snps_down$POS1[iterate(snps_down$R2, tr)]
  }

  pos_down_bad <- snps_down[snps_down$POS1 >= pos_down & snps_down$R2 < tr,c(1,2,5)]
  colnames(pos_down_bad)[2] <- "POS"
  pos_down_tail <- sum(snps_down$R2[(iterate(snps_down$R2, tr)+1):(iterate(snps_down$R2, tr)+5)]< tr)
  
  
  # unite the positions with low R2
  pos_bad <- rbind(pos_up_bad,pos_down_bad)
  
  # create a table with SNPs fall into the inferred haplotype
  on <- snps_up[snps_up$POS2 <= pos_up,c(1,3,5)] %>% setnames(.,"POS2","POS")
  tw <- snps_down[snps_down$POS1 >= pos_down,c(1,2,5)] %>% setnames(.,"POS1","POS")
  whole_t <- rbind(tw[order(tw$POS),],on)
  
  # wrap up all statistics per the inferred haplotype
  statis <- tibble(snp = paste0(args[2],":",args[3]),
                   treshold_r2=tr,
                   chr=unique(data$CHR),
                   str=pos_down,
                   end=pos_up,
                   size=pos_up - pos_down,
                   aSnps_number=nrow(whole_t),
                   num_of_low_r2_in_region = nrow(pos_bad),
                   HowMuchFivePosBefore = pos_down_tail,
                   HowMuchFivePosAfter = pos_up_tail)
  # the last statistic looks at whether the positions that lay right around the chosen haplotype have r2 with our given aSNP less than the threshold
  # if yes, then it will be counted, so 5 is good and means that all was good
  
  name <- paste0(args[2],"_",args[3],".",pos_down,"_",pos_up,".R2_",tr) # create big name
  
  write.table(whole_t,paste0(name,".snps_r2.tsv"), row.names = F, quote = F, sep = "\t") # write the haplotype SNPs (with low R2!) 
  write.table(pos_bad,paste0(name,".low_r2_aSNP_in_the_region.tsv"), row.names = F, quote = F, sep = "\t") # write only low R2 SNPs in the region
  write.table(statis,paste0(name,".stat.tsv"), row.names = F, quote = F, sep = "\t") # write statistics
  
  
  # create a table for final plotting of inferred haplotype
  ready <- data[data$POS1 %in% pos_down:pos_up & data$POS2 %in% pos_down:pos_up,]
  ready <- rbind(ready,fake_str)
  ready$POS1 <- as.character(ready$POS1)
  ready$POS2 <- as.character(ready$POS2)
  
  
  plot_list[[i]] = ggplot(data = ready, aes(x=POS1, y=POS2, fill=R2)) + 
    geom_tile()+
    ggtitle(name)+
    scale_fill_gradient('R2', limits=c(0, 1), breaks = c(0, 0.2, 0.5, 0.8, 1),  low = "yellow", high = "darkblue", na.value = "red")
  i=i+1
}

# https://stackoverflow.com/questions/26034177/save-multiple-ggplots-using-a-for-loop
# save plots as PDF
pdf(paste0(args[2],"_",args[3],".full_longest_haplo_r2_0.5_0.8.pdf"))
for (i in c(1,2,3)) {
  print(plot_list[[i]])
}
dev.off()