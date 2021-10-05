# https://vbaliga.github.io/verify-that-r-packages-are-installed-and-loaded/
## load or install&load all
packages = c("data.table", "ggplot2",
             "dplyr", "tidyr","gridExtra")
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
# args <- c("12.56627300_56753822.final.vcf.gz",56740682, "0|0") # for the debugging

# read data
args <- commandArgs(trailingOnly = TRUE)
name <- args[1]
df <- fread(name, skip = "#CHR")
name <- gsub(".vcf.gz",".seq_comp", name)

# positions for removing
bad_pos <- rbind(df[df$AltaiNeandertal == "./." |
                df$Vindija33.19 == "./." |
                df$Denisova == "./.",c(1,2)],
                
                df[df$AltaiNeandertal == "0/1" |
                     df$Vindija33.19 == "0/1" |
                     df$Denisova == "0/1",c(1,2)]) %>% unique()

# check the difference between two Neanderthals in the region for your own
comp_Neand <- tibble(pos = gsub(".1000G.seq_comp","", name),
                     total = nrow(df),
                     Neand = sum(df$AltaiNeandertal == "1/1"),
                     Vindi = sum(df$Vindija33.19 == "1/1"),
                     diff = sum(df$AltaiNeandertal != df$Vindija33.19))




# two extremely similar function to separate samples into haplotypes according to YRI geno in the target position
# the first gather Neand haplotypes, second - Modern ones

# i=13 # for the debugging
# dt=df # for the debugging
split_haplo <- function(i,dt){
  geno <- dt[[i]][(dt$POS == args[2])]
  # fixme/workaround for case when REF is ancient haplo
  if (args[3] == "1|1"){
     geno <- gsub("2","1",
          gsub("1","0",
               gsub("0","2",geno))) 
  }
  
  x <- data.frame(dt[[i]])
  smpl <- colnames(dt)[i]
  colnames(x) <- smpl
  x[[smpl]] <- as.character(x[[smpl]])
  x[[smpl]] <- gsub("./.","0|0",x[[smpl]])
  
  x <- separate(x,
                colnames(x),
                into = c(paste0(smpl,"_hap1"),
                         paste0(smpl,"_hap2")),
                sep = "\\|")
  x <- data.frame(lapply(x, as.numeric))
  if (geno == "0|1"){
    proper_haplo <- c(F,T)
  } else if (geno == "1|0"){
    proper_haplo <- c(T,F)
  } else if (geno == "1|1"){
    proper_haplo <- c(T,T)
  } else {
    proper_haplo <- c(F,F)
  }
  return(x[,proper_haplo, drop = FALSE])
}

split_haplo_ref <- function(i,dt){
  geno <- dt[[i]][(dt$POS == args[2])]
  #fixme for case when REF is ancient haplo | sorry for this
  if (args[3] == "1|1"){
    geno <- gsub("2","1",
                 gsub("1","0",
                      gsub("0","2",geno))) 
  }
  x <- data.frame(dt[[i]])
  smpl <- colnames(dt)[i]
  colnames(x) <- smpl
  x[[smpl]] <- as.character(x[[smpl]])
  x[[smpl]] <- gsub("./.","0|0",x[[smpl]])
  x <- separate(x,
                colnames(x),
                into = c(paste0(smpl,"_hap1"),
                         paste0(smpl,"_hap2")),
                sep = "\\|")
  x <- data.frame(lapply(x, as.numeric))
  if (geno == "0|1"){
    proper_haplo <- c(T,F)
  } else if (geno == "1|0"){
    proper_haplo <- c(F,T)
  } else if (geno == "0|0"){
    proper_haplo <- c(T,T)
  } else {
    proper_haplo <- c(F,F)
  }
  return(x[,proper_haplo, drop = FALSE])
}

# !!!! 13 column is the first column with 1KG sample
# splitting samples into haplotypes
new_df <- lapply(13:length(df), function(x) split_haplo(x, df)) %>%
  do.call(cbind, .)
new_df_ref <- lapply(13:length(df), function(x) split_haplo_ref(x, df)) %>%
  do.call(cbind, .) 

# I guess the reason why I remove them here is that sometimes target snps can fall into the bad_pos category? #workaround
new_df <- new_df[!(df$`#CHROM` %in% bad_pos$`#CHROM` & df$POS %in% bad_pos$POS),]
new_df_ref <- new_df_ref[!(df$`#CHROM` %in% bad_pos$`#CHROM` & df$POS %in% bad_pos$POS),]
df <- df[!(df$`#CHROM` %in% bad_pos$`#CHROM` & df$POS %in% bad_pos$POS),]

# make empty table for results
stats <- tibble(medn="",all_pos="",rel="",who="")[-1,]
distr <- tibble(diff="",who="",type="",samples="")[-1,]



# count mismatches between modern and archaic humans 

#col="AltaiNeandertal" # for the debugging
for (col in c("AltaiNeandertal","Vindija33.19","Denisova")){

new_df_ins <- cbind(df %>% select(col),new_df)

new_df_ins[[col]][new_df_ins[[col]] == "1/1"] <- 1
new_df_ins[[col]][new_df_ins[[col]] == "0/0"] <- 0


tmp <- data.table(lapply(2:length(new_df_ins), function(x) sum(new_df_ins[[col]] != new_df_ins[[x]])) %>% unlist())
colnames(tmp) <- "diff" #col
tmp$who <- col
tmp$type <- "aHaplo"
tmp$samples <- colnames(new_df)
distr <- rbind(distr, tmp)

stats <- rbind(stats,
               tibble(medn = median(as.numeric(distr$diff[distr$who == col & distr$type == "aHaplo"])),
                      all_pos = nrow(new_df_ins), 
                      rel = median(as.numeric(distr$diff[distr$who == col & distr$type == "aHaplo"]))/nrow(new_df_ins), 
                                                           who = col))

new_df_ref_ins <- cbind(df %>% select(col),new_df_ref)
new_df_ref_ins[[col]][new_df_ref_ins[[col]] == "1/1"] <- 1
new_df_ref_ins[[col]][new_df_ref_ins[[col]] == "0/0"] <- 0

tmp <- data.table(lapply(2:length(new_df_ref_ins), function(x) sum(new_df_ref_ins[[col]] != new_df_ref_ins[[x]])) %>% unlist())
colnames(tmp) <- "diff" #col
tmp$who <- col
tmp$type <- "ref"
tmp$samples <- colnames(new_df_ref)
distr <- rbind(distr, tmp)
}




# adjust, plot and save results
distr$diff <- as.numeric(distr$diff)
distr$who <- factor(distr$who, levels = c("AltaiNeandertal","Vindija33.19","Denisova"))

write.table(distr,paste0(name,".distr.tsv"), quote = F, row.names = F, sep = "\t")
write.table(stats,paste0(name,".stats.tsv"), quote = F, row.names = F, sep = "\t")

###hists plots
jpeg(paste0(name,".hist.jpeg"),width = 1200, height = 839, res = 110)
ggplot(distr, aes(diff, fill = type))+
  geom_histogram(binwidth = 0.5)+
  facet_grid(rows = vars(who))
dev.off()
