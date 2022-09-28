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

setwd("~/Neandertal_trans-eQTLs/TF_prediction_DBs")


load("func_and_data_for_change_gene_names.RData")
# the two tables are generated from this one https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz
# the function make a normalizing merge; you can assign some dt and column directly to variables inside and check what's going on
################################
# all commands here is better to launch one-by-one with thinking about what they return


# Human TFTG processing; different double-checks and merges with manifest and the reference table of synonyms

mart_tgft <- fread("Human_TFTG_DB.tsv")
manifest <- fread("humanTFMotifEntrezMappings.csv", select = 1:2) %>% unique()

tmp <- merge(syn_var, manifest, all.y = T)
sum(syn_var$EntrezID %in% manifest$EntrezID) #check gene below
data.frame(table(tmp$Gene)) %>% nrow() 
data.frame(table(tmp$TF_another_name)) %>% nrow() 

tmp[is.na(tmp$EntrezID),]
tmp[is.na(tmp$TF_another_name),]
tmp <- merge(syn_var, manifest)
sum(syn_var$EntrezID %in% tmp$EntrezID)
data.frame(table(tmp$Gene)) %>% nrow() 
data.frame(table(tmp$TF_another_name)) %>% nrow()


#check for artifacts

colnames(tmp)[3] <- "TF"
data.frame(table(tmp$TF)) %>% nrow() 
data.frame(table(tmp$TF_another_name)) %>% nrow()
mart_tgft$TF %>% unique() %>%  length()

tmp2 <- merge(mart_tgft, tmp, by = "TF")


tmp3 <- unique(tmp2[,-1])
data.frame(table(tmp3$TF)) %>% nrow() 
unique(tmp$TF_another_name)[!unique(tmp$TF_another_name) %in% unique(tmp3$TF)]

length(unique(tmp3$regulated_gene))
sum(unique(tmp3$regulated_gene) %in% syn$Symbol)
sum(unique(tmp3$regulated_gene) %in% syn$Synonyms)
nums <- unique(tmp3$regulated_gene)[!unique(tmp3$regulated_gene) %in% syn$Synonyms]
sum(nums %in% syn$GeneID)
#

tmp4 <- strange_merge(tmp3,"regulated_gene")


tmp3_5 <- tmp3[tmp3$regulated_gene %in% fls_from_func_last,]
unique(tmp3_5$regulated_gene) %>% length()
tmp3_5$regulated_gene <- as.numeric(tmp3_5$regulated_gene)
unique(tmp3_5$regulated_gene) %>% length()
tmp3_6 <- merge(tmp3_5, unique(syn[,c(1,2)]), by.x = "regulated_gene", by.y = "GeneID") %>% unique()
colnames(tmp3_6)[c(1,5)] <- c("GeneID_regulated_gene","regulated_gene")

tmp5 <- rbind(tmp4,tmp3_6) %>% unique()

sum(unique(tmp5$regulated_gene) %in% syn$Symbol)
sum(unique(tmp5$TF_another_name) %in% syn$Symbol)
sum(unique(tmp5$GeneID_regulated_gene) %in% syn$GeneID)
sum(unique(tmp5$EntrezID) %in% syn$GeneID)

z <- data.frame(table(unique(tmp5[,c(1,2)]) %>% select(regulated_gene))) %>% filter(Freq > 1)
z <- tmp5[tmp5$regulated_gene %in% z$Var1,c(1,2)] %>% unique()
tmp3[tmp3$regulated_gene %in% z$regulated_gene,1] %>% unique()

tmp5$GeneID_regulated_gene[tmp5$regulated_gene == "HBD"] <- 3045
tmp5$GeneID_regulated_gene[tmp5$regulated_gene == "MEMO1"] <- 7795
tmp5$GeneID_regulated_gene[tmp5$regulated_gene == "TEC"] <- 7006
tmp5 <- unique(tmp5)

tmp5 <- tmp5[,c(5,4,1,2,3)]
colnames(tmp5)[c(1,2)] <- c("TF","GeneID_TF")
length(unique(tmp5$GeneID_TF))
length(unique(tmp5$TF))
length(unique(tmp5$GeneID_regulated_gene))

write.table(tmp5,"Human_TFTG_DB_v1.tsv", row.names = F, quote = F, sep = "\t")
tgft <- fread("Human_TFTG_DB_v1.tsv")




###############




# Harmonizome processing; different double-checks and merges between subdatabases and the reference table of synonyms
harmony <- fread("Harmonizome_TFtargets.tsv")

harmony$all <- apply(harmony[,3:8], 1, function(x) toString(na.omit(x)))

harmony <- harmony[,c(1,2,9)] %>% unique()
harmony <- separate_rows(harmony, all, sep = ", ")
harmony <- unique(harmony)
table(harmony$all)
table(unique(harmony[,c(1,3)]) %>% select(all))

Motif <- harmony[harmony$all == "MotifMap_Pred",]
other <- harmony[!harmony$all == "MotifMap_Pred",]


# All Harmonizome DBs except MotifMap (since it has some issues)
other2 <- strange_merge(other,"TF")

length(unique(other2$regulated_gene))
sum(unique(other2$regulated_gene) %in% syn$Symbol)
sum(unique(other2$regulated_gene) %in% syn$Synonyms)

other3 <- strange_merge(other2,"regulated_gene")

sum(unique(other$TF) %in% syn$Symbol)
sum(unique(other2$TF) %in% syn$Symbol)
sum(unique(other3$TF) %in% syn$Symbol)
sum(unique(other3$TF) %in% syn$Synonyms)
sum(unique(other3$GeneID_TF) %in% syn$GeneID)

sum(unique(other2$regulated_gene) %in% syn$Symbol)
sum(unique(other3$regulated_gene) %in% syn$Symbol)
sum(unique(other3$regulated_gene) %in% syn$Synonyms)

sum(unique(other3$GeneID_regulated_gene) %in% syn$GeneID)

z <- data.frame(table(unique(other3[,c(1,2)]) %>% select(regulated_gene))) %>% filter(Freq > 1)
z <- other3[other3$regulated_gene %in% z$Var1,c(1,2)] %>% unique()
zz <- z[c(1,3,5,7,9),]

for (i in 1:nrow(zz)){
print(other3$GeneID_regulated_gene[other3$regulated_gene == zz$regulated_gene[i]])  
other3$GeneID_regulated_gene[other3$regulated_gene == zz$regulated_gene[i]] <- zz$GeneID_regulated_gene[i]
print(other3$GeneID_regulated_gene[other3$regulated_gene == zz$regulated_gene[i]])
}
other3 <- unique(other3)

write.table(other3,"Hamonize_Other_v1.tsv", row.names = F, quote = F, sep = "\t")
Harmon <- fread("Hamonize_Other_v1.tsv")



# MotifMap processing: workarounds with case, etc.

Motif <- harmony[harmony$all == "MotifMap_Pred",]

Motif$TF[!(Motif$TF %in% syn$Synonyms) & (toupper(Motif$TF) %in% syn$Synonyms)] <- 
  toupper(Motif$TF)[!(Motif$TF %in% syn$Synonyms) & (toupper(Motif$TF) %in% syn$Synonyms)]

Motif$TF[!(Motif$TF %in% syn$Synonyms) & (tolower(Motif$TF) %in% syn$Synonyms)] <- 
  tolower(Motif$TF)[!(Motif$TF %in% syn$Synonyms) & (tolower(Motif$TF) %in% syn$Synonyms)]

Motif <- unique(Motif)

other2 <- strange_merge(Motif,"TF")

length(unique(other2$regulated_gene))
sum(unique(other2$regulated_gene) %in% syn$Symbol)
sum(unique(other2$regulated_gene) %in% syn$Synonyms)

other3 <- strange_merge(other2,"regulated_gene")

fls_from_func_last

sum(unique(Motif$TF) %in% syn$Symbol)
sum(unique(other2$TF) %in% syn$Symbol)
sum(unique(other3$TF) %in% syn$Symbol)
sum(unique(other3$TF) %in% syn$Synonyms)
sum(unique(other3$GeneID_TF) %in% syn$GeneID)

sum(unique(other2$regulated_gene) %in% syn$Symbol)
sum(unique(other3$regulated_gene) %in% syn$Symbol)
sum(unique(other3$regulated_gene) %in% syn$Synonyms)

sum(unique(other3$GeneID_regulated_gene) %in% syn$GeneID)

z <- data.frame(table(unique(other3[,c(1,2)]) %>% select(regulated_gene))) %>% filter(Freq > 1)
z <- other3[other3$regulated_gene %in% z$Var1,c(1,2)] %>% unique()
zz <- z[c(1,3,5,7),]

for (i in 1:nrow(zz)){
  print(other3$GeneID_regulated_gene[other3$regulated_gene == zz$regulated_gene[i]])  
  other3$GeneID_regulated_gene[other3$regulated_gene == zz$regulated_gene[i]] <- zz$GeneID_regulated_gene[i]
  print(other3$GeneID_regulated_gene[other3$regulated_gene == zz$regulated_gene[i]])
}
other3 <- unique(other3)

write.table(other3,"Motiff_v1.tsv", row.names = F, quote = F, sep = "\t")
M <- fread("Motiff_v1.tsv")


# AnimalTFDB gene list check
#http://bioinfo.life.hust.edu.cn/AnimalTFDB/#!/download

anim <- fread("Homo_sapiens_TF.txt")
anim2 <- strange_merge(anim,"Symbol")

anim2$Symbol[!anim2$Symbol %in% anim$Symbol]

write.table(anim2,"Homo_sapiens_TF_adj.txt", quote = F, row.names = F, sep = "\t")



# merge all together:

## to TFs' table
colnames(tgft)[5] <- "all"
the_whole <- rbind(M,Harmon,tgft)
write.table(the_whole,"TF_database_regGenes_v1.tsv", row.names = F, quote = F, sep = "\t") # it's good to gzip it

## to a TFs list

animal <- fread("Homo_sapiens_TF_adj.txt", select = c(1,2)) %>% unique()
animal$all <- "Animal_TFDB"
colnames(animal)[1:2] <- c("TF","GeneID_TF")
to_r <- the_whole[,c(3,4,5)] %>% unique()
res <- rbind(to_r,animal) %>% unique()


write.table(res,"all_databses_TF_list_v1.txt", quote = F, row.names = F, sep = "\t")





### general overview; plotting intersection of TFs between different databases

ALL <- fread("all_databses_TF_list_v1.txt")

pairs <- data.frame(t(combn(unique(ALL$all), 2)))
pairs <- rbind(pairs,tibble(X1=unique(ALL$all),X2=unique(ALL$all)))

get_int <- function(i){
  x <- pairs[i,]
  return(length(intersect(ALL$TF[ALL$all == x$X1[1]],ALL$TF[ALL$all == x$X2[1]])))
}

inters <- lapply(1:nrow(pairs), function(x) get_int(x)) %>%
  do.call(rbind, .)

pairs_res <- cbind(pairs,inters)


p_tmp <- pairs_res[,c(2,1,3)]
colnames(p_tmp)[1:2] <- c("X1","X2")


all_tfs <- data.frame(table(ALL$all)) 
pairs_res_plot <- rbind(pairs_res, p_tmp) %>% unique()

pairs_res_plot <- merge(pairs_res_plot,all_tfs, by.x = "X1", by.y = "Var1")
pairs_res_plot$proc <- round((pairs_res_plot$inters / pairs_res_plot$Freq)*100,2)
pairs_res_plot$proc <- gsub("^","(",pairs_res_plot$proc)
pairs_res_plot$proc <- gsub("$",")",pairs_res_plot$proc)
pairs_res_plot$proc[pairs_res_plot$X1 == pairs_res_plot$X2] <- ""
pairs_res_plot$inters_N <- paste(pairs_res_plot$inters, pairs_res_plot$proc, sep = " ")

pairs_res_plot$X1 <- factor(pairs_res_plot$X1, levels = c("CHEA","ENCODE","TRANSFAC_Cur","TRANSFAC_Pred",
                                                          "JASPAR_Pred","MotifMap_Pred","Human_TFTG_DB","Animal_TFDB"))
pairs_res_plot$X2 <- factor(pairs_res_plot$X2, levels = c("CHEA","ENCODE","TRANSFAC_Cur","TRANSFAC_Pred",
                                                          "JASPAR_Pred","MotifMap_Pred","Human_TFTG_DB","Animal_TFDB"))

dt <- pairs_res_plot %>% 
  select(X1,X2,inters_N) %>% 
  spread(key=X2,sep="",value=inters_N)
colnames(dt) <- gsub("X2","",colnames(dt))

write.table(dt,"heatmap_table.tsv", sep = "\t", row.names = F, quote = F)

ggplot(pairs_res_plot, aes(X1,X2)) +
  geom_tile(aes(fill = inters)) + 
  geom_text(aes(label = inters)) +
  scale_fill_gradient(low = "orange", high = "midnightblue")+
  theme(axis.title.x =element_blank(),axis.title.y =element_blank())
