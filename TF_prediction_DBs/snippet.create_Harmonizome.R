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


setwd("~/Neandertal_trans-eQTLs/TF_prediction_DBs")

# merge Harmonizome databases, write raw result into a file
# lasts long, better to laucnh as a script (it means in the terminal not in Rstudio)

dt_list <- fread("df_list", header = F)$V1

create_a_row <- function(i,data,name){
  row <- as.character(data[i,3:length(data)])
  row <- row[row != ""]
  data.table(data[i,1],row, name)
}

df <- fread(paste0(dt_list[1],".gmt.gz"), fill = TRUE)
dt <- lapply(1:nrow(df), function(x) create_a_row(x,df,dt_list[1])) %>%
  do.call(rbind, .)
colnames(dt) <- c("TF","regulated_gene",dt_list[1])
total <- dt

for (i in c(2,3,4,5,6)){
df <- fread(paste0(dt_list[i],".gmt.gz"), fill = TRUE)
dt <- lapply(1:nrow(df), function(x) create_a_row(x,df,dt_list[i])) %>%
  do.call(rbind, .)
colnames(dt) <- c("TF","regulated_gene",dt_list[i])
total <- merge(dt,total, all = T, by = c("TF","regulated_gene"))
}


write.table(total,"Harmonizome_TFtargets.tsv", row.names = F, quote = F, sep = "\t")
total$TF %>% unique() %>% length()