####################################
####################################
## Enrichment analyses of predicted targets in different genes sets

####################################
## table with TF target prediction database information
## see scripts in TF_predict_dbs for information how this table was generated
tf1=read.table("Human_TFTG_DB_and_Harmonizome.tsv.gz",head=T,as.is=T,sep="\t")

## > head(tf1,n=4)
##      TF regulated_gene TRANSFAC_Pred TRANSFAC_Cur MotifMap_Pred JASPAR_Pred
## 1 ACAAT            A2M          <NA>         <NA>          <NA>        <NA>
## 2 ACAAT        ABHD17B          <NA>         <NA>          <NA>        <NA>
## 3 ACAAT         ACOT12          <NA>         <NA>          <NA>        <NA>
## 4 ACAAT          ACOT7          <NA>         <NA>          <NA>        <NA>
##   ENCODE CHEA Human_TFTG_DB
## 1   <NA> <NA> Human_TFTG_DB
## 2   <NA> <NA> Human_TFTG_DB
## 3   <NA> <NA> Human_TFTG_DB
## 4   <NA> <NA> Human_TFTG_DB

## Enrichment of a-cTF predicted target genes in archaic deserts
## relates to Fig. 2 panels a and b and result section 'The regulatory reach of predicted trans-eQTLs'

## definition of desert genes
## gene_coordinates_hg19_sort.txt is obtained from https://www.ensembl.org/biomart/martview
## cat  deserts.bed
## 1	104000000	114900000
## 3	76500000	91400000
## 7	108000000	128000000
## 8	54500000	65400000
## 13	49000000	61000000

system("bedtools intersect -a deserts.bed -b gene_coordinates_hg19_sort.txt -wo -F 1 > desert_genes_overlap.txt")

dg=read.table("desert_genes_overlap.txt",as.is=T,sep="\t")

## set of protein coding genes located in archaic deserts
dgg=unique(dg[,7])

## fisher test ORs
or.m=p.m=matrix(NA,ncol=7,nrow=length(unique(tf1[,1])))
rownames(or.m)=rownames(p.m)=unique(tf1[,1])

test.set=dgg

## function to test enrichment using Fisher's exact test
fu=function(db){
    id=tf1[,1]==tttf & !is.na(tf1[,db])
    out=c(NA,NA)
    if(sum(id)>0){
        tgs=unique(tf1[id,2])
        bgg=tdb[[db-2]][!(tdb[[db-2]] %in% tgs)]
        or.p=unlist(fisher.test(matrix(c(sum(test.set %in% tgs),sum(test.set %in% bgg),
                                         length(tgs),length(bgg)),ncol=2))[c("estimate","p.value")])
        out=or.p[1:2]
    }
    out
}

## for all available TFs
library(parallel)
for(tttf in unique(tf1[,1])){
    cat(tttf,"\n")
    o=do.call(cbind,mclapply(3:9,fu,mc.cores=7))
    or.m[tttf,]=o[1,]
    p.m[tttf,]=o[2,]
}

## Analysis of target enrichment in random regions
## relates to results presented in XXX
## 
## generate random regions
## on the same chromosomes as the original deserts
## taking into account centromeres and telomeres
## generate gene lists for each random match set

## original deserts
de=read.table("deserts.bed",as.is=T) # provide file

## telomere/centromere information
tc=read.table("gap.txt",as.is=T)
tc=tc[tc[,8] %in% c("telomere","centromere") & tc[,2] %in% paste0("chr",de[,1]),c(2:4,8)]

tf.es=c()
etf=vector("list",length=100)
for(sim in 1:100){
    cat(sim,"\n")
    co=c()
    ## generating for each of the 5 desert regions a size-matched background region
    ## on the chromosome of the deserts location
    ## excluding centromeric and telomeric regions
    for(chr in de[,1]){
        start=sample(max(tc[tc[,1]==paste0("chr",chr),3]),1)
        if(start<=min(tc[tc[,1]==paste0("chr",chr) & tc[,4]=="telomere",3])) start=min(tc[tc[,1]==paste0("chr",chr) & tc[,4]=="telomere",3])
        if(start>=tc[tc[,1]==paste0("chr",chr) & tc[,4]=="centromere",2] & start<=tc[tc[,1]==paste0("chr",chr) & tc[,4]=="centromere",3]) start=tc[tc[,1]==paste0("chr",chr) & tc[,4]=="centromere",3]
        end=start+de[de[,1]==chr,3]-de[de[,1]==chr,2]
        if(end>=tc[tc[,1]==paste0("chr",chr) & tc[,4]=="centromere",2] & start<tc[tc[,1]==paste0("chr",chr) & tc[,4]=="centromere",2]) end=start+de[de[,1]==chr,3]-de[de[,1]==chr,2]+tc[tc[,1]==paste0("chr",chr) & tc[,4]=="centromere",3]-tc[tc[,1]==paste0("chr",chr) & tc[,4]=="centromere",2]
        if(end>=max(tc[tc[,1]==paste0("chr",chr) & tc[,4]=="telomere",3])) end=min(tc[tc[,1]==paste0("chr",chr) & tc[,4]=="telomere",3])+de[de[,1]==chr,3]-de[de[,1]==chr,2]-(max(tc[tc[,1]==paste0("chr",chr) & tc[,4]=="telomere",3])-start)
        if(end<start){
            co=rbind(co,
                     c(chr,10000,end),
                     c(chr,start,max(tc[tc[,1]==paste0("chr",chr) & tc[,4]=="telomere",3]))
                     )
            
        }
        if(end>start) co=rbind(co,c(chr,start,end))
    }
    write.table(co,col.names=F,row.names=F,quote=F,sep="\t",file="sim_desert_region.bed")

    ## assessing gene content in random regions
    system(paste0("bedtools intersect -a sim_desert_region.bed -b gene_coordinates_hg19_sort.txt -wo -F 1 > sim_region_",sim,"_genes.txt"))

    og=read.table(paste0("sim_region_",sim,"_genes.txt"),as.is=T,sep="\t")

    ## set of protein coding genes located in random region
    dgg=unique(og[,7])

    ## fisher test ORs
    or.m=p.m=matrix(NA,ncol=7,nrow=length(unique(tf1[,1])))
    rownames(or.m)=rownames(p.m)=unique(tf1[,1])
    
    test.set=dgg
    ## calculating enrichment of TF target genes in deserts
    ## compared to the entire genomic background
    for(tttf in unique(tf1[,1])){
        cat(tttf,"\r")
        o=do.call(cbind,mclapply(3:9,fu,mc.cores=7))
        or.m[tttf,]=o[1,]
        p.m[tttf,]=o[2,]
    }
}

