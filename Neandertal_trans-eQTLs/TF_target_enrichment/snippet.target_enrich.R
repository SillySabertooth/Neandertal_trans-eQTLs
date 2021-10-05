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

####################################
## test whether target trans-eQTLs gene sets are enriched for TFs

## table of significant trans-eQTLs from https://www.eqtlgen.org/
a=fread("2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz",data.table=F,sep="\t")

## find associations with aSNPs
p.asnp=o[o %in% c(paste0(a$SNPChr,"_",a$SNPPos,"_",a$AssessedAllele,"_",a$OtherAllele),
                  paste0(a$SNPChr,"_",a$SNPPos,"_",a$OtherAllele,"_",a$AssessedAllele))]

tes=lapply(p.asnp,function(x)
    a[paste0(a$SNPChr,"_",a$SNPPos,"_",a$AssessedAllele,"_",a$OtherAllele)==x | paste0(a$SNPChr,"_",a$SNPPos,"_",a$OtherAllele,"_",a$AssessedAllele)==x,,drop=F])

## definition of two hub gene sets
hub1=unique(unlist(sapply(tes[3:4],function(x) x[,"GeneSymbol"])))
hub2=tes[[14]][,"GeneSymbol"]

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

test.set=dgg # enrichment of desert genes
## test.set=hub1 # enrichment of trans-eQTL genes with link to aSNPs on chr3
## test.set=hub2 # enrichment of trans-eQTL genes with link to aSNPs on chr18

## for all available TFs
for(tttf in unique(tf1[,1])){
    cat(tttf,"\n")
    ## for all prediction databases
    for(db in 3:9){
        id=tf1[,1]==tttf & !is.na(tf1[,db])
        if(sum(id)>0){
            tgs=unique(tf1[id,2])
            bgg=tdb[[db-2]][!(tdb[[db-2]] %in% tgs)]
            or.p=unlist(fisher.test(matrix(c(sum(test.set %in% tgs),sum(hub1 %in% bgg),
                                             length(tgs),length(bgg)),ncol=2))[c("estimate","p.value")])
            or.m[tttf,db-2]=or.p[1]
            p.m[tttf,db-2]=or.p[2]
        }
    }
}
