####################################
####################################
## correlation between Z-scores (eQTLGen) and slopes (GTEx) for shared cTFs
## cis-eQTLs in eQTLGen with aSNP associations
eg1=read.table("aSNP_cisreg_TF_FDR_441.tsv",as.is=T,head=T) # Table S11

aa=read.table("GTEx_aSNP_TFs_cis-eQTLs_q0.05_new.tab",as.is=T,head=T) # Table S12

## overlapping cTFs between eQTLGen and GTEx
ooo=sort(intersect(aa[,21],eg1[,1]))

## function to select largest abs(Z-score) per cTF
di=function(x){
    y=eg1[eg1$GeneSymbol==x,,drop=F]
    y=y[which.max(abs(y[,"Zscore"])),,drop=F]
    na=ifelse(y$YRI=="0/0",y[,4],y[,3])
    if(y$AssessedAllele==na) out=y$Zscore else out=y$Zscore*(-1)
    out
}

eqtlgen.tf=sort(sapply(unique(eg1$GeneSymbol),di))

## generate table with Z-scores (eQTLGen) and slopes (GTEx) for shared cTFs between both databases
cor.beta=c()
for(g in ooo){
    cor.beta=rbind(cor.beta,cbind(rep(eqtlgen.tf[g],sum(aa[,21]==g)),-aa[aa[,21]==g,"slope"],aa[aa[,21]==g,20]=="Whole_Blood"))
}
