# here is the code to infer a haplotype around the given aSNP

chr=12 
pos=56740682

#wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
#wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi

#file=ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

#### ${chr}_$pos.aSNP.region.chr_pos.txt #### is list of aSNPs of interest from desired region
# if your vcf have already prefiltered for aSNPs, you can just provide some range around your target aSNP to count R2 (as we do)
# you can also try to get the region (or SNP list) on the fly without downloading files from the 1kg server,
# for that, please change the file name to the file url



# as a toy dataset, we provide in this folder ld map of archaic SNPs in 500kb around the given aSNP 

#bcftools view -R ${chr}_$pos.aSNP.region.chr_pos.txt $file -Oz -o ${chr}_$pos.for_r2.tmp.vcf.gz
#vcftools --gzvcf ${chr}_$pos.for_r2.tmp.vcf.gz --geno-r2 --out ${chr}_$pos.for_r2.tmp
Rscript script.infer_haplotype_based_on_R2.R ${chr}_$pos.for_r2.tmp.geno.ld $chr $pos


# https://samtools.github.io/bcftools/bcftools.html
# http://vcftools.sourceforge.net/man_latest.html
