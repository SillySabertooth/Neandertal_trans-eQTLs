# for looking into
# target aSNP chr and position - rs of interest

chr=12 
pos=56740682

wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi

file=ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

# ${chr}_$pos.aSNP.region.chr_pos.txt is list of aSNPs of interest from desired region
# if your vcf have already prefilterred for aSNPs, you can just provide some range around your taget aSNP to look into for R2
# you can also try to get region on the fly
bcftools view -R ${chr}_$pos.aSNP.region.chr_pos.txt $file -Oz -o ${chr}_$pos.for_r2.tmp.vcf.gz
vcftools --gzvcf ${chr}_$pos.for_r2.tmp.vcf.gz --geno-r2 --out ${chr}_$pos.for_r2.tmp
Rscript script.find_R2.R ${chr}_$pos.for_r2.tmp.geno.ld $chr $pos # && rm *tmp* # if working, can uncomment this


# https://samtools.github.io/bcftools/bcftools.html
# http://vcftools.sourceforge.net/man_latest.html
