
###### one cat get all chr from 1kG using the links
#for chr in {1..22}; do
#wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
#wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
#done
######

#path=~/genomes/ #prove some path to 1kG files
file=cTFs_GTEx.aSNPs.to_collapse.tsv #

tail -n+2 $file | sort -k1,1 -k2,2 > variants.chr_pos # in the resulting file it should be chr and pos

#chr=1 #debugging
#subset all the chromosomes
for chr in {1..22}; do
bcftools view -R variants.chr_pos -S ../EUR.samples -v snps \
  $path/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -Oz -o $chr.vcf.gz
done
  
bcftools concat {1..22}.vcf.gz -Oz -o all.vcf.gz && rm {1..22}.vcf.gz # concat and remove

vcftools --gzvcf all.vcf.gz --geno-r2 --out background # count R2

# https://samtools.github.io/bcftools/bcftools.html
# http://vcftools.sourceforge.net/man_latest.html
