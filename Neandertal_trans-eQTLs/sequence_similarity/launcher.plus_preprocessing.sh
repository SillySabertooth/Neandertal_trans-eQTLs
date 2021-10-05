chr=12
pos=56740682
geno="0|0"

name=56627300_56753822
reg=`echo $name | sed "s/_/-/"`

# you can get Neanderthals from any place where do you want
# you can download them directly and then process or subset on the fly (below), it can take time
# for example http://cdna.eva.mpg.de/neandertal/Vindija/VCF/
# remove index since it has the same name acros samples

# download all + merge Neandertals + index 
bcftools view -r $chr:$reg \
   http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/chr${chr}_mq25_mapab100.vcf.gz \
   -Oz -o N.Vinjia.$chr.$name.vcf.gz && rm chr${chr}_mq25_mapab100.vcf.gz.tbi   

bcftools view -r $chr:$reg \
   http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Altai/chr${chr}_mq25_mapab100.vcf.gz \
   -Oz -o N.Altai.$chr.$name.vcf.gz && rm chr${chr}_mq25_mapab100.vcf.gz.tbi
   
bcftools view -r $chr:$reg \
   http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Denisova/chr${chr}_mq25_mapab100.vcf.gz \
   -Oz -o N.Denisova.$chr.$name.vcf.gz && rm chr${chr}_mq25_mapab100.vcf.gz.tbi   

bcftools view -r $chr:$reg \
  http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
  -Oz -o KG.$chr.$name.vcf.gz  
bcftools index -t KG.$chr.$name.vcf.gz
for i in N*.$chr.$name.vcf.gz; do bcftools index -t $i; done

bcftools merge N.Vinjia.$chr.$name.vcf.gz N.Altai.$chr.$name.vcf.gz N.Denisova.$chr.$name.vcf.gz |\
  bcftools annotate -x FORMAT -Oz -o Neand.$chr.$name.vcf.gz && bcftools index -t Neand.$chr.$name.vcf.gz
  
bcftools merge Neand.$chr.$name.vcf.gz KG.$chr.$name.vcf.gz -Oz -o Neand.KG.$chr.$name.vcf.gz &&\
  bcftools index -t Neand.KG.$chr.$name.vcf.gz 
  
  
  
  
# get a list for subset

bcftools query -l Neand.$chr.$name.vcf.gz > ${chr}_$pos.samples # Neandertals first - that's important!!
# in the script, 1KG starts from 13 columns, preassuming that the first 3 are Neanderthal...
   
# get samples that have Neand allele(-s) in the position 
bcftools view -r $chr:$pos KG.$chr.$name.vcf.gz | bcftools query -f '[%CHROM:%POS:%REF:%ALT %SAMPLE %GT\n]' |\
  grep -v $geno | cut -f2 -d " " >> ${chr}_$pos.samples

# extract random two times more samples without Neand allele  
num=$(((`cat ${chr}_$pos.samples | wc -l`-3) *2)) 
bcftools query -l KG.$chr.$name.vcf.gz | grep -v -f ${chr}_$pos.samples |\
  shuf | head -n $num >> ${chr}_$pos.samples # sample with no haplo
  
cat ${chr}_$pos.samples | wc -l

  

# make the final vcf with 3 Neand, no multi positions, subset 1kg samples; run the script  
# yes, I might subset samples before this step, but I don't

bcftools view -r $chr:$reg -S ${chr}_$pos.samples -m2 -M2 -v snps Neand.KG.$chr.$name.vcf.gz -Oz -o ${chr}.$name.final.vcf.gz

# it can work a bit long in the loop, I suggest to remove it from here and feed to the "parallel" command on your own  
Rscript script.plot_seq_comp.R ${chr}.$name.final.vcf.gz $pos $geno  
