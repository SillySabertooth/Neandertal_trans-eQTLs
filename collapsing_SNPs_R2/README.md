**snippet.count_R2_between_SNPs.sh** commands to collapse sets of SNPs in linkage disequilibrium (based on R2) in EUR populations of 1,000Genomes dataset using bcftools + vcftools <br>

**snippet.find_dependent_SNPs_by_R2.R** algorithm for collapsing SNPs into independent sets of variants in linkage disequilibrium <br>

The auxiliary files are i) list of aSNPs to collapse ii) pairwise R2 information for the aSNPs <br> 

The implementation directly relates to the enrichment test in *aSNPs_proportion_enrichment* folder, as well it was used to collapse dependent aSNPs that are linked to a-cTFs <br>
