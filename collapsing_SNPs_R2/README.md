**snippet.count_R2_between_SNPs.sh** includes commands to count R2 among the list of given SNPs using 1KG and bcftools + vcftools <br>

**snippet.find_dependent_SNPs_by_R2.R** has an algorithm for collapsing SNPs into independent sets of interlinked variants <br>

The auxiliary files are i) list of aSNPs to collapse ii) pairwise R2 information for the aSNPs <br> 

The implementation directly relates to the enrichment test in *aSNPs_proportion_enrichment* folder, as well it was used to collapse dependent aSNPs that are linked to a-cTFs (the folder example) <br>
