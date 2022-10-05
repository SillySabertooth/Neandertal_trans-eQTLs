**snippet.haplotype_inferring.sh** includes commands to count R2 among the list of given SNPs (or in the region) using 1KG and bcftools + vcftools, also it has a laucher to run the script <br>

**script.infer_haplotype_based_on_R2.R** has an algorithm for building an archaic haplotype around the target aSNP  <br>
the code was used to infer haplotypes which in turn i) help to gather putative molecular consequences linked to target aSNPs ii) necessary to reject ILS origin of target aSNPs
Mainly, haplotypes were used to elaborate on regulatory potential of a-cTFs aSNPs in the **Discussion**, the 4th paragraph that starts with <br>
*However, our results are primarily based on information from association analyses and computational prediction algorithms.*(c)<br>

the auxiliary file has pairwise R2 linkage information of aSNPs around 500kb of the target example aSNP
