The folder consists of the following files:  <br>
1. **snippet.create_TFTG_DB.plusMotifmap.sh** - mutating of Human TFTG database  <br>
2. **snippet.process_TFTG.R** - processing of the obtained table  <br>
3. **snippet.create_Harmonizome.R** - merging Harmonizome databases  <br>
4. **workflow.normalising_TFDBs.R** - the main long workflow of normalizing all DBs and the upfollowing merging  <br>

The auxiliary files are raw and intermediate databases, Ensembl BioMart gene name / ID table, a list of Harmonizome databases to unite, and a function made to normalize gene names regarding synonyms <br>

*TF_database_regGenes_v1.tsv.gz* is the resulting database that was used to analyze the reach of a-cTF in the Neandertal deserts, as well as to define TFs. In general, the database was used throughout the **Results** section starting from *Prediction of Neandertal-linked trans-eQTL effects* part
