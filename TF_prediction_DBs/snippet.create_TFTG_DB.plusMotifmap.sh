#download data from here http://tfbsdb.systemsbiology.net/download

# I guess this only works when you launch the code as a script ">bash snippet.process_TFTG_DB.plusMotifmap.sh" 
# previously commented the rest of the code
exec &>> logfile2.txt


# there will be some losses that should be treated semi-manually; you can append them to the resulting table
tail -n2 humanTFMotifEntrezMappings.csv | while read line; do
gene_TF=`echo $line | cut -f1 -d ","`
doc_id=`echo $line | cut -f3 -d ","`
cat geneHits_minus5Kbp_plus5Kbp/fullGenome_motifHits_$doc_id.csv | while read i; do
echo $gene_TF `echo $i | cut -f1 -d","` >> tftg_tf.tsv 
done
done

# you need to download all DBs that are mentioned in the df_list as gene set library, rename them accordingly with df_list names
# https://maayanlab.cloud/Harmonizome/

# a command for adjusting MotifMap_Pred_TFtargets.gmt.gz
# gunzip -c MotifMap_Pred_TFtargets.gmt.gz | awk '{ print length, $0 }' | sort -n -s -r | \
#  cut -d" " -f2- | bgzip -c > MotifMap_Pred_TFtargets_adj.gmt.gz
