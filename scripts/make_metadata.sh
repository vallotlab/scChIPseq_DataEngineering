#!/bin/bash

## Pacome Prompsy

## Create metadata info
make_metadata() {
odir=$1
prefix=$2
scChIPseq_logfile=$3
flagged_rmPCR_RT_rmDup_count=$4

calc() { awk "BEGIN{print $*}"; }

#Retrieve data from scChIPseq log file:
total_reads=$(grep -e "Number of input reads " $scChIPseq_logfile | sed 's/.*|//g' | grep -o -e '[0-9]*\.*[0-9]*')
avg_input_read_length=$(grep -e "Average input read length " $scChIPseq_logfile | sed 's/.*|//g' | grep -o -e '[0-9]*\.*[0-9]*')
uniquely_mapped=$(grep -e "Uniquely mapped reads number " $scChIPseq_logfile | sed 's/.*|//g' | grep -o -e '[0-9]*\.*[0-9]*')
uniquely_mapped_percent=$(grep -e "Uniquely mapped reads % " $scChIPseq_logfile | sed 's/.*|//g' | grep -o -e '[0-9]*\.*[0-9]*')
avg_mapped_read_length=$(grep -e "Average mapped length " $scChIPseq_logfile | sed 's/.*|//g' | grep -o -e '[0-9]*\.*[0-9]*')
num_splices=$(grep -e "Number of splices: Total " $scChIPseq_logfile | sed 's/.*|//g' | grep -o -e '[0-9]*\.*[0-9]*')
num_annotated_splices=$(grep -e "Number of splices: Annotated (sjdb) " $scChIPseq_logfile | sed 's/.*|//g' | grep -o -e '[0-9]*\.*[0-9]*')
num_GTAG_splices=$(grep -e "Number of splices: GT/AG " $scChIPseq_logfile | sed 's/.*|//g' | grep -o -e '[0-9]*\.*[0-9]*')
num_GCAG_splices=$(grep -e "Number of splices: GC/AG " $scChIPseq_logfile | sed 's/.*|//g' | grep -o -e '[0-9]*\.*[0-9]*')
num_ATAC_splices=$(grep -e "Number of splices: AT/AC " $scChIPseq_logfile | sed 's/.*|//g' | grep -o -e '[0-9]*\.*[0-9]*')
num_noncanonical_splices=$(grep -e "Number of splices: Non-canonical " $scChIPseq_logfile | sed 's/.*|//g' | grep -o -e '[0-9]*\.*[0-9]*')
mismatch_rate=$(grep -e "Mismatch rate per base, % " $scChIPseq_logfile | sed 's/.*|//g' | grep -o -e '[0-9]*\.*[0-9]*')
deletion_rate=$(grep -e "Deletion rate per base " $scChIPseq_logfile | sed 's/.*|//g' | grep -o -e '[0-9]*\.*[0-9]*')
deletion_length=$(grep -e "Deletion average length " $scChIPseq_logfile | sed 's/.*|//g' | grep -o -e '[0-9]*\.*[0-9]*')
insertion_rate=$(grep -e "Insertion rate per base " $scChIPseq_logfile | sed 's/.*|//g' | grep -o -e '[0-9]*\.*[0-9]*')
insertion_length=$(grep -e "Insertion average length " $scChIPseq_logfile | sed 's/.*|//g' | grep -o -e '[0-9]*\.*[0-9]*')
multimapped=$(grep -e "Number of reads mapped to multiple loci " $scChIPseq_logfile | sed 's/.*|//g' | grep -o -e '[0-9]*\.*[0-9]*')
multimapped_percent=$(grep -e "% of reads mapped to multiple loci " $scChIPseq_logfile | sed 's/.*|//g' | grep -o -e '[0-9]*\.*[0-9]*')
multimapped_toomany=$(grep -e "Number of reads mapped to too many loci " $scChIPseq_logfile | sed 's/.*|//g' | grep -o -e '[0-9]*\.*[0-9]*')
multimapped_toomany_percent=$(grep -e "% of reads mapped to too many loci " $scChIPseq_logfile | sed 's/.*|//g' | grep -o -e '[0-9]*\.*[0-9]*')
unmapped_mismatches_percent=$(grep -e "% of reads unmapped: too many mismatches " $scChIPseq_logfile | sed 's/.*|//g' | grep -o -e '[0-9]*\.*[0-9]*')
unmapped_tooshort_percent=$(grep -e "% of reads unmapped: too short " $scChIPseq_logfile | sed 's/.*|//g' | grep -o -e '[0-9]*\.*[0-9]*')
unmapped_other_percent=$(grep -e "% of reads unmapped: other " $scChIPseq_logfile | sed 's/.*|//g' | grep -o -e '[0-9]*\.*[0-9]*')
match_index_1=$(grep -e "## Number of matched indexes 1:" $scChIPseq_logfile | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
match_index_2=$(grep -e "## Number of matched indexes 2:" $scChIPseq_logfile | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
match_index_1_2=$(grep -e "## Number of matched indexes 1 and 2:" $scChIPseq_logfile | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
match_index_3=$(grep -e "## Number of matched indexes 3:" $scChIPseq_logfile | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
match_barcode=$(grep -e "## Number of matched barcodes:" $scChIPseq_logfile | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
uniquely_mapped_and_barcoded=$(grep -e "## Number of reads mapped and barcoded:" $scChIPseq_logfile | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
pcr_duplicates=$(grep -e "## Number of pcr duplicates:" $scChIPseq_logfile | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
rt_duplicates=$(grep -e "## Number of rt duplicates:" $scChIPseq_logfile | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
R1_mapped_R2_unmapped=$(grep -e "## Number of R1 mapped but R2 unmapped:" $scChIPseq_logfile | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
reads_after_pcr_rt_rm=$(grep -e "## Number of reads after PCR and RT removal (not R1 unmapped R2):" $scChIPseq_logfile | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
R2_unmapped_duplicates=$(grep -e "## Number of duplicates:" $scChIPseq_logfile | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
unique_reads=$(grep -e "## Number of reads after duplicates removal:" $scChIPseq_logfile | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')


mkdir -p ${odir}

total_mapped=$( calc $uniquely_mapped + $multimapped + $multimapped_toomany)
unmapped_count=$( calc $total_reads - $total_mapped)
total_unmapped_percent=$( calc $unmapped_mismatches_percent + $unmapped_tooshort_percent + $unmapped_other_percent)
  
uniquely_mapped_unbarcoded=$( calc $uniquely_mapped - $uniquely_mapped_and_barcoded)
multimapped=$( calc $multimapped + $multimapped_toomany)
unmapped=$unmapped_count
              
#Data for the barcode matching graph
reads_after_pcr_rt_rm=$( calc $reads_after_pcr_rt_rm - $R1_mapped_R2_unmapped)
index_1_2_not_3=$( calc $match_index_1_2 - $match_barcode)
index_1_not_2_not_3=$( calc $match_index_1 - $index_1_2_not_3 - $match_barcode)
index_2_not_1_3=$( calc $match_index_2 - $match_index_1_2)
index_3_not_1_2=$( calc $match_index_3 - $match_barcode)
no_index_found=$( calc $total_reads - $match_barcode - $index_1_2_not_3 - $index_1_not_2_not_3 - $index_2_not_1_3 - $index_3_not_1_2)
uniquely_mapped_and_barcoded_percent=$( calc 100*$uniquely_mapped_and_barcoded / $total_reads)
unique_reads_percent=$( calc 100*$unique_reads/ $total_reads)


if [ $total_unmapped_percent -eq 0 ]
then
   unmapped_mismatches=0
   unmapped_other=0
   unmapped_tooshort=0
else
  unmapped_mismatches=$(echo "scale=4;$unmapped_count*$unmapped_mismatches_percent / $total_unmapped_percent" | bc)
  unmapped_tooshort=$(echo "scale=4;$unmapped_count*$unmapped_tooshort_percent / $total_unmapped_percent" | bc)
  unmapped_other=$(echo "scale=4;$unmapped_count*$unmapped_other_percent / $total_unmapped_percent" | bc)
fi


echo $total_mapped 
echo $unmapped_count
echo $total_unmapped_percent
echo $multimapped
echo $unmapped
echo $reads_after_pcr_rt_rm
echo $index_1_2_not_3
echo $index_1_not_2_not_3
echo $index_2_not_1_3
echo $no_index_found
echo $index_3_not_1_2
echo $uniquely_mapped_and_barcoded_percent
echo $unique_reads_percent
echo $unmapped_mismatches
echo $unmapped_tooshort
echo $unmapped_other


echo "Sample_Name,Barcoded,Index 1 and 2 found not 3,Index 1 found not 2 and 3,Index 2 found not 1 and 3,Index 3 found not 1 and 2,No Index Found ~ genomic DNA" > ${odir}/scChIPseq_barcode.csv
echo "$prefix,$match_barcode,$index_1_2_not_3,$index_1_not_2_not_3,$index_2_not_1_3,$index_3_not_1_2,$no_index_found" >> ${odir}/scChIPseq_barcode.csv

echo "Sample_Name,Deduplicated reads, Window duplicates,RT duplicates,PCR duplicates,Uniquely mapped not barcoded,Mapped to multiple loci,Unmapped" > ${odir}/scChIPseq_alignments.csv
echo "$prefix,$unique_reads,$R2_unmapped_duplicates,$rt_duplicates,$pcr_duplicates,$uniquely_mapped_unbarcoded,$multimapped,$unmapped" >> ${odir}/scChIPseq_alignments.csv

n100=$( sed 's/^\s*//g' ${flagged_rmPCR_RT_rmDup_count} | awk -v limit=100 '$1>=limit && NR>1{c++} END{print c+0}')
n500=$( sed 's/^\s*//g' ${flagged_rmPCR_RT_rmDup_count} | awk -v limit=500 '$1>=limit && NR>1{c++} END{print c+0}' )
n1000=$( sed 's/^\s*//g' ${flagged_rmPCR_RT_rmDup_count} | awk -v limit=1000 '$1>=limit && NR>1{c++} END{print c+0}')
n1500=$( sed 's/^\s*//g' ${flagged_rmPCR_RT_rmDup_count} | awk -v limit=1500 '$1>=limit && NR>1{c++} END{print c+0}')

echo "Sample_Name,%Aligned,%Aligned_Barcoded,%Unique_Reads" > ${odir}/scChIPseq_table.csv
echo "$prefix,$uniquely_mapped_percent,$uniquely_mapped_and_barcoded_percent,$unique_reads_percent" >> ${odir}/scChIPseq_table.csv

##Create R report
# Rscript --vanilla ${curdir}/read_stats.R ${odir}/${prefix}_flagged.ReadName_Barcode ${odir} ${prefix}_flagged $MIN_COUNT_PER_BARCODE_BEFORE_RMDUP > ${logdir}/read_stats.Rout 2>&1
# Rscript --vanilla ${curdir}/read_stats.R ${odir}/${prefix}_flagged_rmDup.ReadName_Barcode ${odir} ${prefix}_flagged_rmDup $MIN_COUNT_PER_BARCODE_AFTER_RMDUP > ${logdir}/read_stats_rmdup.Rout 2>&1
# Rscript --vanilla ${curdir}/plot_before_after.R ${odir} ${prefix} > ${logdir}/plot_before_after.Rout 2>&1

desc="Number of barcodes with more than 100 unique reads = $n100 <br>Number of barcodes with more than 500 unique reads = $n500 <br>Number of barcodes with more than 1000 unique reads = $n1000 <br>Number of barcodes with more than 1500 unique reads = $n1500" 
sed -i "s|{desc}|$desc|g" ${odir}/multiqc_config.yaml

echo "make metadata for mqc plots, done !"

}
