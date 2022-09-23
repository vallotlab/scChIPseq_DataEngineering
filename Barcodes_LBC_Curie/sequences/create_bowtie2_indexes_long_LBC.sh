#!/bin/bash
annot_path=/home/pprompsy/Documents/Data/Annotation/Barcodes_LBC_Curie
#Index B --> Column 3 (Forward) - link left = TGAC
cut -f 3 ${annot_path}/LBC_Barcode_Oligos.csv | awk 'NR>1{if(NR<11){print ">B0"NR-1"\n"$0"TGAC"}else{print ">B"NR-1"\n"$0"TGAC"}}' > ${annot_path}/index_B_long.fasta

#Index C --> Column 5 (Forward) - link left = ACCA
cut -f 5 ${annot_path}/LBC_Barcode_Oligos.csv | awk 'NR>1{if(NR<11){print ">C0"NR-1"\n"$0"ACCA"}else{print ">C"NR-1"\n"$0"ACCA"}}' > ${annot_path}/index_C_long.fasta

#Index D --> Column 7 (Forward) - link left = CAAC
cut -f 7 ${annot_path}/LBC_Barcode_Oligos.csv | awk 'NR>1{if(NR<11){print ">D0"NR-1"\n"$0"CAAC"}else{print ">D"NR-1"\n"$0"CAAC"}}' > ${annot_path}/index_D_long.fasta

#Build indexes naming them ref_index_x
bowtie2_dir_path=${annot_path}/bowtie_2_index_long
mkdir -p ${bowtie2_dir_path}

bowtie2-build ${annot_path}/index_B_long.fasta ${bowtie2_dir_path}/ref_index_1
bowtie2-build ${annot_path}/index_C_long.fasta ${bowtie2_dir_path}/ref_index_2
bowtie2-build ${annot_path}/index_D_long.fasta ${bowtie2_dir_path}/ref_index_3
