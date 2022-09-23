#!/bin/bash
annot_path=/home/pprompsy/Documents/Data/Annotation/Barcodes_HiFiBio/sequences/

#Index B -->  - link left = TGAC
 awk '{split($1,a,"_"); if(a[2]=="B"){print $2}}' ${annot_path}/AllPossibleBarcodes.txt | tr -d '\r'  | awk '{if(NR<10){print ">B0"NR"\n"$0"TGAC"}else{print ">B"NR"\n"$1"TGAC"}}'  > ${annot_path}/index_B_long.fasta

#Index C --> Column 5 (Forward) - link left = TCCC
awk '{split($1,a,"_"); if(a[2]=="C"){print $2}}' ${annot_path}/AllPossibleBarcodes.txt | tr -d '\r'  | awk '{if(NR<10){print ">C0"NR"\n"$0"TCCC"}else{print ">C"NR"\n"$1"TCCC"}}'  >${annot_path}/index_C_long.fasta

#Index D --> Column 7 (Forward) - link left = CAAC
awk '{split($1,a,"_"); if(a[2]=="D"){print $2}}' ${annot_path}/AllPossibleBarcodes.txt | tr -d '\r'  | awk '{if(NR<10){print ">D0"NR"\n"$0"CAAC"}else{print ">D"NR"\n"$1"CAAC"}}'  >${annot_path}/index_D_long.fasta

#Build indexes naming them ref_index_x

bowtie2_dir_path=${annot_path}/../index_barcode_bowtie2/bowtie_2_index_long
mkdir -p ${bowtie2_dir_path}

bowtie2-build ${annot_path}/index_B_long.fasta ${bowtie2_dir_path}/ref_index_1
bowtie2-build ${annot_path}/index_C_long.fasta ${bowtie2_dir_path}/ref_index_2
bowtie2-build ${annot_path}/index_D_long.fasta ${bowtie2_dir_path}/ref_index_3
