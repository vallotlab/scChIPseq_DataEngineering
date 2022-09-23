#!/bin/bash
annot_path=/data/users/pprompsy/Annotation/Barcodes_10X_bowtie2_737K-cratac-v1

#Build indexes naming them ref_index_x
bowtie2_dir_path=${annot_path}/bowtie_2_index/

/bioinfo/local/build/Centos/bowtie2/bowtie2-2.2.9/bowtie2-build ${annot_path}/FASTA/737K-cratac-v1.fasta ${bowtie2_dir_path}/ref_index_1

