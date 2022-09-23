#!/bin/bash
sample_sheet=$1

get_mem () {
local READ1=$1
file_size_kb=$(du -Lk "$READ1" | cut -f1 | sed 's/[a-zA-Z]*//g' | sed 's/\..*//g')


if [ $file_size_kb == "" ]
then
file_size_kb=0
fi

if [ $file_size_kb -le 1000000 ]
then
mem=34
elif [ $file_size_kb -le 5000000 ]
then
mem=40
elif [  $file_size_kb -le 10000000 ]
then
mem=60
elif [ $file_size_kb -le 15000000 ]
then
mem=80
elif [  $file_size_kb -le 20000000 ]
then
mem=100
elif [ $file_size_kb -le 25000000 ]
then
mem=120
else 
mem=120
fi
echo "$mem"

}

########################################################
#############     HUMAN IP K27 CURIE      ##################
########################################################


#CREATE CONFIG FILE : HUMAN, BEADS CURIE, IP, BED = 20k +50k
cd ~/GitLab/ChIP-seq_single-cell_LBC_4.0/

while IFS= read -r line
do
  DATASET_NAME=$(echo "$line" | cut -d',' -f1)
  DATASET_NUMBER=$(echo "$line" | cut -d',' -f2)
  FINAL_NAME=$(echo "$line" | cut -d',' -f3)
  ASSEMBLY=$(echo "$line" | cut -d',' -f4)
  MARK=$(echo "$line" | cut -d',' -f5)
  if [[ $MARK == "h3k27me3" && $ASSEMBLY == "hg38" ]]
  then
  OUTPUT_CONFIG=/data/tmp/pprompsy/results/CONFIG_HUMAN_LBC_K27
  fi
  if [[ $MARK == "h3k4me3" && $ASSEMBLY == "hg38" ]]
  then
  OUTPUT_CONFIG=/data/tmp/pprompsy/results/CONFIG_HUMAN_LBC_K4
  fi
  if [[ $MARK == "unbound" && $ASSEMBLY == "hg38" ]]
  then
  OUTPUT_CONFIG=/data/tmp/pprompsy/results/CONFIG_HUMAN_LBC_UNBOUND
  fi
  DESIGN_TYPE=LBC
  if [[ $MARK == "h3k27me3" && $ASSEMBLY == "mm10" ]]
  then
  OUTPUT_CONFIG=/data/tmp/pprompsy/results/CONFIG_MOUSE_LBC_H3K27ME3
  fi
  if [[ $MARK == "h3k4me3" && $ASSEMBLY == "mm10" ]]
  then
  OUTPUT_CONFIG=/data/tmp/pprompsy/results/CONFIG_MOUSE_LBC_H3K4ME3
  fi

  #./schip_processing.sh GetConf --template  CONFIG_TEMPLATE --configFile species_design_configs.csv --designType ${DESIGN_TYPE} --genomeAssembly ${ASSEMBLY} --outputConfig ${OUTPUT_CONFIG} --mark ${MARK}
 
  OUTPUT_DIR=/data/kdi_prod/project_result/1184/02.00/results/${ASSEMBLY}/${FINAL_NAME}
  DOWNSTREAM_DIR=/data/kdi_prod/project_result/1184/02.00/results/${ASSEMBLY}/
  
  READ1=/data/tmp/pprompsy/tmp.R1.fastq.gz
  READ2=/data/tmp/pprompsy/tmp.R2.fastq.gz
  touch $READ1
  touch $READ2
  echo "cd ~/GitLab/ChIP-seq_single-cell_LBC_4.0/; ./schip_processing.sh Counting -f ${READ1} -r ${READ2} -c ${OUTPUT_CONFIG}  -o ${OUTPUT_DIR} --name ${FINAL_NAME}" | qsub -l "nodes=1:ppn=8,mem=20gb" -N job_${FINAL_NAME}_${ASSEMBLY}_counting

done < "$sample_sheet"


