#!/bin/bash
#Example file to create different configurations and run multiple samples 

get_mem () {
local READ1=$1
file_size_kb=$(du -Lk "$READ1" | cut -f1 | sed 's/[a-Z]*//g' | sed 's/\..*//g')

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
cd ~/GitLab/ChIP-seq_single-cell_LBC_PAIRED_END_3.4/
ASSEMBLY=hg38
OUTPUT_CONFIG=/data/tmp/pprompsy/results/CONFIG_HUMAN_LBC_K27
DESIGN_TYPE=LBC
MARK=h3k27me3

./schip_processing.sh GetConf --template  CONFIG_TEMPLATE --configFile species_design_configs.csv --designType ${DESIGN_TYPE} --genomeAssembly ${ASSEMBLY} --outputConfig ${OUTPUT_CONFIG} --mark ${MARK}

#RUN PIPELINE -> D180
list="MM468_Rcells_5FUR3_K27me3_reseq MM468_Rcells_5FUR3_K27me3_backup"
seq_names=(D180C04 D180C09)
dataset=2009252
n=0
for name in $list
do

NAME=$name

OUTPUT_DIR=/data/kdi_prod/project_result/1184/02.00/results/${ASSEMBLY}/${NAME}
DOWNSTREAM_DIR=/data/kdi_prod/project_result/1184/02.00/results/${ASSEMBLY}/	

READ1=/data/kdi_prod/dataset_all/${dataset}/export/user/${seq_names[$n]}/${seq_names[$n]}.R1.fastq.gz 
READ2=/data/kdi_prod/dataset_all/${dataset}/export/user/${seq_names[$n]}/${seq_names[$n]}.R2.fastq.gz 

mem=$(get_mem ${READ1})

echo $name ${seq_names[$n]} $mem
echo "cd ~/GitLab/ChIP-seq_single-cell_LBC_PAIRED_END_3.4/; ./schip_processing.sh All -f ${READ1} -r ${READ2} -c ${OUTPUT_CONFIG}  -o ${OUTPUT_DIR} --name ${NAME}"  | qsub -l "nodes=1:ppn=10,mem=${mem}gb" -N job_${NAME}
echo "cd ~/GitLab/ChIP-seq_single-cell_LBC_PAIRED_END_3.4/; ./schip_processing.sh All -f ${READ1} -r ${READ2} -c ${OUTPUT_CONFIG}  -o ${OUTPUT_DIR} --name ${NAME}" >> justine_cmds.sh
n=$((n+1))

done



#RUN PIPELINE -> V535
list="MM468_Rcells_5FUR5_K27me3 MM468_Rcells_5FUR5_K27me3_lowcells BC1176_m08_CPT11_K27me3 BC1176_m20_CPT11_K27me3 BC1176_m22_UNT_K27me3"
seq_names=(V535C01 V535C02 V535C03 V535C05 V535C07)
dataset=2009255
n=0
for name in $list
do

NAME=$name

OUTPUT_DIR=/data/kdi_prod/project_result/1184/02.00/results/${ASSEMBLY}/${NAME}
DOWNSTREAM_DIR=/data/kdi_prod/project_result/1184/02.00/results/${ASSEMBLY}/	

READ1=/data/kdi_prod/dataset_all/${dataset}/export/user/${seq_names[$n]}/${seq_names[$n]}.R1.fastq.gz 
READ2=/data/kdi_prod/dataset_all/${dataset}/export/user/${seq_names[$n]}/${seq_names[$n]}.R2.fastq.gz 

mem=$(get_mem ${READ1})

echo $name ${seq_names[$n]} $mem
echo "cd ~/GitLab/ChIP-seq_single-cell_LBC_PAIRED_END_3.4/; ./schip_processing.sh All -f ${READ1} -r ${READ2} -c ${OUTPUT_CONFIG}  -o ${OUTPUT_DIR} --name ${NAME} -s ${DOWNSTREAM_DIR}"  | qsub -l "nodes=1:ppn=10,mem=${mem}gb" -N job_${NAME}
echo "cd ~/GitLab/ChIP-seq_single-cell_LBC_PAIRED_END_3.4/; ./schip_processing.sh All -f ${READ1} -r ${READ2} -c ${OUTPUT_CONFIG}  -o ${OUTPUT_DIR} --name ${NAME} -s ${DOWNSTREAM_DIR}" >> justine_cmds.sh
n=$((n+1))

done

########################################################
#############     MOUSE IP K27 CURIE      ##################
########################################################


#CREATE CONFIG FILE : MOUSE, BEADS CURIE, IP, BED = 20k +50k
cd ~/GitLab/ChIP-seq_single-cell_LBC_PAIRED_END_3.4/
ASSEMBLY=mm10
OUTPUT_CONFIG=/data/tmp/pprompsy/results/CONFIG_MOUSE_LBC_K27
DESIGN_TYPE=LBC
MARK=h3k27me3

./schip_processing.sh GetConf --template  CONFIG_TEMPLATE --configFile species_design_configs.csv --designType ${DESIGN_TYPE} --genomeAssembly ${ASSEMBLY} --outputConfig ${OUTPUT_CONFIG} --mark ${MARK}

#RUN PIPELINE -> D180
list="BC1176_m08_CPT11_K27me3 BC1176_m20_CPT11_K27me3 BC1176_m22_UNT_K27me3"
seq_names=(V535C03 V535C05 V535C07)
dataset=2009255
n=0
for name in $list
do

NAME=$name

OUTPUT_DIR=/data/kdi_prod/project_result/1184/02.00/results/${ASSEMBLY}/${NAME}
DOWNSTREAM_DIR=/data/kdi_prod/project_result/1184/02.00/results/${ASSEMBLY}/	

READ1=/data/kdi_prod/dataset_all/${dataset}/export/user/${seq_names[$n]}/${seq_names[$n]}.R1.fastq.gz 
READ2=/data/kdi_prod/dataset_all/${dataset}/export/user/${seq_names[$n]}/${seq_names[$n]}.R2.fastq.gz 

mem=$(get_mem ${READ1})

echo $name ${seq_names[$n]} $mem
echo "cd ~/GitLab/ChIP-seq_single-cell_LBC_PAIRED_END_3.4/; ./schip_processing.sh All -f ${READ1} -r ${READ2} -c ${OUTPUT_CONFIG}  -o ${OUTPUT_DIR} --name ${NAME} -s ${DOWNSTREAM_DIR}"  | qsub -l "nodes=1:ppn=10,mem=${mem}gb" -N job_${NAME}
echo "cd ~/GitLab/ChIP-seq_single-cell_LBC_PAIRED_END_3.4/; ./schip_processing.sh All -f ${READ1} -r ${READ2} -c ${OUTPUT_CONFIG}  -o ${OUTPUT_DIR} --name ${NAME} -s ${DOWNSTREAM_DIR}" >> justine_cmds.sh
n=$((n+1))

done

########################################################
#############     MOUSE IP K4 CURIE      ##################
########################################################


#CREATE CONFIG FILE : MOUSE, BEADS CURIE, IP, BED = 20k +50k
cd ~/GitLab/ChIP-seq_single-cell_LBC_PAIRED_END_3.4/
ASSEMBLY=mm10
OUTPUT_CONFIG=/data/tmp/pprompsy/results/CONFIG_MOUSE_LBC_K4
DESIGN_TYPE=LBC
MARK=h3k4me3

./schip_processing.sh GetConf --template  CONFIG_TEMPLATE --configFile species_design_configs.csv --designType ${DESIGN_TYPE} --genomeAssembly ${ASSEMBLY} --outputConfig ${OUTPUT_CONFIG} --mark ${MARK}

#RUN PIPELINE -> D180
list="BC1176_m08_CPT11_K4me3 BC1176_m20_CPT11_K4me3"
seq_names=(V535C04 V535C06)
dataset=2009255
n=0
for name in $list
do

NAME=$name

OUTPUT_DIR=/data/kdi_prod/project_result/1184/02.00/results/${ASSEMBLY}/${NAME}
DOWNSTREAM_DIR=/data/kdi_prod/project_result/1184/02.00/results/${ASSEMBLY}/	

READ1=/data/kdi_prod/dataset_all/${dataset}/export/user/${seq_names[$n]}/${seq_names[$n]}.R1.fastq.gz 
READ2=/data/kdi_prod/dataset_all/${dataset}/export/user/${seq_names[$n]}/${seq_names[$n]}.R2.fastq.gz 

mem=$(get_mem ${READ1})

echo $name ${seq_names[$n]} $mem
echo "cd ~/GitLab/ChIP-seq_single-cell_LBC_PAIRED_END_3.4/; ./schip_processing.sh All -f ${READ1} -r ${READ2} -c ${OUTPUT_CONFIG}  -o ${OUTPUT_DIR} --name ${NAME} -s ${DOWNSTREAM_DIR}"  | qsub -l "nodes=1:ppn=10,mem=${mem}gb" -N job_${NAME}
echo "cd ~/GitLab/ChIP-seq_single-cell_LBC_PAIRED_END_3.4/; ./schip_processing.sh All -f ${READ1} -r ${READ2} -c ${OUTPUT_CONFIG}  -o ${OUTPUT_DIR} --name ${NAME} -s ${DOWNSTREAM_DIR}" >> justine_cmds.sh
n=$((n+1))

done

########################################################
#############     HUMAN IP K4 CURIE      ##################
########################################################


#CREATE CONFIG FILE : MOUSE, BEADS CURIE, IP, BED = 20k +50k
cd ~/GitLab/ChIP-seq_single-cell_LBC_PAIRED_END_3.4/
ASSEMBLY=hg38
OUTPUT_CONFIG=/data/tmp/pprompsy/results/CONFIG_HUMAN_LBC_K4
DESIGN_TYPE=LBC
MARK=h3k4me3

./schip_processing.sh GetConf --template  CONFIG_TEMPLATE --configFile species_design_configs.csv --designType ${DESIGN_TYPE} --genomeAssembly ${ASSEMBLY} --outputConfig ${OUTPUT_CONFIG} --mark ${MARK}

#RUN PIPELINE
list="BC1176_m08_CPT11_K4me3 BC1176_m20_CPT11_K4me3"
seq_names=(V535C04 V535C06)
dataset=2009255
n=0
for name in $list
do

NAME=$name

OUTPUT_DIR=/data/kdi_prod/project_result/1184/02.00/results/${ASSEMBLY}/${NAME}
DOWNSTREAM_DIR=/data/kdi_prod/project_result/1184/02.00/results/${ASSEMBLY}/	

READ1=/data/kdi_prod/dataset_all/${dataset}/export/user/${seq_names[$n]}/${seq_names[$n]}.R1.fastq.gz 
READ2=/data/kdi_prod/dataset_all/${dataset}/export/user/${seq_names[$n]}/${seq_names[$n]}.R2.fastq.gz 

mem=$(get_mem ${READ1})

echo $name ${seq_names[$n]} $mem
echo "cd ~/GitLab/ChIP-seq_single-cell_LBC_PAIRED_END_3.4/; ./schip_processing.sh All -f ${READ1} -r ${READ2} -c ${OUTPUT_CONFIG}  -o ${OUTPUT_DIR} --name ${NAME} -s ${DOWNSTREAM_DIR}"  | qsub -l "nodes=1:ppn=10,mem=${mem}gb" -N job_${NAME}
echo "cd ~/GitLab/ChIP-seq_single-cell_LBC_PAIRED_END_3.4/; ./schip_processing.sh All -f ${READ1} -r ${READ2} -c ${OUTPUT_CONFIG}  -o ${OUTPUT_DIR} --name ${NAME} -s ${DOWNSTREAM_DIR}" >> justine_cmds.sh
n=$((n+1))

done


########################################################
#############     HUMAN IP K27 HIFIBIO      ##################
########################################################


#CREATE CONFIG FILE : HUMAN, BEADS HIFIBIO, IP, BED = 20k +50k
cd ~/GitLab/ChIP-seq_single-cell_LBC_PAIRED_END_3.4/
ASSEMBLY=hg38
OUTPUT_CONFIG=/data/tmp/pprompsy/results/CONFIG_HUMAN_Hifibio_K27
DESIGN_TYPE=Hifibio
MARK=h3k27me3

./schip_processing.sh GetConf --template  CONFIG_TEMPLATE --configFile species_design_configs.csv --designType ${DESIGN_TYPE} --genomeAssembly ${ASSEMBLY} --outputConfig ${OUTPUT_CONFIG} --mark ${MARK}

#RUN PIPELINE
list="HBCx95_m00_UNT_K27me3_reseq HBCx95_m40_CAPAR_K27me3_reseq HBCx22_m99_TAMR_K27me3_reseq HBCx95_m00_UNT_K27me3_backup HBCx95_m40_CAPAR_K27me3_backup HBCx22_m99_TAMR_K27me3_backup HBCx95_m43_UNT_K27me3_backup"
seq_names=(D180C01 D180C02 D180C03 D180C05 D180C06 D180C07 D180C08)
dataset=2009252
n=0
for name in $list
do

NAME=$name

OUTPUT_DIR=/data/kdi_prod/project_result/1184/02.00/results/${ASSEMBLY}/${NAME}
DOWNSTREAM_DIR=/data/kdi_prod/project_result/1184/02.00/results/${ASSEMBLY}/	

READ1=/data/kdi_prod/dataset_all/${dataset}/export/user/${seq_names[$n]}/${seq_names[$n]}.R1.fastq.gz 
READ2=/data/kdi_prod/dataset_all/${dataset}/export/user/${seq_names[$n]}/${seq_names[$n]}.R2.fastq.gz 

mem=$(get_mem ${READ1})

echo $name ${seq_names[$n]} $mem
echo "cd ~/GitLab/ChIP-seq_single-cell_LBC_PAIRED_END_3.4/; ./schip_processing.sh All -f ${READ1} -r ${READ2} -c ${OUTPUT_CONFIG}  -o ${OUTPUT_DIR} --name ${NAME}"  | qsub -l "nodes=1:ppn=10,mem=${mem}gb" -N job_${NAME}
echo "cd ~/GitLab/ChIP-seq_single-cell_LBC_PAIRED_END_3.4/; ./schip_processing.sh All -f ${READ1} -r ${READ2} -c ${OUTPUT_CONFIG}  -o ${OUTPUT_DIR} --name ${NAME}" >> justine_cmds.sh
n=$((n+1))

done

######################################################################################################################
######################################################################################################################

########################################################
#############        MOUSE IP K27 HIFIBIO          ##################
########################################################


#CREATE CONFIG FILE : MOUSE, BEADS CURIE, IP, BED = 20k +50k
cd ~/GitLab/ChIP-seq_single-cell_LBC_PAIRED_END_3.4/
ASSEMBLY=mm10
OUTPUT_CONFIG=/data/tmp/pprompsy/results/CONFIG_MOUSE_Hifibio_K27
DESIGN_TYPE=Hifibio
MARK=h3k27me3

./schip_processing.sh GetConf --template  CONFIG_TEMPLATE --configFile species_design_configs.csv --designType ${DESIGN_TYPE} --genomeAssembly ${ASSEMBLY} --outputConfig ${OUTPUT_CONFIG}  --mark ${MARK}

#RUN PIPELINE
list="HBCx95_m00_UNT_K27me3_reseq HBCx95_m40_CAPAR_K27me3_reseq HBCx22_m99_TAMR_K27me3_reseq HBCx95_m00_UNT_K27me3_backup HBCx95_m40_CAPAR_K27me3_backup HBCx22_m99_TAMR_K27me3_backup HBCx95_m43_UNT_K27me3_backup"
seq_names=(D180C01 D180C02 D180C03 D180C05 D180C06 D180C07 D180C08)
dataset=2009252
n=0
for name in $list
do

NAME=$name

OUTPUT_DIR=/data/kdi_prod/project_result/1184/02.00/results/${ASSEMBLY}/${NAME}
DOWNSTREAM_DIR=/data/kdi_prod/project_result/1184/02.00/results/${ASSEMBLY}/	


READ1=/data/kdi_prod/dataset_all/${dataset}/export/user/${seq_names[$n]}/${seq_names[$n]}.R1.fastq.gz 
READ2=/data/kdi_prod/dataset_all/${dataset}/export/user/${seq_names[$n]}/${seq_names[$n]}.R2.fastq.gz 

mem=$(get_mem ${READ1})

echo $name ${seq_names[$n]} $mem
echo "cd ~/GitLab/ChIP-seq_single-cell_LBC_PAIRED_END_3.4/; ./schip_processing.sh All -f ${READ1} -r ${READ2} -c ${OUTPUT_CONFIG}  -o ${OUTPUT_DIR} --name ${NAME} " | qsub -l "nodes=1:ppn=10,mem=${mem}gb" -N job_${NAME}
echo "cd ~/GitLab/ChIP-seq_single-cell_LBC_PAIRED_END_3.4/; ./schip_processing.sh All -f ${READ1} -r ${READ2} -c ${OUTPUT_CONFIG}  -o ${OUTPUT_DIR} --name ${NAME} " >> justine_cmds.sh
n=$((n+1))

done


######################################################################################################################
######################################################################################################################

########################################################
#############        HUMAN UNBOUND CURIE          ##################
########################################################

#CREATE CONFIG FILE : MOUSE, BEADS CURIE, IP, BED = 20k +50k
cd ~/GitLab/ChIP-seq_single-cell_LBC_PAIRED_END_3.4/
ASSEMBLY=hg38
OUTPUT_CONFIG=/data/tmp/pprompsy/results/CONFIG_HUMAN_LBC_UNBOUND
MARK=unbound
DESIGN_TYPE=LBC

./schip_processing.sh GetConf --template  CONFIG_TEMPLATE --configFile species_design_configs.csv --designType ${DESIGN_TYPE} --genomeAssembly ${ASSEMBLY} --outputConfig ${OUTPUT_CONFIG} --mark ${MARK}

#RUN PIPELINE
list="Mix_UB1_RT_N6 Mix_UB2_RT_Top Mix_UB3_RT_Top_BotRead1 Mix_UB4_RT_TP53Read1"
seq_names=(V536C16 V536C17 V536C18 V536C19)
dataset=2009484
n=0
for name in $list
do

NAME=$name

OUTPUT_DIR=/data/kdi_prod/project_result/1184/02.00/results/${ASSEMBLY}/${NAME}
DOWNSTREAM_DIR=/data/kdi_prod/project_result/1184/02.00/results/${ASSEMBLY}/	

READ1=/data/kdi_prod/dataset_all/${dataset}/export/user/${seq_names[$n]}/${seq_names[$n]}.R1.fastq.gz 
READ2=/data/kdi_prod/dataset_all/${dataset}/export/user/${seq_names[$n]}/${seq_names[$n]}.R2.fastq.gz 

mem=$(get_mem ${READ1})

echo $name ${seq_names[$n]} $mem
echo "cd ~/GitLab/ChIP-seq_single-cell_LBC_PAIRED_END_3.4/; ./schip_processing.sh All -f ${READ1} -r ${READ2} -c ${OUTPUT_CONFIG}  -o ${OUTPUT_DIR} --name ${NAME}" | qsub -l "nodes=1:ppn=10,mem=${mem}gb" -N job_${NAME}
echo "cd ~/GitLab/ChIP-seq_single-cell_LBC_PAIRED_END_3.4/; ./schip_processing.sh All -f ${READ1} -r ${READ2} -c ${OUTPUT_CONFIG}  -o ${OUTPUT_DIR} --name ${NAME}" >> justine_cmds.sh
n=$((n+1))

done

echo "Sent all jobs !"
