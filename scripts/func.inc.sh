## Mapping Research Pipeline
## Copyleft 2018 Institut Curie
## Author(s): Nicolas Servant, Pacôme Prompsy
## Contact: nicolas.servant@curie.fr; pacome.prompsy@curie.fr
## This software is distributed without any guarantee under the terms of the CECILL License
## See the LICENCE file for details


##GetConf function <- from EWOK (https://gitlab.curie.fr/data-analysis/ewok)

# getConf_func(){
# if [[ $# != "6" && $# != "7" ]]; then echo "All arguments must be filled : getConf_func TEMPLATE CONFIGS DESIGN_TYPE ASSEMBLY OUT_FILE MARK [TARGET_BED] "; echo; help_func GetConf;
# else 	
# 	local TEMPLATE=${1}
# 	local CONFIGS=${2}
# 	local DESIGN_TYPE=${3}
# 	local ASSEMBLY=${4}
# 	local TYPE=${5}
# 	local OUT_FILE=${6}
# 	local TARGET_BED=($(echo ${7} | sed 's/,/ /g'))
# 
# 	mkdir -p $(dirname ${OUT_FILE})
# 
# 	CONFIG_LINE=$(awk -v d=${DESIGN_TYPE} -v g=${ASSEMBLY} -v m=${MARK} '{FS="\t"; OFS="\t"}{if($1 == d && $2 == g && $3 == m) print NR }' ${CONFIGS})
# 	cp ${TEMPLATE} ${OUT_FILE}
# 
# 	for i in $(seq 1 $(awk '{FS="\t"; OFS="\t"}{if(NR == 1) print NF}' ${CONFIGS})); do 
# 		KEY=$(awk -v i=$i '{FS="\t"; OFS="\t"}{if(NR == 1)print "{"$i"}"}' ${CONFIGS})
# 		VALUE=$(awk -v i=$i -v c=${CONFIG_LINE} '{FS="\t"; OFS="\t"}{if(NR == c) print $i}' ${CONFIGS})
# 		if [[ $(echo $VALUE |grep -c " ") > 0 ]] ; then 
# 			sed -i "s|${KEY}|\"${VALUE}\"|g" ${OUT_FILE} 
# 		else 
# 			sed -i "s|${KEY}|${VALUE}|g" ${OUT_FILE} 
# 		fi 
# 	done 
#   sed -i 's/"//g' ${OUT_FILE}
# 	#change script dir
# 	if [ ! ${#TARGET_BED[@]} -eq 0 ]
# 	then
#   	for bed in ${TARGET_BED[@]}; do
#   		if [[ -f $bed ]]
#     	then
#     	  sed -i "s|BED_FEATURES = .*|BED_FEATURES = ${bed}|g"  ${OUT_FILE}
#       else 
#         echo "!! !! !! !!! !!"
#         echo "WARNING, file $bed doesn't exist !"
#         echo "!! !! !! !!! !!"
#         
#       fi
#     done
#   fi
# 	echo -e "COMMENT: Configuration file is ${OUT_FILE} \n"
# fi 
# }

getConf_func(){
	if [[ $# != "6" && $# != "7" ]]; then 
		echo "All arguments must be filled : getConf_func TEMPLATE CONFIGS DESIGN_TYPE ASSEMBLY OUT_FILE MARK [TARGET_BED] "; 
		echo; 
		help_func GetConf;
	else 	
		local TEMPLATE=${1}
		local CONFIGS=${2}
		local DESIGN_TYPE=${3}
		local ASSEMBLY=${4}
		local TYPE=${5}
		local OUT_FILE=${6}
		local TARGET_BED=($(echo ${7} | sed 's/,/ /g'))

		mkdir -p $(dirname ${OUT_FILE})

		CONFIG_LINE=$(awk -v d=${DESIGN_TYPE} -v g=${ASSEMBLY} -v m=${MARK} '{FS="\t"; OFS="\t"}{if($1 == d && $2 == g && $3 == m) print NR }' ${CONFIGS})
		echo "CONFIG_LINE: ${CONFIG_LINE}"
		cp ${TEMPLATE} ${OUT_FILE}

		for i in $(seq 1 $(awk '{FS="\t"; OFS="\t"}{if(NR == 1) print NF}' ${CONFIGS})); do 
			KEY=$(awk -v i=$i '{FS="\t"; OFS="\t"}{if(NR == 1)print "{"$i"}"}' ${CONFIGS})
			VALUE=$(awk -v i=$i -v c=${CONFIG_LINE} '{FS="\t"; OFS="\t"}{if(NR == c) print $i}' ${CONFIGS})
			echo "KEY: ${KEY}, VALUE: ${VALUE}"
			if [[ $(echo $VALUE |grep -c " ") > 0 ]] ; then 
				sed -i "s|${KEY}|\"${VALUE}\"|g" ${OUT_FILE} 
			else 
				sed -i "s|${KEY}|${VALUE}|g" ${OUT_FILE} 
			fi 
		done 
		sed -i 's/"//g' ${OUT_FILE}
		echo "OUT_FILE contents: "
		cat ${OUT_FILE}
		#change script dir
		if [ ! ${#TARGET_BED[@]} -eq 0 ]
		then
			for bed in ${TARGET_BED[@]}; do
				if [[ -f $bed ]]
				then
					sed -i "s|BED_FEATURES = .*|BED_FEATURES = ${bed}|g"  ${OUT_FILE}
				else 
					echo "!! !! !! !!! !!"
					echo "WARNING, file $bed doesn't exist !"
					echo "!! !! !! !!! !!"
				fi
			done
		fi
		echo -e "COMMENT: Configuration file is ${OUT_FILE} \n"
	fi 
}



## Create FASTQs from BCLs and reverse the index
## 
## $1 = BCL directory
## $2 = Output directory
## $3 = SampleSheet (SampleSheet.csv obtained from the KDI)
## $4 = NGS Sample Name
## $5 Location of the reads and indexes: Y101,I8,Y69,Y51 --> First 101bp of genomic DNA (first part of read),
##  then 8bp of AR adapter, then 69bp of barcode DNA, then 51bp of genomic DNA (second part of read)
## $6 = log directory
reverse_fastq_func()
{

    output_dir=$1
    index=$2
    final_name=$3


    local log=$4/reverse_fastq_func.log
    echo -e "Running Reverse FASTQ ..."
    echo -e "Logs: $log"
    echo

    rm -rf ${output_dir}/fastqs/

    mkdir -p ${output_dir}/fastqs/

    echo "Running reversing FASTQ ..." > ${log}
    echo >> ${log}

    # Reverse Index containing cell barcodes for correct mapping
    cmd="/bioinfo/local/build/fastx_toolkit_0.0.13/fastx_reverse_complement -Q33 -i <(gzip -cd $index) -z -o $output_dir/fastqs/${final_name}.R2.fastq.gz"
    exec_cmd ${cmd} >> ${log} 2>&1

   echo
   echo "Finished reversing fastqs !"

}

## Trim reads
## -m <mode>
## $1 = fastq
## $2 = length from 5' to keep
## $3 = out dire
## $4 = log dir
fastx_trimmer_func()
{
    ##Logs                                                                                                                                                                                                                               
    local log=$4/trim.log
    echo -e "Running Fastx trimmer ..."
    echo -e "Logs: $log"
    echo
    
    ##Mode
    mode=$5
    prefix=$6
    
    ## -Q 33
    ## Centos version ?
    local out=$3
    mkdir -p ${out}
    out_prefix=${out}/${prefix}
	
	if [ $mode == 'barcode' ] 
		then
		local ofile=${out_prefix}_trimmed_BC.R2.fastq
    		if [[ $1 =~ \.gz ]]; then
           			cmd="${FASTX_PATH}fastx_trimmer -Q 33 -l $2 -i <(gzip -cd $1) -o ${ofile}"
      	 		else
           			cmd="${FASTX_PATH}fastx_trimmer -Q 33 -l $2 -i $1 -o ${ofile}"
       		fi
    		exec_cmd ${cmd} > ${log} 2>&1	
		else			#Change option -l (last) to -f (first) to remove the barcode and linker from R2 reads in preparation for genome alignment
			local ofile=${out_prefix}_trimmed_G.R2.fastq
                	if [[ $1 =~ \.gz ]]; then
                        	cmd="${FASTX_PATH}fastx_trimmer -Q 33 -f $2 -i <(gzip -cd $1) -o ${ofile}"
                        else
                        	cmd="${FASTX_PATH}fastx_trimmer -Q 33 -f $2 -i $1 -o ${ofile}"
                	fi		
    		exec_cmd ${cmd} > ${log} 2>&1	
		fi

    cmd="gzip ${ofile}"
    exec_cmd ${cmd} >> ${log} 2>&1
}

#BWA MEM mapping function 
bwa_func()
{
   ## Logs
   mkdir -p $3
   local log=$3/mapping_bwa.log
   echo -e "Running BWA mapping ..."
   echo -e "Logs: $log"
   echo

   odir=$2  
   prefix=$4
   out_prefix=${odir}/${prefix}
   ## input type [SE/PE]
   inputs=($1)
   if [[ ${#inputs[@]} -eq 1 ]]; then
       if [[ ${inputs[0]} =~ \.gz ]]; then
           cmd_in=" <(gzip -cd ${inputs[0]})"
       else
           cmd_in="${inputs[0]}"
       fi
   elif [[ ${#inputs[@]} -eq 2 ]]; then
       if [[ ${inputs[0]} =~ \.gz ]]; then
           cmd_in="<(gzip -cd ${inputs[0]}) <(gzip -cd ${inputs[1]})"
       else
           cmd_in="${inputs[0]} ${inputs[1]}"
       fi
   fi

   ## Run Mapping
   local out=$2
   mkdir -p ${out}

   cmd="bwa mem -t ${NB_PROC} -M ${GENOME_MAPPING_OPTS_BWA} ${GENOME_IDX_PATH_BWA}  ${cmd_in} > ${out_prefix}.sam"
   exec_cmd ${cmd} > ${log} 2>&1
   
   #Remove multimappers
   cmd="samtools view -H ${out_prefix}.sam | sed '/^@CO/ d' > ${out_prefix}_header.sam"
   exec_cmd ${cmd} >> ${log} 2>&1
   #Remove secondary alignments that perturb the line order 
   cmd="samtools view -F 256 ${out_prefix}.sam | awk -v OFS='\t' 'NR%2==1{r1=\$0;mapped=\$3; if( \$0 ~ /XA:Z:|SA:Z:/){p=\"no\"}else{p=\"yes\"}} NR%2==0{if( \$0 ~ /XA:Z:|SA:Z:/){} else {if(p==\"yes\" && mapped != \"*\"){print r1\"\n\"\$0}} }' > ${out_prefix}_noMultimappers.sam"
   exec_cmd ${cmd} >> ${log} 2>&1
   
   cmd="cat ${out_prefix}_header.sam ${out_prefix}_noMultimappers.sam > ${out_prefix}.sam"
   exec_cmd ${cmd} >> ${log} 2>&1

   #Reconvert to BAM
   cmd="samtools view -@ ${NB_PROC} -bS ${out_prefix}.sam > ${out_prefix}.bam"
   exec_cmd ${cmd} >> ${log} 2>&1

	 #Sort
   cmd="samtools sort -n -@ ${NB_PROC} ${out_prefix}.bam -o ${out_prefix}_nsorted.bam && mv ${out_prefix}_nsorted.bam ${out_prefix}.bam"
   exec_cmd ${cmd} >> ${log} 2>&1

   ## Clean sam file
  # cmd="rm ${out_prefix}.sam"
  # exec_cmd ${cmd} >> ${log} 2>&1
}

#BOWTIE mapping function 
bowtie_func()
{
   ## Logs
   mkdir -p $3
   local log=$3/mapping_bowtie.log
   echo -e "Running Bowtie mapping ..."
   echo -e "Logs: $log"
   echo
   star_mapping=$5
   
   ## input type 
   inputs=($1)
   if [[ ${#inputs[@]} -eq 1 ]]; then
       if [[ ${inputs[0]} =~ \.gz ]]; then
           cmd_in=" <(gzip -cd ${inputs[0]})"
       else
           cmd_in="${inputs[0]}"
       fi
   elif [[ ${#inputs[@]} -eq 2 ]]; then
       if [[ ${inputs[0]} =~ \.gz ]]; then
           cmd_in="-1 <(gzip -cd ${inputs[0]}) -2 <(gzip -cd ${inputs[1]})"
       else
           cmd_in="-1 ${inputs[0]} -2 ${inputs[1]}"
       fi
   fi

   ## sample_id
   if [ ! -z ${SAMPLE_ID} ]; then
       cmd_id="--sam-RG ID:${SAMPLE_ID} --sam-RG SM:${SAMPLE_ID} --sam-RG LB:${SAMPLE_ID} --sam-RG PU:${SAMPLE_ID} --sam-RG PL:ILLUMINA"
   else
       cmd_id=""
   fi

   ## Run
   local out=$2
   mkdir -p ${out}
   local out_prefix=${out}/$4

   cmd="bowtie -p ${NB_PROC} -S -m1 --sam ${MAPPING_IDX_BOWTIE1} ${cmd_in} > ${out_prefix}.sam"
   exec_cmd ${cmd} > $log 2>&1
   
  
   cmd="grep -v '^@' ${out_prefix}.sam  | awk -v OFS=\"\t\" '\$3!=\"*\"{print \$1,\$3,\$4}' | sort -T ${TMP_DIR} --parallel=${NB_PROC} -k1,1 >  ${out_prefix}.nsorted.txt"
   exec_cmd ${cmd} >> ${log} 2>&1
   
   cmd="rm -f ${out_prefix}.sam"
   exec_cmd ${cmd} >> ${log} 2>&1

   cmd="wc -l  ${out_prefix}.nsorted.txt >  ${out_prefix}.bowtie_uniquely_mapped"
   exec_cmd ${cmd} >> ${log} 2>&1
   
   cmd="samtools view -H ${star_mapping}  > ${out_prefix}.header.sam"
   exec_cmd ${cmd} >> ${log} 2>&1
   
   cmd="samtools view -F 256 ${star_mapping}  | sort -T ${TMP_DIR} --parallel=${NB_PROC} -k1,1  > ${out_prefix}.star.nsorted.sam"
   exec_cmd ${cmd} >> ${log} 2>&1
  
   cmd="rm -f  ${star_mapping} "
   exec_cmd ${cmd} >> ${log} 2>&1

   ##Join STAR & Bowtie -m1 bam files
   cmd="join -1 1  -2 1 ${out_prefix}.nsorted.txt ${out_prefix}.star.nsorted.sam  > ${out_prefix}.joined.txt"
   exec_cmd ${cmd} >> ${log} 2>&1

   cmd="rm -f  ${out_prefix}.nsorted.txt ${out_prefix}.star.nsorted.sam"
   exec_cmd ${cmd} >> ${log} 2>&1

   cmd="awk -v OFS=\"\t\" -v FS=' ' 'NR%2==1{a=\$1\"\t\"\$4\"\t\"\$5\"\t\"\$6\"\t\"\$7\"\t\"\$8\"\t\"\$9\"\t\"\$10\"\t\"\$11\"\t\"\$12\"\t\"\$13\"\t\"\$14\"\t\"\$15\"\t\"\$16\"\t\"\$17\"\t\"\$18\"\t\"\$19\"\t\"\$20; next} {if(\$2==\$5 && (\$3>=\$6-200 && \$3 <= \$6+200) ){print a\"\n\"\$1,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13,\$14,\$15,\$16,\$17,\$18,\$19,\$20}}'  ${out_prefix}.joined.txt | sed -e 's/\s*$//g' > ${out_prefix}.joined.sam"
   exec_cmd ${cmd} >> ${log} 2>&1
   
   cmd="cat ${out_prefix}.header.sam ${out_prefix}.joined.sam > ${out_prefix}.joined.2.sam"
   exec_cmd ${cmd} >> ${log} 2>&1
   
   cmd="rm -f  ${out_prefix}.joined.sam"
   exec_cmd ${cmd} >> ${log} 2>&1

   cmd="samtools sort -n -@ ${NB_PROC} ${out_prefix}.joined.2.sam -o ${out_prefix}_nsorted.bam && mv ${out_prefix}_nsorted.bam ${out_prefix}.bam"
   exec_cmd ${cmd} >> ${log} 2>&1
   
   cmd="rm -f  ${out_prefix}.joined.2.sam"
   exec_cmd ${cmd} >> ${log} 2>&1
}

## STAR MAPPING FUNCTION
star_func()
{
   ## Logs
   mkdir -p $3
   local log=$3/mapping_star.log
   echo -e "Running STAR mapping ..."
   echo -e "Logs: $log"
   echo

   odir=$2  
   prefix=$4
   out_prefix=${odir}/${prefix}
   ## input type [SE/PE]
   inputs=($1)
   if [[ ${#inputs[@]} -eq 1 ]]; then
       if [[ ${inputs[0]} =~ \.gz ]]; then
           cmd_in=" <(gzip -cd ${inputs[0]})"
       else
           cmd_in="${inputs[0]}"
       fi
  elif [[ ${#inputs[@]} -eq 2 ]]; then
       if [[ ${inputs[0]} =~ \.gz ]]; then
           cmd_in="<(gzip -cd ${inputs[0]})"
       else
           cmd_in="${inputs[0]}"
       fi
       if [[ ${inputs[1]} =~ \.gz ]]; then
           cmd_in="$cmd_in <(gzip -cd ${inputs[1]})"
       else
           cmd_in="$cmd_in ${inputs[1]}"
       fi
   fi

   ## Run Mapping
   local out=$2
   mkdir -p ${out}
   
     cmd="STAR --runMode alignReads --runThreadN ${NB_PROC} ${MAPPING_OPTS_STAR} --readFilesIn ${cmd_in} --genomeDir ${MAPPING_INDEX_STAR} --outFileNamePrefix ${out}/"
   exec_cmd ${cmd} > ${log} 2>&1
   
   #Reconvert to BAM
   cmd="samtools view -@ ${NB_PROC} -bS ${out}/Aligned.out.sam > ${out_prefix}.bam"
   exec_cmd ${cmd} >> ${log} 2>&1

	 #Sort
   cmd="samtools sort -n -@ ${NB_PROC} ${out_prefix}.bam -o ${out_prefix}_nsorted.bam && mv ${out_prefix}_nsorted.bam ${out_prefix}.bam"
   exec_cmd ${cmd} >> ${log} 2>&1

   ## Clean sam file
   cmd="rm ${out}/Aligned.out.sam"
   exec_cmd ${cmd} >> ${log} 2>&1
}

barcode_index_mapping_func()
{
##Function mappin barcodes from single-cell ChIP-seq data, index by index (3 different indexes) to the reference (forward) indexes 
#Input : Read 2 (.fastq.gz) , 		out,  		${PREFIX} ,log_dir
#Input :     $1             ,              $2      ,       $3        , 	$4 
read2=$1
out=$2
prefix=$3
## logs
    mkdir -p $4 
    mkdir -p $out
    local log=$4/index_mapping_bowtie2.log
    echo -e "Running Bowtie2 index (barcode) mapping ..."
    echo -e "Logs: $log"
    echo
  
  #Beads type [Hifibio | LBC] -> different Index lengths [20 + 4 | 16 + 4]
  if [[ ${BARCODE_LENGTH} -eq 60 ]]
  	then
  	echo -e "Short barcode :LBC"
      start_index_1=6
  	  start_index_2=26
  	  start_index_3=46
  #     start_index_1=5
  # 	  start_index_2=25
  # 	  start_index_3=45
  	  size_index=16
         else
  	echo -e "Long barcode : Hifibio "
      start_index_1=1
    	start_index_2=25
    	start_index_3=49
    	size_index=20
  fi

##Extract three indexes from reads : 1 - 16 = index 1 ; 21 - 36 = index 2; 41 - 56 = index 3
if [[ ${read2} =~ \.gz ]]; then
         cmd="gzip -cd  $read2 | awk -v start_index_1=$start_index_1 -v size_index=$size_index  'NR%4==1{print \">\"substr(\$0,2)}; NR%4==2{print substr(\$0,start_index_1,size_index)}' > ${out}/read_indexes_1.fasta"
     else
         cmd="cat $read2 | awk -v start_index_1=$start_index_1 -v size_index=$size_index  'NR%4==1{print \">\"substr(\$0,2)}; NR%4==2{print substr(\$0,start_index_1,size_index)}' > ${out}/read_indexes_1.fasta"
 fi
  exec_cmd ${cmd} > ${log} 2>&1
  
 if [[ ${read2} =~ \.gz ]]; then
         cmd="gzip -cd $read2 | awk -v start_index_2=$start_index_2 -v size_index=$size_index 'NR%4==1{print  \">\"substr(\$0,2)}; NR%4==2{print substr(\$0,start_index_2,size_index)}' > ${out}/read_indexes_2.fasta"
     else
         cmd="cat $read2 | awk -v start_index_2=$start_index_2 -v size_index=$size_index 'NR%4==1{print  \">\"substr(\$0,2)}; NR%4==2{print substr(\$0,start_index_2,size_index)}' > ${out}/read_indexes_2.fasta"
 fi
  exec_cmd ${cmd} >> ${log} 2>&1

 if [[ ${read2} =~ \.gz ]]; then
       cmd="gzip -cd  $read2 | awk -v start_index_3=$start_index_3 -v size_index=$size_index 'NR%4==1{print \">\"substr(\$0,2)}; NR%4==2{print substr(\$0,start_index_3,size_index)}' > ${out}/read_indexes_3.fasta"
    else
 cmd="cat $read2 | awk -v start_index_3=$start_index_3 -v size_index=$size_index 'NR%4==1{print \">\"substr(\$0,2)}; NR%4==2{print substr(\$0,start_index_3,size_index)}' > ${out}/read_indexes_3.fasta"    
fi  
  exec_cmd ${cmd} >> ${log} 2>&1
  
  #Map INDEXES 1 against Index1 library
  cmd="bowtie2 -x ${BARCODE_BOWTIE_IDX_PATH}1 -f ${out}/read_indexes_1.fasta ${BARCODE_MAPPING_OPTS} -p ${NB_PROC} > ${out}/index_1_bowtie2.sam"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  #Keep only reads that were matched by a unique index 1 + counting matched index1
  cmd="awk -v out=\"${out}/\" '/XS/{next} \$2!=4{print \$1,\$3;count++} ;END{print count > out\"/count_index_1\"}' ${out}/index_1_bowtie2.sam > ${out}/reads_matching_index_1.txt"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  #Map INDEXES 2 against Index2 library
  cmd="bowtie2 -x ${BARCODE_BOWTIE_IDX_PATH}2 -f ${out}/read_indexes_2.fasta ${BARCODE_MAPPING_OPTS} -p ${NB_PROC} > ${out}/index_2_bowtie2.sam"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  #Keep only reads that were matched by a unique index 2 + counting matched index2
  cmd="awk -v out=\"${out}/\" '/XS/{next} \$2!=4{print \$1,\$3;count++} ;END{print count > out\"count_index_2\"}' ${out}/index_2_bowtie2.sam > ${out}/reads_matching_index_2.txt"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  #Map INDEXES 3 against Index3 library
  cmd="bowtie2 -x ${BARCODE_BOWTIE_IDX_PATH}3 -f ${out}/read_indexes_3.fasta ${BARCODE_MAPPING_OPTS} -p ${NB_PROC} > ${out}/index_3_bowtie2.sam"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  #Keep only reads that were matched by a unique index 3 + counting matched index3
  cmd="awk -v out=\"${out}/\" '/XS/{next} \$2!=4{print \$1,\$3;count++} ;END{print count > out\"count_index_3\"}' ${out}/index_3_bowtie2.sam > ${out}/reads_matching_index_3.txt"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  ##Sort indexes by read name: 
  cmd="sort -T ${TMP_DIR} --parallel=${NB_PROC} -k1,1 ${out}/reads_matching_index_1.txt > ${out}/reads_matching_index_1_sorted.txt"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  cmd="rm ${out}/reads_matching_index_1.txt"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  cmd="sort -T ${TMP_DIR} --parallel=${NB_PROC} -k1,1 ${out}/reads_matching_index_2.txt > ${out}/reads_matching_index_2_sorted.txt"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  cmd="rm ${out}/reads_matching_index_2.txt"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  cmd="sort -T ${TMP_DIR} --parallel=${NB_PROC} -k1,1 ${out}/reads_matching_index_3.txt > ${out}/reads_matching_index_3_sorted.txt"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  cmd="rm ${out}/reads_matching_index_3.txt"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  #Join indexes 1 & 2 together (inner join)
  cmd="join -t$' ' -1 1 -2 1 ${out}/reads_matching_index_1_sorted.txt ${out}/reads_matching_index_2_sorted.txt > ${out}/tmp"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  #Count matched index 1 & 2
  cmd="echo \$(wc -l ${out}/tmp) | cut -d' ' -f1 > ${out}/count_index_1_2"
  exec_cmd ${cmd} >> ${log} 2>&1	
  
  #Join indexes (1 & 2) & 3 together to recompose full barcode (inner join)
  cmd="join -t$' ' -1 1 -2 1 ${out}/tmp ${out}/reads_matching_index_3_sorted.txt > ${out}/final"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  #Reformat & count matched index (1 & 2 & 3) <=> barcode
  cmd="awk -v out=\"${out}/\" '{print substr(\$1,1)\"\tBC\"substr(\$2,2)substr(\$3,2)substr(\$4,2);count++} ;END{print count > out\"count_index_1_2_3\"}' ${out}/final > ${out}/${prefix}_read_barcodes.txt"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  ##Write logs
  cmd="n_index_1=\$(cat ${out}/count_index_1)"
  exec_cmd ${cmd} >> ${log} 2>&1
  cmd="n_index_2=\$(cat ${out}/count_index_2)"
  exec_cmd ${cmd} >> ${log} 2>&1
  cmd="n_index_3=\$(cat ${out}/count_index_3)"
  exec_cmd ${cmd} >> ${log} 2>&1
  cmd="n_index_1_2=\$(cat ${out}/count_index_1_2)"
  exec_cmd ${cmd} >> ${log} 2>&1
  cmd="n_index_1_2_3=\$(cat ${out}/count_index_1_2_3)"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  echo "## Number of matched indexes 1: $n_index_1" >> ${log}
  echo "## Number of matched indexes 2: $n_index_2" >> ${log}
  echo "## Number of matched indexes 1 and 2: $n_index_1_2" >> ${log}
  echo "## Number of matched indexes 3: $n_index_3" >> ${log}
  echo "## Number of matched barcodes: $n_index_1_2_3" >> ${log}
  
  ## Remove all non used files
  cmd="rm -f ${out}/*.sam ${out}/read* ${out}/tmp ${out}/final ${out}/count"
  exec_cmd ${cmd} >> ${log} 2>&1
}

## bowtie2
bowtie2_func()
{

    ## logs
    local log=$4/mapping_bowtie2.log
    echo -e "Running Bowtie2 mapping ..."
    echo -e "Logs: $log"
    echo

    ## input type
    inputs=($1)
    if [[ ${#inputs[@]} -eq 1 ]]; then
        cmd_in="-U $1"
    elif [[ ${#inputs[@]} -eq 2 ]]; then
        cmd_in="-1 ${inputs[0]} -2 ${inputs[1]}"
    fi

    ## sample_id
    if [ ! -z ${SAMPLE_ID} ]; then
        cmd_id="--rg-id ${SAMPLE_ID} --rg PL:ILLUMINA"
    else
        cmd_id=""
    fi

    ## Run Mapping
    local out=$2
    local prefix=$3
    mkdir -p ${out}
    local out_prefix=${out}/${prefix}

    cmd="bowtie2 -t -p ${NB_PROC} -x ${MAPPING_INDEX} ${MAPPING_OPTS} ${cmd_id} ${cmd_in} > ${out_prefix}.bam"
    exec_cmd ${cmd} > ${log} 2>&1

    cmd="samtools sort -n -@ ${NB_PROC} ${out_prefix}.bam -o ${out_prefix}_sorted.bam && mv ${out_prefix}_sorted.bam ${out_prefix}.bam"
    exec_cmd ${cmd} >> ${log} 2>&1
}

##BigWig function
bw_func()
{
    
    #Load Deeptools 3.1.0
    source /bioinfo/local/build/Centos/miniconda/miniconda3-4.4.6/bin/activate /bioinfo/local/build/Centos/envs_conda/deeptools_3.1.0

    local out=$2
    mkdir -p ${out}
    local log=$3
    mkdir -p ${log}
    log=$log/bamCoverage.log

    echo -e "Generate bigwig file(s) ..."
    echo -e "Logs: $log"
    echo
 	
    name=$(basename $1 ".bam")
    if [[ ! -z ${ENCODE_BLACKLIST} && -e ${ENCODE_BLACKLIST} ]]; then
        local cmd="bamCoverage --bam $1 --outFileName $out/${name}.bw --numberOfProcessors ${NB_PROC} --normalizeUsing CPM --ignoreForNormalization chrX --binSize 50 --smoothLength 500 --extendReads 150 --blackListFileName ${ENCODE_BLACKLIST}"
    else
        local cmd="bamCoverage --bam $1 --outFileName $out/${name}.bw --numberOfProcessors ${NB_PROC} --normalizeUsing CPM --ignoreForNormalization chrX --binSize 50 --smoothLength 500 --extendReads 150"
    fi
    exec_cmd ${cmd} > ${log} 2>&1

    source /bioinfo/local/build/Centos/miniconda/miniconda3-4.4.6/bin/deactivate
}



## add_cellBarcode_func on Paired End Bam file and transform into Single End Bam
#usage : add_cellBarcode_func ${GENOME_BAM} ${BARCODE_BAM} ${ODIR}/mapping/ ${LOGDIR}
add_cellBarcode_func () {
## logs
  local log=$4/add_cellBarcode_func.log 
  echo -e "Running add_cellBarcode_func..."
  echo -e "Logs: $log"
  echo
  local out=$3
  local out_prefix=${out}/$(basename $1 | sed -e 's/.bam$//')

  #Remove secondary aligned reads (256 <=> "not primary alignment") & If read is unmapped or multimapped (NH != 1), tag Read with flag "4" <=> "unmapped" & "chr" = '*'
  cmd="samtools view -F 256 $1 | awk -v OFS='\t' '{if(\$12==\"NH:i:1\"){print \$0} else{\$2=4;\$3=\"*\";\$4=0;\$6=\"*\";print \$0}}'> ${out_prefix}.sam"
  exec_cmd ${cmd} > ${log} 2>&1

  
  #Remove pair of reads that are both unmapped
  cmd="cat ${out_prefix}.sam | awk -v OFS='\t' 'NR%2==1{if(\$3==\"*\"){mapped=0} else{mapped=1} R1=\$0}; NR%2==0{if(\$3==\"*\" && mapped==0){} else {print R1\"\\n\"\$0}}'> ${out_prefix}_1.sam"
  exec_cmd ${cmd} > ${log} 2>&1

  #If read2 & unmapped -> set R2 position as '2147483647'
  cmd="cat ${out_prefix}_1.sam | awk -v OFS='\t' 'NR%2==1{print \$0} NR%2==0{if(\$3==\"*\"){\$4=2147483647;print \$0} else{print \$0} }' > ${out_prefix}_2.sam"
  exec_cmd ${cmd} >> ${log} 2>&1

  #If read 1 is unmapped R1 -> set R1 position as '2147483646'
  cmd="cat ${out_prefix}_2.sam | awk -v OFS='\t' 'NR%2==0{print \$0} NR%2==1{if(\$3==\"*\"){\$4=2147483646;print \$0} else{print \$0} }' > ${out_prefix}_3.sam"
  exec_cmd ${cmd} >> ${log} 2>&1
 
  #Remove comments from the header that produce bugs in the count phase
  cmd="samtools view -H $1 | sed '/^@CO/ d' > ${out_prefix}_header.sam"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  cmd="cat ${out_prefix}_3.sam >> ${out_prefix}_header.sam && mv ${out_prefix}_header.sam ${out_prefix}.sam && samtools view -b -@ ${NB_PROC} ${out_prefix}.sam > ${out_prefix}_unique.bam"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  cmd="rm -f ${out_prefix}_*.sam ${out_prefix}.sam"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  # Keeping R1 aligned + R2 start as tag 'XS' (Switch from Paired End Bam to Single End Bam)
  cmd="samtools view ${out_prefix}_unique.bam | awk '{OFS = \"\t\" ; if(NR%2==1 && !(\$3==\"*\")) {R1=\$0; R1_mapped=1} else if(NR%2==1){R1=\$10; R1_mapped=0; tagR1Pos=\"XS:i:\"\$4}; if(NR%2==0 && !(R1_mapped==0)){tagR2Seq=\"XD:Z:\"\$10; tagR2Pos=\"XS:i:\"\$4;print R1,tagR2Pos,tagR2Seq} if(NR%2==0 && R1_mapped==0){print \$0,tagR1Pos,\"XD:Z:\"R1 } }' > ${out_prefix}_unique.sam"
  exec_cmd ${cmd} >> ${log} 2>&1
   
  #Sort and join on read names reads barcoded and reads mapped to genome (barcode as tag 'XB') --> filter out unbarcoded OR unmapped reads
  cmd="sort -T ${TMP_DIR} --parallel=${NB_PROC} -k1,1 ${out_prefix}_unique.sam > ${out_prefix}_unique_sorted.sam"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  cmd="join -1 1  -2 1  ${out_prefix}_unique_sorted.sam <(awk -v OFS=\"\t\" '{print \$1,\"XB:Z:\"\$2}' $2) > ${out_prefix}_flagged.sam"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  cmd="sed -i 's/ /\t/g' ${out_prefix}_flagged.sam"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  #Remove comments from the header that produce bugs in the count phase
  cmd="samtools view -H ${out_prefix}_unique.bam | sed '/^@CO/ d' > ${out_prefix}_header.sam"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  cmd="cat ${out_prefix}_flagged.sam >> ${out_prefix}_header.sam && mv ${out_prefix}_header.sam ${out_prefix}_flagged.sam && samtools view -@ ${NB_PROC} -b ${out_prefix}_flagged.sam > ${out_prefix}_flagged.bam"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  #Cleaning
  cmd="rm -f ${out_prefix}_unique.bam ${out_prefix}_flagged.sam ${out_prefix}_unique_sorted.sam"
  exec_cmd ${cmd} >> ${log} 2>&1
}

##remove_PCR_RT_duplicates_func (Sort and remove consecutive reads if identified as #PCR or #RT)
#usage : remove_PCR_RT_duplicates_func ${GENOME_BAM_FLAGGED} ${ODIR}/mapping ${LOGDIR}
remove_PCR_RT_duplicates_func(){
  ## logs
  local log=$3/remove_PCR_RT_duplicates_func.log 
  echo -e "Running remove_PCR_RT_duplicates_func..."
  echo -e "Logs: $log"
  echo
  local out=$2
  local out_prefix=${out}/$(basename $1 | sed -e 's/_flagged.bam$//')
    
  ##Sort by barcode then chromosome then position R2
  #Find the column containing the barcode tag XB
  cmd="barcode_field=\$(samtools view ${out_prefix}_flagged.bam  | sed -n \"1 s/XB.*//p\" | sed 's/[^\t]//g' | wc -c)"
  exec_cmd ${cmd} > ${log} 2>&1
  
  echo "The barcode field is $barcode_field" >> ${log}
  
  cmd="printf '@HD\tVN:1.4\tSO:unsorted\n' > ${out_prefix}_header.sam"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  cmd="samtools view -H ${out_prefix}_flagged.bam | sed '/^@HD/ d' >> ${out_prefix}_header.sam"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  #Sort by barcode then chromosome then position R1 (for the PCR removal) 
  #It is important to sort by R1 pos also because the removal is done by comparing consecutive lines ! 
  cmd="samtools view ${out_prefix}_flagged.bam | LC_ALL=C sort -T ${TMP_DIR} --parallel=${NB_PROC} -t $'\t' -k \"$barcode_field.6\" -k 3.4,3g -k 4,4n >> ${out_prefix}_header.sam && samtools view -@ ${NB_PROC} -b ${out_prefix}_header.sam > ${out_prefix}_flagged.sorted.bam"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  #Create count Table from flagged (already sorted by barcode)
  cmd="samtools view ${out_prefix}_flagged.sorted.bam | awk -v bc_field=$barcode_field '{print substr(\$bc_field,6)}' |  uniq -c > ${out_prefix}_flagged.count"
  exec_cmd ${cmd} >> ${log} 2>&1
 
  echo -e "Removing PCR duplicates using R1 POSITION"
  #Remove PCR duplicates = read having same barcode, R1 position, same R2 POSITION, same chr ("exactly equal")
  cmd="samtools view ${out_prefix}_flagged.sorted.bam | awk -v bc_field=$barcode_field -v out=${out}/ 'BEGIN{countPCR=0};NR==1{print \$0;lastChrom=\$3;lastBarcode=\$bc_field; lastR1Pos=\$4} ; NR>=2{R1Pos=\$4; if( (R1Pos==lastR1Pos) && ( \$3==lastChrom ) && (\$bc_field==lastBarcode) ){countPCR++;next} {print \$0;lastR1Pos=\$4;lastChrom=\$3;lastBarcode=\$bc_field }} END {print countPCR > out\"count_PCR_duplicates\"}' > ${out_prefix}_flagged_rmPCR.sam"
  exec_cmd ${cmd} >> $log 2>&1
 
  cmd="samtools view -H ${out_prefix}_flagged.sorted.bam  | sed '/^@CO/ d' > ${out_prefix}_header.sam"
  exec_cmd ${cmd} >> ${log} 2>&1
  cmd="cat ${out_prefix}_flagged_rmPCR.sam >> ${out_prefix}_header.sam && samtools view -@ ${NB_PROC} -b ${out_prefix}_header.sam > ${out_prefix}_flagged_rmPCR.bam"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  #Create count Table from flagged - PCR dups (already sorted by barcode)
  cmd="samtools view ${out_prefix}_flagged_rmPCR.bam | awk -v bc_field=$barcode_field '{print substr(\$bc_field,6)}' |  uniq -c > ${out_prefix}_flagged_rmPCR.count"
  exec_cmd ${cmd} >> $log 2>&1
  
  ## Sort flagged_rmPCR file
  cmd="samtools sort -@ ${NB_PROC} ${out_prefix}_flagged_rmPCR.bam > ${out_prefix}_flagged_rmPCR_sorted.bam"
  exec_cmd ${cmd} >> $log 2>&1
  
  ## Rename flagged_rmPCR file
  cmd="mv ${out_prefix}_flagged_rmPCR_sorted.bam ${out_prefix}_flagged_rmPCR.bam"
  exec_cmd ${cmd} >> $log 2>&1
  
  ## Index flagged_rmPCR file
  cmd="samtools index ${out_prefix}_flagged_rmPCR.bam"
  exec_cmd ${cmd} >> $log 2>&1
  
  echo -e "Not removing RT duplicates"
  ## Copy flagged_rmPCR to flagged_rmPCR_RT
  cmd="cp ${out_prefix}_flagged_rmPCR.bam ${out_prefix}_flagged_rmPCR_RT.bam"
  exec_cmd ${cmd} >> $log 2>&1
  cmd="cp ${out_prefix}_flagged_rmPCR.count ${out_prefix}_flagged_rmPCR_RT.count"
  exec_cmd ${cmd} >> $log 2>&1
  ## Set RT duplicate count to 0
  cmd="echo 0 > ${out}/count_RT_duplicates"
  exec_cmd ${cmd} >> $log 2>&1
   
  ## Index flagged_rmPCR_RT file
  cmd="samtools index ${out_prefix}_flagged_rmPCR_RT.bam"
  exec_cmd ${cmd} >> $log 2>&1
  
  #Write logs
  cmd="n_mapped_barcoded=\$(samtools view -c  ${out_prefix}_flagged.bam)"
  exec_cmd ${cmd} >> $log 2>&1
  cmd="n_pcr_duplicates=\$(cat ${out}/count_PCR_duplicates)"
  exec_cmd ${cmd} >> $log 2>&1
  cmd="n_rt_duplicates=\$(cat ${out}/count_RT_duplicates)"
  exec_cmd ${cmd} >> $log 2>&1
  
  ## Rename flagged.sorted -> flagged
  cmd="mv ${out_prefix}_flagged.sorted.bam ${out_prefix}_flagged.bam"
  exec_cmd ${cmd} >> $log 2>&1
 
  cmd="n_unique_except_R1_unmapped_R2=\$(($n_mapped_barcoded - $n_pcr_duplicates - $n_rt_duplicates))"
  exec_cmd ${cmd} >> $log 2>&1
  
  echo "## Number of reads mapped and barcoded: $n_mapped_barcoded" >> ${log}
  echo "## Number of pcr duplicates: $n_pcr_duplicates" >> ${log}
  echo "## Number of rt duplicates: $n_rt_duplicates" >> ${log}
  echo "## Number of R1 mapped but R2 unmapped: 0" >> ${log}
  echo "## Number of reads after PCR and RT removal (not R1 unmapped R2): $n_unique_except_R1_unmapped_R2" >> ${log}

  ## Remove all non used files
  cmd="rm -f ${out}/count* ${out}/*.sam"
  exec_cmd ${cmd} >> $log 2>&1
}

## Remove Duplicates by "Window"
remove_duplicates()
{
    ## logs
    local log=$3/rmDup.log
    echo -e "Removing window duplicates ..."
    echo -e "Logs: $log"
    echo

    local odir=$2
    mkdir -p ${odir}
    local prefix=${odir}/$(basename $1 | sed -e 's/.bam$//')

    if [ ! -z ${DUPLICATES_WINDOW} ]; then
	cmd="${PYTHON_PATH}/python ${SCRIPTS_PATH}/rmDup.py -i ${prefix}.bam -o ${prefix}_rmDup.bam -d ${DUPLICATES_WINDOW} -v "
    else
	cmd="${PYTHON_PATH}/python ${SCRIPTS_PATH}/rmDup.py -i ${prefix}.bam -o ${prefix}_rmDup.bam -v "
    fi
    exec_cmd ${cmd} > ${log} 2>&1
    #Create count Table from flagged - PCR dups - RT dups and window-based rmDup (need to sort by b arcode)
    cmd="barcode_field=\$(samtools view ${prefix}_rmDup.bam  | sed -n \"1 s/XB.*//p\" | sed 's/[^\t]//g' | wc -c)"
    exec_cmd ${cmd} >> $log 2>&1
    cmd="samtools view ${prefix}_rmDup.bam | awk -v bc_field=$barcode_field '{print substr(\$bc_field,6)}' | sort | uniq -c > ${prefix}_rmDup.count"		
    exec_cmd ${cmd} >> $log 2>&1

    ## Index BAM file
    cmd="samtools index ${prefix}_rmDup.bam"
    exec_cmd ${cmd} >> $log 2>&1
    
}

## filter_black_regions 
#Filter out dark regions in bam file
#bam_to_bedGraph ${FLAGGED_RM_DUP_BAM} $ODIR ${LOGS}
filter_black_regions() {
  bam_in=$1
  local odir=$2
  log=local log=$3/filter_black_regions.log

  local prefix=${odir}/$(basename $bam_in | sed -e 's/.bam$//')
  
  echo -e "Filtering out dark regions in bam file..."
  echo -e "Logs: $log"
  echo
  
  mkdir -p ${odir}
  
  if [[ ! -z ${BIN_PATH}/${ENCODE_BLACKLIST} && -e ${BIN_PATH}/${ENCODE_BLACKLIST} ]]; then
    cmd="bedtools intersect -v -abam ${bam_in} -b ${BIN_PATH}/${ENCODE_BLACKLIST} > ${bam_in}.2 && mv ${bam_in}.2 ${bam_in}"
    exec_cmd ${cmd} > $log 2>&1
    
    cmd="samtools index ${bam_in}"
    exec_cmd ${cmd} >> $log 2>&1
  fi
  ## Index BAM file

 echo "Filtered blacklisted regions !"
  
}

## bam_to_bedGraph 
#Generates a bedgraph for sushi plots from a BAM file
#bam_to_bedGraph ${FLAGGED_RM_DUP_BAM} ${FLAGGED_RM_DUP_COUNT} ${FLAGGED_RM_DUP_BEDGRAPH} ${LOGS}
bam_to_bedGraph() {
  bam_in=$1
  bc_count=$2
  local odir=$3
  log=local log=$4/BamToBedGraph.log

  local prefix=${odir}/$(basename $bam_in | sed -e 's/.bam$//')
  
  echo -e "Creating bedgraph from Mapped & Dedup BAM ..."
  echo -e "Logs: $log"
  echo
  
  mkdir -p ${odir}
  
  #Get barcode field & read length
  barcode_field=$(samtools view $bam_in  | sed -n "1 s/XB.*//p" | sed 's/[^\t]//g' | wc -c)
  read_length=$(samtools view $bam_in | awk 'NR==1{print length($10)}')
  
  #Create header
  cmd="samtools view -H $bam_in | sed '/^@HD/ d' > ${prefix}_tmp_header.sam"
  exec_cmd ${cmd} > ${log} 2>&1
    
  #Sort by Barcode, Chr, Pos R1 :
  cmd="samtools view $bam_in | LC_ALL=C sort -T ${TMP_DIR} --parallel=${NB_PROC} -t $'\t' -k \"$barcode_field.8,$barcode_field\"n -k 3.4,3g -k 4,4n >> ${prefix}_tmp_header.sam"
  exec_cmd ${cmd} >> ${log} 2>&1
    
  cmd="samtools view -@ 4 -b ${prefix}_tmp_header.sam > ${prefix}_tmp.sorted.bam"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  #Convert to bedgraph: Input must be sorted by barcode, chr, position R1
  samtools view ${prefix}_tmp.sorted.bam | awk -v rl=$read_length -v bc_field=$barcode_field -v OFS=" " -v count=1 '
  NR==1{
    start=$4;
    end=$4+rl;
    lastBC=substr($bc_field,6,15);
    lastPos=$4;lastChr=$3
    
  }
  NR>1{
  if(lastBC==substr($bc_field,6,15)){
    
    if(lastChr==$3 && $4<lastPos+rl)
    {
      count++;
    } else{
      print lastChr,start,end,count,lastBC;
      count=1;
      start=$4;
      
    }
    }
    else{
       print lastChr,start,end,count,lastBC;
       count=1;
       start=$4;
    }
      lastChr=$3;
      lastPos=$4;
      end=$4+rl;
      lastBC=substr($bc_field,6,15);
  
  }
  ' > ${prefix}_tmp.bedgraph
  
  #Add barcode count information
  cmd="join -t$' ' -o 1.1,1.2,1.3,1.4,1.5,2.1 -1 5 -2 2 ${prefix}_tmp.bedgraph <(sed -s 's/ //g' $bc_count | sed -s 's/BC/ BC/g') | sed 's/ /\\t/g' > ${prefix}_tmp.join.bedgraph"
  exec_cmd ${cmd} >> ${log} 2>&1
  #Sort
  cmd="bedtools sort -i ${prefix}_tmp.join.bedgraph > ${prefix}.bedgraph"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  cmd="rm -f ${prefix}_tmp.bedgraph ${prefix}_tmp_header.sam ${prefix}_tmp.sorted.bam ${prefix}_tmp.join.bedgraph"
  exec_cmd ${cmd} >> ${log} 2>&1

}

## bam_to_bedGraph 
#Generates a bedgraph for sushi plots from a BAM file
#bam_to_bedGraph ${FLAGGED_RM_DUP_BAM} ${FLAGGED_RM_DUP_COUNT} ${FLAGGED_RM_DUP_BEDGRAPH} ${LOGS}
bam_to_sc_bed() {
  bam_in=$1
  count=$2
  local odir=$3
  log=local log=$4/BamToScBed.log


  local prefix=${odir}/$(basename $bam_in | sed -e 's/.bam$//')
  
  echo -e "Creating sc bed from Mapped & Dedup BAM ..."
  echo -e "Logs: $log"
  echo
  
  mkdir -p ${odir}
  for i in $(echo $2 | sed 's/,/ /g'); do mkdir -p ${odir}/scBed_$i/ ; done
  
  #Get barcode field & read length
  barcode_field=$(samtools view $bam_in  | sed -n "1 s/XB.*//p" | sed 's/[^\t]//g' | wc -c)

  #Create header
  cmd="samtools view -H $bam_in | sed '/^@HD/ d' > ${prefix}_tmp_header.sam"
  exec_cmd ${cmd} > ${log} 2>&1
    
  #Sort by Barcode, Chr, Pos R1 :
  cmd="samtools view $bam_in | LC_ALL=C sort -T ${TMP_DIR} --parallel=${NB_PROC} -t $'\t' -k \"$barcode_field\" -k 3.4,3g -k 4,4n >> ${prefix}_tmp_header.sam"
  exec_cmd ${cmd} >> ${log} 2>&1
    
  cmd="samtools view -@ ${NB_PROC} -b ${prefix}_tmp_header.sam > ${prefix}_tmp.sorted.bam"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  #Convert to bedgraph: Input must be sorted by barcode, chr, position R1
  samtools view ${prefix}_tmp.sorted.bam | awk -v odir=${odir}/scBed -v bc_field=$barcode_field -v OFS="\t" -v count=$count '
  BEGIN{
    split(count,min_counts,",")
  }
  NR==1{
    lastBC=substr($bc_field,6,10000);
    i=1
    chr[i] = $3
    start[i] = $4
    end[i] = $4 +1
  }
  NR>1{
  if(lastBC==substr($bc_field,6,10000)){
    i = i +1
    chr[i] = $3
    start[i] = $4
    end[i] = $4 +1
    }
    else{
    for(m=1; m<=length(min_counts);m++){
      if(i > min_counts[m]){
        for (x=1; x<=i; x++){
          out = odir"_"min_counts[m]"/"lastBC".bed"
          print chr[x],start[x],end[x] >> out
        }
      }
    }
    i=0
    }
     lastBC=substr($bc_field,6,10000);
}
'

  #Gzip
  files=$(ls $odir/scBed_${count}/*.bed | head -n1)
  if [ -f $files ];then
        cmd="for i in $odir/scBed_${count}/*.bed; do gzip -9 \$i; done"
        exec_cmd ${cmd} >> ${log} 2>&1
  fi

  cmd="rm -f ${prefix}_tmp_header.sam ${prefix}_tmp.sorted.bam"
  exec_cmd ${cmd} >> ${log} 2>&1

}

## Generate Fragment Files count table
# Create Fragment File
bam_to_fragment_file(){

  in_prefix=$1
  odir=$2
  local prefix=$(basename $in_prefix | sed -e 's/.bam$//')
  out_prefix=${odir}/${prefix}

  ##Sort by barcode then chromosome then position R2
  #Find the column containing the barcode tag XB
  barcode_field=$(/bioinfo/local/build/Centos/samtools/samtools-0.1.19/samtools view $in_prefix |sed -n "1 s/XB.*//p" |sed 's/[^\t]//g' | wc -c)

  #Sort by barcode then chromosome then read position
  samtools view ${in_prefix} | grep -E "XB:Z" | awk -v bc_field=$barcode_field -v OFS="\t" '{gsub("XB:Z:","",$bc_field); print $3,$4,$4+100,$bc_field,1}' > ${out_prefix}.fragments.tsv

  #Compress
  /bioinfo/local/build/Centos/samtools/samtools-1.9/bin/bgzip -@ 8 -f -l 9 ${out_prefix}.fragments.tsv

  ## Index flagged_rmPCR_RT file
  /bioinfo/local/build/Centos/samtools/samtools-1.9/bin/tabix -p bed ${out_prefix}.fragments.tsv.gz


 }


## Generate genomic count table
make_counts(){

    ## logs
    local log=$4/make_counts.log
    echo -e "Generating counts table ..."
    echo -e "Logs: $log"
    echo
	
    local odir=$2
    mkdir -p ${odir}
    local prefix=${odir}/$(basename $1 | sed -e 's/.bam$//')
    local mapping_dir=$3

    #Counting unique BCs
    local bc_prefix=$(basename $1 | sed -e 's/.bam$//')
    cmd="barcodes=\$(wc -l ${mapping_dir}/${bc_prefix}.count | awk '{print \$1}')"
    exec_cmd ${cmd} > ${log} 2>&1

    echo "Barcodes found = $barcodes" >> ${log} 
    for bsize in ${BIN_SIZE}
    do
	echo -e "at ${bsize} resolution ..."
	opts="-b ${bsize} "
        if [ ! -z ${MIN_COUNT_PER_BARCODE_AFTER_RMDUP} ]; then
	    opts="${opts} -f ${MIN_COUNT_PER_BARCODE_AFTER_RMDUP} "
	fi

	cmd="${PYTHON_PATH}/python ${SCRIPTS_PATH}/sc2sparsecounts.py -i $1 -o ${prefix} ${opts} -s $barcodes -t 'XB' -v"
        exec_cmd ${cmd} >> ${log} 2>&1
    done

    for bed in ${BED_FEATURES}
    do
        echo -e "at ${bed} resolution ..."
        opts="-B ${bed} "
        if [ ! -z ${MIN_COUNT_PER_BARCODE_AFTER_RMDUP} ]; then
            opts="${opts} -f ${MIN_COUNT_PER_BARCODE_AFTER_RMDUP} "
        fi
	osuff=$(basename ${bed} | sed -e 's/.bed//')
        cmd="${PYTHON_PATH}/python ${SCRIPTS_PATH}/sc2sparsecounts.py -i $1 -o ${prefix}_${osuff} ${opts} -s $barcodes -t 'XB' -v"
        exec_cmd ${cmd} >> ${log} 2>&1

    done 
}

add_info_to_log(){

genomedir=$1
logdir=$2
odir=$3
prefix=$4
bin_path=$5
arguments=$6
bin_name=$7
R_downstream_dir="$(dirname ${R_DOWNSTREAM})"
echo ${R_downstream_dir}

 ## logs
    local log=$2/add_info_to_log.log
    echo -e "add_info_to_log ..."
    echo -e "Logs: $log"
    echo

cmd="cat ${genomedir}/Log.final.out > ${odir}/${prefix}_scChIPseq_logs.txt"
exec_cmd ${cmd} > ${log} 2>&1
cmd="echo \"Barcode mapping\" >> ${odir}/${prefix}_scChIPseq_logs.txt"
exec_cmd ${cmd} >> ${log} 2>&1

n=$(samtools view -c -@ ${NB_PROC} ${odir}/mapping/${prefix}_flagged.bam)

echo "## Number of matched indexes 1: $n"  >> ${odir}/${prefix}_scChIPseq_logs.txt
echo "## Number of matched indexes 2: $n" >> ${odir}/${prefix}_scChIPseq_logs.txt
echo "## Number of matched indexes 1 and 2: $n" >> ${odir}/${prefix}_scChIPseq_logs.txt
echo "## Number of matched indexes 3: $n" >> ${odir}/${prefix}_scChIPseq_logs.txt
echo "## Number of matched barcodes: $n" >> ${odir}/${prefix}_scChIPseq_logs.txt

cmd="echo \"Barcode Flagging And Exact Duplicate Removal\"  >> ${odir}/${prefix}_scChIPseq_logs.txt"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="grep '##' ${logdir}/remove_PCR_RT_duplicates_func.log >> ${odir}/${prefix}_scChIPseq_logs.txt"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="echo \"Duplicate Removal Window (rmDup)\"  >> ${odir}/${prefix}_scChIPseq_logs.txt"
exec_cmd ${cmd} >> ${log} 2>&1

echo '## Number of duplicates: 0' >> ${odir}/${prefix}_scChIPseq_logs.txt
echo "## Number of reads after duplicates removal: $(samtools view -c ${odir}/mapping/${prefix}_flagged_rmPCR_RT.bam)" >> ${odir}/${prefix}_scChIPseq_logs.txt

make_metadata ${odir} ${prefix} ${odir}/${prefix}_scChIPseq_logs.txt "$odir/mapping/${prefix}_flagged_rmPCR_RT_rmDup.count" 

cmd="mkdir -p ${odir}/.data_engineering_run/"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="cp ${bin_path}/scripts/func.inc.sh ${odir}/.data_engineering_run/"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="cp ${bin_path}/schip_processing.sh ${odir}/.data_engineering_run/"
exec_cmd ${cmd} >> ${log} 2>&1	


if [[ -d ${DOWNSTREAM_ODIR} ]]; then
  
  cmd="mkdir -p ${odir}/.data_analysis_run/"
  exec_cmd ${cmd} >> ${log} 2>&1

  cmd="cp ${R_downstream_dir}/R_scChIP_seq_analysis.R ${odir}/.data_analysis_run/"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  cmd="cp -r ${R_downstream_dir}/Modules/ ${odir}/.data_analysis_run/"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  cmd="cp ${R_downstream_dir}/report_scChIPseq_analysis.Rmd ${odir}/.data_analysis_run/"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  cmd="echo '${bin_name} ${arguments}' > ${odir}/.data_engineering_run/call.sh"
  exec_cmd ${cmd} >> ${log} 2>&1
fi
}


##Write metadata

#write_metadata ${ODIR} ${PREFIX} "${CMD_LINE}" ${CONF} ${LOGDIR} ${GENOME_BAM_FLAGGED} ${GENOME_BAM_FLAGGED_RMDUP} ${BARCODE_BAM}
write_metadata(){

    ## logs
    local log=$5/metadata.log
    echo -e "Generating metadata table ..."
    echo -e "Logs: $log"
    echo
 
    cmd="bash ${SCRIPTS_PATH}/make_metadata.sh $1/metadata $2 \"$3\" $4 $5 $6 $7 $8"
    exec_cmd ${cmd} > ${log} 2>&1
}

get_command_line(){
    line=$(history 1)
    line=${line#*[0-9]  }
    echo "$line"
}
