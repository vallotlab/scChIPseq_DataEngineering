## Mapping Research Pipeline
## Copyleft 2018 Institut Curie
## Author(s): Nicolas Servant, Pacôme Prompsy
## Contact: nicolas.servant@curie.fr; pacome.prompsy@curie.fr
## This software is distributed without any guarantee under the terms of the CECILL License
## See the LICENCE file for details


##GetConf function <- from EWOK (https://gitlab.curie.fr/data-analysis/ewok)
getConf_func(){
if [[ $# != "6" && $# != "7" ]]; then echo "All arguments must be filled : getConf_func TEMPLATE CONFIGS DESIGN_TYPE ASSEMBLY OUT_FILE MARK [TARGET_BED] "; echo; help_func GetConf;
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
	cp ${TEMPLATE} ${OUT_FILE}

	for i in $(seq 1 $(awk '{FS="\t"; OFS="\t"}{if(NR == 1) print NF}' ${CONFIGS})); do 
		KEY=$(awk -v i=$i '{FS="\t"; OFS="\t"}{if(NR == 1)print "{"$i"}"}' ${CONFIGS})
		VALUE=$(awk -v i=$i -v c=${CONFIG_LINE} '{FS="\t"; OFS="\t"}{if(NR == c) print $i}' ${CONFIGS})
		if [[ $(echo $VALUE |grep -c " ") > 0 ]] ; then 
			sed -i "s|${KEY}|\"${VALUE}\"|g" ${OUT_FILE} 
		else 
			sed -i "s|${KEY}|${VALUE}|g" ${OUT_FILE} 
		fi 
	done 
  sed -i 's/"//g' ${OUT_FILE}
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
    		cmd="${FASTX_PATH}fastx_trimmer -Q 33 -l $2 -i <(gzip -cd $1) -o ${ofile}"
    		exec_cmd ${cmd} > ${log} 2>&1	
		else			#Change option -l (last) to -f (first) to remove the barcode and linker from R2 reads in preparation for genome alignment
		local ofile=${out_prefix}_trimmed_G.R2.fastq
		cmd="${FASTX_PATH}fastx_trimmer -Q 33 -f $2 -i <(gzip -cd $1) -o ${ofile}" 
    	exec_cmd ${cmd} > ${log} 2>&1	
	fi

    cmd="gzip ${ofile}"
    exec_cmd ${cmd} > ${log} 2>&1
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
           cmd_in="<(gzip -cd ${inputs[0]}) <(gzip -cd ${inputs[1]})"
       else
           cmd_in="${inputs[0]} ${inputs[1]}"
       fi
   fi

   ## sample_id
   if [ ! -z ${SAMPLE_ID} ]; then
       cmd_id="--outSAMattrRGline ID:${SAMPLE_ID} PL:ILLUMINA"
   else
       cmd_id=""
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
  if [[ ${BARCODE_LENGTH} -eq 56 ]]
  	then
  	echo -e "Short barcode :LBC"
      start_index_1=1
  	  start_index_2=21
  	  start_index_3=41
  	  size_index=16
         else
  	echo -e "Long barcode : Hifibio "
      start_index_1=1
    	start_index_2=25
    	start_index_3=49
    	size_index=20
  fi
  
  ##Extract three indexes from reads : 1 - 16 = index 1 ; 21 - 36 = index 2; 41 - 56 = index 3
  cmd="gzip -cd  $read2 | awk -v start_index_1=$start_index_1 -v size_index=$size_index  'NR%4==1{print \">\"substr(\$0,2)}; NR%4==2{print substr(\$0,start_index_1,size_index)}' > ${out}/read_indexes_1.fasta"
  exec_cmd ${cmd} > ${log} 2>&1
  
  cmd="gzip -cd  $read2 | awk -v start_index_2=$start_index_2 -v size_index=$size_index 'NR%4==1{print  \">\"substr(\$0,2)}; NR%4==2{print substr(\$0,start_index_2,size_index)}' > ${out}/read_indexes_2.fasta"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  cmd="gzip -cd  $read2 | awk -v start_index_3=$start_index_3 -v size_index=$size_index 'NR%4==1{print \">\"substr(\$0,2)}; NR%4==2{print substr(\$0,start_index_3,size_index)}' > ${out}/read_indexes_3.fasta"
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
        local cmd="bamCoverage --bam $1 --outFileName $out/${name}.bw --numberOfProcessors ${NB_PROC} --normalizeUsing RPKM --blackListFileName ${ENCODE_BLACKLIST}"
    else
        local cmd="bamCoverage --bam $1 --outFileName $out/${name}.bw --numberOfProcessors ${NB_PROC} --normalizeUsing RPKM "
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

  #Remove secondary aligned reads (256 <=> "not primary alignment") & If R1 is unmapped or multimapped (NH != 1), tag R1 & R2 with flag "4" <=> "unmapped" & "chr" = '*'
  cmd="samtools view -F 256 $1 | awk -v OFS='\t' 'NR%2==1{if(\$12==\"NH:i:1\"){mapped=1;print \$0} else{mapped=0;\$2=4;\$3=\"*\";\$4=0;\$6=\"*\";print \$0}} NR%2==0{if(mapped==1){print \$0} else{\$2=4;\$3=\"*\";\$4=0;\$6=\"*\";print \$0} }' > ${out_prefix}.sam"
  exec_cmd ${cmd} > ${log} 2>&1
  
  #If read is mapped R1 & unmapped R2 -> set R2 position as '2147483647'
  cmd="cat ${out_prefix}.sam | awk -v OFS='\t' 'NR%2==1{print \$0} NR%2==0{if(\$3==\"*\"){\$4=2147483647;print \$0} else{print \$0} }' > ${out_prefix}_2.sam"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  #Remove comments from the header that produce bugs in the count phase
  cmd="samtools view -H $1 | sed '/^@CO/ d' > ${out_prefix}_header.sam"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  cmd="cat ${out_prefix}_2.sam >> ${out_prefix}_header.sam && mv ${out_prefix}_header.sam ${out_prefix}.sam && samtools view -b -@ ${NB_PROC} ${out_prefix}.sam > ${out_prefix}_unique.bam"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  cmd="rm -f ${out_prefix}_2.sam ${out_prefix}.sam"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  #Keeping R1 aligned + R2 start as tag 'XS' (Switch from Paired End Bam to Single End Bam)
  cmd="samtools view ${out_prefix}_unique.bam | awk '{OFS = \"\t\" ; if(NR%2==1 && !(\$3==\"*\")) {R1=\$0} else if(NR%2==1){R1=0}; if(NR%2==0 && !(R1==0)){tagR2Seq=\"XD:Z:\"\$10; tagR2Pos=\"XS:i:\"\$4;print R1,tagR2Pos,tagR2Seq}}' > ${out_prefix}_unique.sam"
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
  exec_cmd ${cmd} >> ${log} 2>&1
  
  echo "The barcode field is $barcode_field" >> ${log}
  
  #Find the column containing the position R2 tag XS
  cmd="posR2_field=\$(samtools view ${out_prefix}_flagged.bam  | sed -n \"1 s/XS.*//p\" | sed 's/[^\t]//g' | wc -c)"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  echo "The position R2 field is $posR2_field" >> ${log}
  
  cmd="printf '@HD\tVN:1.4\tSO:unsorted\n' > ${out_prefix}_header.sam"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  cmd="samtools view -H ${out_prefix}_flagged.bam | sed '/^@HD/ d' >> ${out_prefix}_header.sam"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  #Sort by barcode then chromosome then position R2 then Position R1 (for the PCR removal) 
  #It is important to sort by R1 pos also	because the removal is done by comparing consecutive lines ! 
  cmd="samtools view ${out_prefix}_flagged.bam | LC_ALL=C sort -T ${TMP_DIR} --parallel=${NB_PROC} -t $'\t' -k \"$barcode_field.8,$barcode_field\"n -k 3.4,3g -k \"$posR2_field.6,$posR2_field\"n -k 4,4n >> ${out_prefix}_header.sam && samtools view -@ ${NB_PROC} -b ${out_prefix}_header.sam > ${out_prefix}_flagged.sorted.bam"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  #Create count Table from flagged (already sorted by barcode)
  cmd="samtools view ${out_prefix}_flagged.sorted.bam | awk -v bc_field=$barcode_field '{print substr(\$bc_field,6)}' |  uniq -c > ${out_prefix}_flagged.count"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  if [ ${REMOVE_BY} == 'position' ]
  then
  echo -e "Removing PCR duplicates using R2 POSITION"
  #Remove PCR duplicates = read having same barcode, R1 position, same R2 POSITION, same chr ("exactly equal")
  cmd="samtools view ${out_prefix}_flagged.sorted.bam | awk -v bc_field=$barcode_field -v R2_field=$posR2_field -v out=${out}/ 'BEGIN{countR1unmappedR2=0;countPCR=0};NR==1{print \$0;lastChrom=\$3;lastBarcode=\$bc_field; split( \$R2_field,lastR2Pos,\":\"); lastR1Pos=\$4} ; NR>=2{split( \$R2_field,R2Pos,\":\");R1Pos=\$4; if(R2Pos[3]==2147483647){print \$0;countR1unmappedR2++; next}; if( (R1Pos==lastR1Pos) && (R2Pos[3]==lastR2Pos[3]) && ( \$3==lastChrom ) && (\$bc_field==lastBarcode) ){countPCR++;next} {print \$0;lastR1Pos=\$4;lastChrom=\$3;lastBarcode=\$bc_field; split( \$R2_field,lastR2Pos,\":\") }} END {print countPCR > out\"count_PCR_duplicates\";print countR1unmappedR2 > out\"countR1unmappedR2\"}' > ${out_prefix}_flagged_rmPCR.sam"
  exec_cmd ${cmd} >> $log 2>&1
  else
  echo -e "Removing PCR duplicates using R2 SEQUENCE"
  cmd="seqR2_field=\$(samtools view ${out_prefix}_flagged.bam  | sed -n \"1 s/XD.*//p\" | sed 's/[^\t]//g' | wc -c)"
  exec_cmd ${cmd} >> ${log} 2>&1
  #Remove PCR duplicates = read having same barcode, R1 position, same R2 SEQUENCE, same chr
  cmd="samtools view ${out_prefix}_flagged.sorted.bam | awk -v bc_field=$barcode_field -v R2_field=$seqR2_field -v out=${out}/ 'BEGIN{countR1unmappedR2=0;countPCR=0};NR==1{print \$0;lastChrom=\$3;lastBarcode=\$bc_field; split( \$R2_field,lastR2Seq,\":\"); lastR1Pos=\$4} ; NR>=2{split( \$R2_field,R2Seq,\":\");R1Pos=\$4; if( (R1Pos==lastR1Pos) && (R2Seq[3]==lastR2Seq[3]) && ( \$3==lastChrom ) && (\$bc_field==lastBarcode) ){countPCR++;next} {print \$0;lastR1Pos=\$4;lastChrom=\$3;lastBarcode=\$bc_field; split( \$R2_field,lastR2Seq,\":\") }} END {print countPCR > out\"count_PCR_duplicates\";print countR1unmappedR2 > out\"countR1unmappedR2\"}' > ${out_prefix}_flagged_rmPCR.sam"
  exec_cmd ${cmd} >> $log 2>&1
  fi
  
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
  
  if [ ${REMOVE_RT_DUPLICATES} == 'TRUE' ] 
  then
  echo -e "Removing RT duplicates"
  #Remove RT duplicates (if two consecutive reads have the same barcode and same R2 chr&start) but not same R1 
  cmd="cat ${out_prefix}_flagged_rmPCR.sam | awk -v bc_field=$barcode_field -v R2_field=$posR2_field -v out=${out}/ 'BEGIN{count=0};NR==1{print \$0;lastChrom=\$3;lastBarcode=\$bc_field; split( \$R2_field,lastR2Pos,\":\")} ; NR>=2{split( \$R2_field,R2Pos,\":\");if((R2Pos[3]==lastR2Pos[3]) && (R2Pos[3]!=2147483647) && (lastR2Pos[3]!=2147483647)  && ( \$3==lastChrom ) && (\$bc_field==lastBarcode) ){count++;next} {print \$0;lastChrom=\$3;lastBarcode=\$bc_field; split( \$R2_field,lastR2Pos,\":\") }} END {print count > out\"count_RT_duplicates\"}' > ${out_prefix}_flagged_rmPCR_RT.sam"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  cmd="samtools view -H ${out_prefix}_flagged.bam  | sed '/^@CO/ d' > ${out_prefix}_header.sam"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  cmd="cat ${out_prefix}_flagged_rmPCR_RT.sam >> ${out_prefix}_header.sam && samtools view -@ ${NB_PROC} -b ${out_prefix}_header.sam > ${out_prefix}_flagged_rmPCR_RT.bam"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  #Create count Table from flagged - PCR dups - RT dups  (already sorted by barcode)
  cmd="samtools view ${out_prefix}_flagged_rmPCR_RT.bam | awk -v bc_field=$barcode_field '{print substr(\$bc_field,6)}' | uniq -c > ${out_prefix}_flagged_rmPCR_RT.count"
  exec_cmd ${cmd} >> ${log} 2>&1
  
  ## Sort flagged_rmPCR_RT file
  cmd="samtools sort -@ ${NB_PROC} ${out_prefix}_flagged_rmPCR_RT.bam > ${out_prefix}_flagged_rmPCR_RT_sorted.bam"
  exec_cmd ${cmd} >> $log 2>&1
  
  ## Rename flagged_rmPCR_RT file
  cmd="mv ${out_prefix}_flagged_rmPCR_RT_sorted.bam ${out_prefix}_flagged_rmPCR_RT.bam"
  exec_cmd ${cmd} >> $log 2>&1
  
  else
  echo -e "Not removing RT duplicates"
  ## Copy flagged_rmPCR to flagged_rmPCR_RT
  cmd="cp ${out_prefix}_flagged_rmPCR.bam ${out_prefix}_flagged_rmPCR_RT.bam"
  exec_cmd ${cmd} >> $log 2>&1
  cmd="cp ${out_prefix}_flagged_rmPCR.count ${out_prefix}_flagged_rmPCR_RT.count"
  exec_cmd ${cmd} >> $log 2>&1
  ## Set RT duplicate count to 0
  cmd="echo 0 > ${out}/count_RT_duplicates"
  exec_cmd ${cmd} >> $log 2>&1
  fi
  
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
  cmd="n_R1_mapped_R2_unmapped=\$(cat ${out}/countR1unmappedR2)"
  exec_cmd ${cmd} >> $log 2>&1
  
  ## Rename flagged.sorted -> flagged
  cmd="mv ${out_prefix}_flagged.sorted.bam ${out_prefix}_flagged.bam"
  exec_cmd ${cmd} >> $log 2>&1
 
  cmd="n_unique_except_R1_unmapped_R2=\$(($n_mapped_barcoded - $n_pcr_duplicates - $n_rt_duplicates))"
  exec_cmd ${cmd} >> $log 2>&1
  
  echo "## Number of reads mapped and barcoded: $n_mapped_barcoded" >> ${log}
  echo "## Number of pcr duplicates: $n_pcr_duplicates" >> ${log}
  echo "## Number of rt duplicates: $n_rt_duplicates" >> ${log}
  echo "## Number of R1 mapped but R2 unmapped: $n_R1_mapped_R2_unmapped" >> ${log}
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
    exec_cmd ${cmd} >> ${log} 2>&1
    #Create count Table from flagged - PCR dups - RT dups and window-based rmDup (need to sort by b arcode)
    cmd="barcode_field=\$(samtools view ${prefix}_rmDup.bam  | sed -n \"1 s/XB.*//p\" | sed 's/[^\t]//g' | wc -c)"
    exec_cmd ${cmd} >> $log 2>&1
    cmd="samtools view ${prefix}_rmDup.bam | awk -v bc_field=$barcode_field '{print substr(\$bc_field,6)}' | sort | uniq -c > ${prefix}_rmDup.count"		
    exec_cmd ${cmd} >> $log 2>&1

    ## Index BAM file
    cmd="samtools index ${prefix}_rmDup.bam"
    exec_cmd ${cmd} >> $log 2>&1
    
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
  exec_cmd ${cmd} >> ${log} 2>&1
    
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
    exec_cmd ${cmd} >> ${log} 2>&1

    echo "Barcodes found = $barcodes" >> ${log} 
    for bsize in ${BIN_SIZE}
    do
	echo -e "at ${bsize} resolution ..."
	opts="-b ${bsize} "
        if [ ! -z ${MIN_COUNT_PER_BARCODE_AFTER_RMDUP} ]; then
	    opts="${opts} -f ${MIN_COUNT_PER_BARCODE_AFTER_RMDUP} "
	fi
        cmd="${PYTHON_PATH}/python ${SCRIPTS_PATH}/sc2counts.py -i $1 -o ${prefix}_counts_${bsize}.tsv ${opts} -s $barcodes -v"
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
        cmd="${PYTHON_PATH}/python ${SCRIPTS_PATH}/sc2counts.py -i $1 -o ${prefix}_counts_${osuff}.tsv ${opts} -s $barcodes -v"
        exec_cmd ${cmd} >> ${log} 2>&1
    done
    echo

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

cmd="grep '##' ${logdir}/barcode/index_mapping_bowtie2.log >> ${odir}/${prefix}_scChIPseq_logs.txt"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="echo \"Barcode Flagging And Exact Duplicate Removal\"  >> ${odir}/${prefix}_scChIPseq_logs.txt"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="grep '##' ${logdir}/remove_PCR_RT_duplicates_func.log >> ${odir}/${prefix}_scChIPseq_logs.txt"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="echo \"Duplicate Removal Window (rmDup)\"  >> ${odir}/${prefix}_scChIPseq_logs.txt"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="grep '## Number' ${logdir}/rmDup.log >> ${odir}/${prefix}_scChIPseq_logs.txt"
exec_cmd ${cmd} >> ${log} 2>&1

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
