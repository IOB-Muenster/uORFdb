#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node 16
#SBATCH --mem=92GB
#SBATCH --time=48:00:00
#SBATCH --partition normal

# ------------------------------------------------------------------------------------------
# INFO
# ------------------------------------------------------------------------------------------
#
# AUTHOR
#
# Lynn Ogoniak
#
#
# COPYRIGHT
#
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#    
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation and/or
# other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS
# OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
# SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
# IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
# OF SUCH DAMAGE. 
# 
#
# DESCRIPTION
#
# Alignment with BWA
#
# ==========================================================================================

# ---------------------------------------------------------------------------------
# VARIABLES
# ---------------------------------------------------------------------------------
while getopts i:t:p:k: flag
do
    case "${flag}" in
        i) UUID=${OPTARG};;
        t) TYPE=${OPTARG};;
        p) PATIENT=${OPTARG};;
        k) TISSUE=${OPTARG}
    esac
done

SCRATCH_DIR="/path/to/working/dir/"
CLOUD_DIR=/path/to/cloud/dir/${TYPE}/${UUID}/

SEQKIT=/path/to/seqkit

mkdir ${SCRATCH_DIR}/Alignment_QC/${UUID}
WORKING_DIR=${SCRATCH_DIR}/Alignment_QC/${UUID}/
cd ${WORKING_DIR}

DOWNLOADED_BAM=$(ls /path/to/downloaded/BAMs/${TYPE}/${UUID}/*.bam))
REF=${SCRATCH_DIR}/Reference_genome/ncbi-genomes-2021-07-15/GCF_000001405.39_GRCh38.p13_genomic.fna

PROCESS="02_bwa"
START_TIME=$(date "+%s")

echo $TYPE
echo $UUID
echo $PATIENT $TISSUE
echo $(date "+%d/%m/%Y %H:%M:%S")
echo "--------------------------------------------------------"

# ---------------------------------------------------------------------------------
# INDEX REFERENCE GENOME
# ---------------------------------------------------------------------------------
if [ ! -e ${REF}.bwt ]; then
	echo $(date "+%d/%m/%Y %H:%M:%S") "Reference genome not indexed yet. BWA index started."
	ml palma/2019a GCC/8.2.0-2.31.1 BWA/0.7.17 
	bwa index ${REF}
	echo $(date "+%d/%m/%Y %H:%M:%S") "BWA index finished."
else
	echo $(date "+%d/%m/%Y %H:%M:%S") "Reference genome already indexed."
fi

# ---------------------------------------------------------------------------------
# CHECK IF STRANGE DATA
# ---------------------------------------------------------------------------------
ml palma/2019a GCC/8.2.0-2.31.1 BWA/0.7.17 SAMtools/1.9

READGROUPS=$(samtools view -H ${DOWNLOADED_BAM} | grep '^@RG' | awk '{gsub(/\t/,"\\t");print $0}' | awk '{gsub(" ","\\b");print $0}')
RG_NUMBER=$(samtools view -H ${DOWNLOADED_BAM} | grep '^@RG' | wc -l)

echo $(date "+%d/%m/%Y %H:%M:%S") "Found ${RG_NUMBER} read group(s):"
for READGROUP in ${READGROUPS}
do
	echo $READGROUP
	RG_ID=$(echo ${READGROUP} | sed 's/\\t/\t/g ; s/ID://1' | awk '{print $2}')
		
done

# ---------------------------------------------------------------------------------
# CALCULATE MEANREADLENGTH
# ---------------------------------------------------------------------------------
MEANREADLENGTH=""

echo $(date "+%d/%m/%Y %H:%M:%S") "Calculating average read length."

cd ${CLOUD_DIR}
FILES=$(ls *.fq.gz)
cd ${WORKING_DIR}

SUM=0
COUNTER=0

for FILE in ${FILES}; do

	if [ "${MEANREADLENGTH}" = "" ]; then
	
		FASTQ=${CLOUD_DIR}/${FILE}
		MEANREADLENGTH=$(${SEQKIT} stats ${FASTQ} | grep -v "avg_len" | awk '{print$8}')
		USED_FASTQ=${FASTQ}
		
	else
		
		echo $(date "+%d/%m/%Y %H:%M:%S") "Used ${USED_FASTQ} to determine meanreadlength: ${MEANREADLENGTH}"
		break
	
	fi
done

# ---------------------------------------------------------------------------------
# GUESS ENCODING BASED ON ALL READ GROUPS
# ---------------------------------------------------------------------------------
# 1,2,o1,o2,s are not relevant here

ENCODING="undefined"

if [ "${ENCODING}" = "undefined" ]; then

	ENCODING="aln -l"
	
	FASTQS=$(ls ${CLOUD_DIR}/*.fq.gz)
	
	for FASTQ in ${FASTQS}
	do
		
		ENCODING1=$(less ${FASTQ} | head -n1000 | \
		awk '{if(NR%4==0) printf("%s",$0);}' | \
		od -A n -t u1 | \
		awk 'BEGIN{min=100;max=0;}{for(i=1;i<=NF;i++) {if($i>max) max=$i; if($i<min) min=$i;}}END{if(max<=74 && min<59) print "aln"; else if(max>73 && min>=64) print "aln -l"; else if(min>=59 && min<64 && max>73) print "aln"; else print "aln";}')
		
		echo $(date "+%d/%m/%Y %H:%M:%S") "${FASTQ}: Encoding was guessed. Fastq1 encoding: ${ENCODING1}."

		if [ "${ENCODING1}" = "aln" ]; then
			ENCODING="aln"
		fi
	
	done
	
	echo $(date "+%d/%m/%Y %H:%M:%S") "Encoding was set to ${ENCODING} for all read groups."
	
fi

# ---------------------------------------------------------------------------------
# PROCESS EACH READGROUP
# ---------------------------------------------------------------------------------
for READGROUP in ${READGROUPS}; do

	RG_ID=$(echo ${READGROUP} | sed 's/\\t/\t/g ; s/ID://1' | awk '{print $2}')
	echo ""
	echo ${RG_ID}
	
	# ---------------------------------------------------------------------------------
	# CHECK IF ALREADY PROCESSED
	# ---------------------------------------------------------------------------------
	if [ -e ${CLOUD_DIR}/${RG_ID}_aligned.bam ]; then
		echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}: File exists already in Cloud. ${PROCESS} canceled."
		PROCESS="02_bwa_subset"
		#exit 0
	
	elif [ -e ${CLOUD_DIR}/${UUID}_qc.bam ]; then
		echo $(date "+%d/%m/%Y %H:%M:%S") "Alignment Pipeline seems to be finished already. QC file exists already in Cloud. ${PROCESS} canceled."
		exit 1
	
	elif [ -e ${WORKING_DIR}/${RG_ID}_aligned.bam ] || [ -e ${WORKING_DIR}/intermediates/${RG_ID}_aligned.bam ]; then
		echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}: File exists already in Scratch. ${PROCESS} canceled. Check for errors in previous processes."
		echo $(date "+%d/%m/%Y %H:%M:%S")",${PROCESS},${UUID},${RG_ID},ERROR_already_in_scratch" >> ${SCRATCH_DIR}/bwa_logs/bwa_error_messages.csv
		exit 1
		
	else
		
		echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}: File does not exist yet. ${PROCESS} starts."
		echo $(date "+%d/%m/%Y %H:%M:%S") "${PROCESS} started."
		echo "Read group (line): ${READGROUP}"
		echo "Read group ID: ${RG_ID}"
		
		# ---------------------------------------------------------------------------------
		# CHECK DISK SPACE
		# ---------------------------------------------------------------------------------
		DISK_SPACE=$(du -sb ${SCRATCH_DIR}/ | awk '{print $1}')
		SPACE_MAX=4000000000000
		
		while [ ${DISK_SPACE} -ge ${SPACE_MAX} ]; do
			echo $(date "+%d/%m/%Y %H:%M:%S") "Not enough disk space (${DISK_SPACE} B). Try again in 15 min."
			sleep 450
			DISK_SPACE=$(du -sb ${SCRATCH_DIR}/ | awk '{print $1}')
		done
		
		echo $(date "+%d/%m/%Y %H:%M:%S") "Disk space sufficient (${DISK_SPACE} B occupied). Operation started."
		
		# ---------------------------------------------------------------------------------
		# BWA ALIGNMENT
		# ---------------------------------------------------------------------------------
		
		# https://ucdavis-bioinformatics-training.github.io/2017-August-Variant-Analysis-Workshop/wednesday/alignment.html
		
		FASTQ1=${CLOUD_DIR}/${RG_ID}_1.fq.gz
		FASTQ2=${CLOUD_DIR}/${RG_ID}_2.fq.gz
		FASTQS=${CLOUD_DIR}/${RG_ID}_s.fq.gz
		
		NUM_SEQS_1="0"
		NUM_SEQS_S="0"
		
		#check if both s and 1+2 exist, decide which one contains the reads
		if [ -e ${FASTQS} ] && [ -e ${FASTQ1} ] && [ -e ${FASTQ2} ]; then
		
			echo "Read group with single end data found. Checking if paired end data exists for this read group."
			
			SEQKIT_1=$(${SEQKIT} stats ${FASTQ1} | grep -v "avg_len")
			NUM_SEQS_1=$(echo ${SEQKIT_1} | awk '{print $4}')
			
			SEQKIT_S=$(${SEQKIT} stats ${FASTQS} | grep -v "avg_len")
			NUM_SEQS_S=$(echo ${SEQKIT_S} | awk '{print $4}')
			
			echo "Num_seqs in ${FASTQ1}: ${NUM_SEQS_1}."
			echo "Num_seqs in ${FASTQS}: ${NUM_SEQS_S}."
			
		fi
		
		if [ -e ${FASTQS} ] && [ "${NUM_SEQS_1}" == "0" ]; then
		
			if [ ${MEANREADLENGTH} -ge 70 ]; then
		
				PROCESS="02_bwa-mem_s"
				echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}: BWA-MEM S started."
				
				bwa mem \
				-t 16 \
				-T 0 \
				-R ${READGROUP} \
				${REF} \
				${FASTQS} | \
				samtools view -Shb \
				-o ${WORKING_DIR}/${RG_ID}_aligned.bam || \
					{ echo "BWA-MEM failed"; exit 1; }
			
			elif [ ${MEANREADLENGTH} -lt 70 ]; then
				
				#http://bio-bwa.sourceforge.net/bwa.shtml
				PROCESS="02_bwa-aln_s"
				echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}: BWA-ALN S started."
				
				bwa ${ENCODING} -t 16 -R ${READGROUP} ${REF} ${FASTQS} > ${WORKING_DIR}/${RG_ID}_s.sai && \
				bwa samse -r ${READGROUP} ${REF} ${WORKING_DIR}/${RG_ID}_s.sai ${FASTQS} | \
				samtools view -Shb -o ${WORKING_DIR}/${RG_ID}_aligned.bam || \
					{ echo "BWA-ALN failed"; exit 1; }
			
			fi
		
		#normal paired end data
		elif [ -e ${FASTQ1} ] && [ -e ${FASTQ2} ] && [ "${NUM_SEQS_S}" == "0" ]; then
			
			if [ ${MEANREADLENGTH} -ge 70 ]; then
			
				PROCESS="02_bwa-mem"
				echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}: BWA-MEM started."
				
				bwa mem \
				-t 16 \
				-T 0 \
				-R ${READGROUP} \
				${REF} \
				${FASTQ1} \
				${FASTQ2} | \
				samtools view -Shb \
				-o ${WORKING_DIR}/${RG_ID}_aligned.bam || \
					{ echo "BWA-MEM failed"; exit 1; }
				
			elif [ ${MEANREADLENGTH} -lt 70 ]; then
				
				PROCESS="02_bwa-aln"
				echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}: BWA-ALN started."
				
				bwa ${ENCODING} -t 16 -R ${READGROUP} ${REF} ${FASTQ1} > ${WORKING_DIR}/${RG_ID}_1.sai && \
				bwa ${ENCODING} -t 16 -R ${READGROUP} ${REF} ${FASTQ2} > ${WORKING_DIR}/${RG_ID}_2.sai && \
				bwa sampe -r ${READGROUP} ${REF} ${WORKING_DIR}/${RG_ID}_1.sai ${WORKING_DIR}/${RG_ID}_2.sai ${FASTQ1} ${FASTQ2} | \
				samtools view -Shb -o ${WORKING_DIR}/${RG_ID}_aligned.bam || \
					{ echo "BWA-ALN failed"; exit 1; }
			
			fi
		
		else
			echo "Check data. Paired end and single data contain reads."
			exit 1
			
		fi
		
		# ---------------------------------------------------------------------------------
		# CHECK IF ALIGNMENT COMPLETED
		# ---------------------------------------------------------------------------------
		if [ -e ${WORKING_DIR}/${RG_ID}_aligned.bam ]; then
			VALIDATE=$(samtools quickcheck -v ${WORKING_DIR}/${RG_ID}_aligned.bam)
			
			if [ "${VALIDATE}" = "" ]; then
				echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}: Validated ${WORKING_DIR}/${RG_ID}_aligned.bam."
				
				mkdir ${WORKING_DIR}/intermediates
				mv ${WORKING_DIR}/${RG_ID}* ${WORKING_DIR}/intermediates/
				echo "${TYPE}" > ${SCRATCH_DIR}/bwa_logs/intermediates_info/${UUID}
				
			else
				echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}: Samtools validation message:"
				echo ${VALIDATE}
				exit 1
			fi
		else
			echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}: Error. Alignment file not in Scratch."
			exit 1
		fi
		
	fi
	
done

echo $(date "+%d/%m/%Y %H:%M:%S") "${PROCESS} completed for all read groups. Checking if uploaded to Cloud."

# ---------------------------------------------------------------------------------
# CHECK RSYNC
# ---------------------------------------------------------------------------------
for READGROUP in ${READGROUPS}; do

	RG_ID=$(echo ${READGROUP} | sed 's/\\t/\t/g ; s/ID://1' | awk '{print $2}')
	
	echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}: Checking if uploaded to Cloud."
	
	while [ -e ${WORKING_DIR}/intermediates/${RG_ID}_aligned.bam ]; do
		sleep 180
	done
	
	if [ -e ${CLOUD_DIR}/${RG_ID}_aligned.bam ]; then
		echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}: Successfully uploaded to Cloud."
	else
		echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}: Missing. Check for errors."
		echo $(date "+%d/%m/%Y %H:%M:%S")",${PROCESS},${UUID},${RG_ID},ERROR_missing" >> ${SCRATCH_DIR}/bwa_logs/bwa_error_messages.csv
		exit 1
	fi

done

# ---------------------------------------------------------------------------------
# TIME LOG
# ---------------------------------------------------------------------------------
END_TIME=$(date "+%s")
RUN_TIME=$((END_TIME-START_TIME))
echo "${PROCESS},${RUN_TIME}" >> ${SCRATCH_DIR}/bwa_logs/runtimes/${UUID}_runtimes.log
			
echo "End of script."
