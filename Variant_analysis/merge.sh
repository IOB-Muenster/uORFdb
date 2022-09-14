#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node 16
#SBATCH --mem=16GB
#SBATCH --time=24:00:00
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
# Picard Merge
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

mkdir ${SCRATCH_DIR}/Alignment_QC/${UUID}
WORKING_DIR=${SCRATCH_DIR}/Alignment_QC/${UUID}/
cd ${WORKING_DIR}

DOWNLOADED_BAM=$(ls /path/to/downloaded/BAMs/${TYPE}/${UUID}/*.bam)
MERGED_BAM=${WORKING_DIR}/${UUID}_merged.bam

PICARD="/path/to/picard/picard.jar"

PROCESS="05_merge"
START_TIME=$(date "+%s")

echo $TYPE
echo $UUID
echo $PATIENT $TISSUE
echo $(date "+%d/%m/%Y %H:%M:%S")
echo "--------------------------------------------------------"

# ---------------------------------------------------------------------------------
# CHECK IF ALREADY MERGED
# ---------------------------------------------------------------------------------
if [ -e ${CLOUD_DIR}/${UUID}_merged.bam ]; then
	echo $(date "+%d/%m/%Y %H:%M:%S") "File exists already in Cloud. ${PROCESS} canceled."
	exit 0
	
elif [ -e ${CLOUD_DIR}/${UUID}_qc.bam ]; then
	echo $(date "+%d/%m/%Y %H:%M:%S") "Alignment Pipeline seems to be finished already. QC file exists already in Cloud. ${PROCESS} canceled."
	exit 1

elif [ -e ${MERGED_BAM} ] || [ -e ${WORKING_DIR}/intermediates/${UUID}_merged.bam ]; then
	echo $(date "+%d/%m/%Y %H:%M:%S") "File exists already in Scratch. ${PROCESS} canceled. Check for errors in previous processes."
	echo $(date "+%d/%m/%Y %H:%M:%S")",${PROCESS},${UUID},ERROR_already_in_scratch" >> ${SCRATCH_DIR}/bwa_logs/bwa_error_messages.csv
	exit 1
	
else
	echo $(date "+%d/%m/%Y %H:%M:%S") "File does not exist yet."
fi

# ---------------------------------------------------------------------------------
# CHECK IF MERGING NEEDED
# ---------------------------------------------------------------------------------
ml palma/2020b  GCC/10.2.0 SAMtools/1.11

RG_IDS=$(samtools view -H ${DOWNLOADED_BAM} | grep "^@RG" | awk '{print $2}' | sed 's/ID://1')
RG_NUMBER=$(samtools view -H ${DOWNLOADED_BAM} | grep "^@RG" | wc -l)
echo $(date "+%d/%m/%Y %H:%M:%S") "Found ${RG_NUMBER} read group(s): ${READGROUPS}"

if [ ${RG_NUMBER} -gt 1 ]; then
	echo $(date "+%d/%m/%Y %H:%M:%S") "Sorted files have to be merged."
	
	#check if all read groups are there
	for RG_ID in ${RG_IDS}
	do
		echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}."
		
		if [ -e ${CLOUD_DIR}/${RG_ID}_sorted.bam ]; then
			echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}: ${CLOUD_DIR}/alignment_qc/${RG_ID}_sorted.bam found."
		else
			echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}: Missing. Check previous processes."
			echo $(date "+%d/%m/%Y %H:%M:%S")",${PROCESS},${UUID},${RG_ID},ERROR_missing" >> ${SCRATCH_DIR}/bwa_logs/bwa_error_messages.csv
			exit 1
		fi
	done

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
	# MERGE WITH PICARD TOOLS
	# ---------------------------------------------------------------------------------
	ml palma/2020b Java/11.0.2
	
	cd ${CLOUD_DIR}
	BAMFILES=$(ls -1 ${CLOUD_DIR}/*_sorted.bam | sed "s/^/-I /")
	cd ${WORKING_DIR}
	
	java -jar ${PICARD} MergeSamFiles \
		${BAMFILES} \
		-O ${MERGED_BAM} \
		-ASSUME_SORTED false \
		-CREATE_INDEX true \
		-SO coordinate \
		-MSD false \
		-USE_THREADING true \
		-VALIDATION_STRINGENCY STRICT || \
    		{ echo "Merging failed"; exit 1; }
    		
    # ---------------------------------------------------------------------------------
	# CHECK IF MERGING COMPLETED
	# ---------------------------------------------------------------------------------
	ml palma/2020b  GCC/10.2.0 SAMtools/1.11
	
	if [ -e ${MERGED_BAM} ]; then
		VALIDATE=$(samtools quickcheck -v ${MERGED_BAM})
		
		if [ "${VALIDATE}" = "" ]; then
			echo $(date "+%d/%m/%Y %H:%M:%S") "Validated ${MERGED_BAM}."
	
			mkdir ${WORKING_DIR}/intermediates
			mv ${MERGED_BAM} ${WORKING_DIR}/intermediates/
			mv ${WORKING_DIR}/*_merged.bai ${WORKING_DIR}/intermediates/
			echo "${TYPE}" > ${SCRATCH_DIR}/bwa_logs/intermediates_info/${UUID}
			
		else
			echo $(date "+%d/%m/%Y %H:%M:%S") "Samtools validation message:"
			echo ${VALIDATE}
			exit 1
		fi
	else
		echo $(date "+%d/%m/%Y %H:%M:%S") "Error. Merged file not in Scratch."
		exit 1
	fi
	
	# ---------------------------------------------------------------------------------
	# CHECK RSYNC
	# ---------------------------------------------------------------------------------
	echo $(date "+%d/%m/%Y %H:%M:%S") "Checking if uploaded to Cloud."
	
	while [ -e ${WORKING_DIR}/intermediates/${UUID}_merged.bam ]; do
		sleep 180
	done
	
	if [ -e ${CLOUD_DIR}/${UUID}_merged.bam ]; then
		echo $(date "+%d/%m/%Y %H:%M:%S") "Successfully uploaded to Cloud."
	else
		echo $(date "+%d/%m/%Y %H:%M:%S") "Missing. Check for errors."
		echo $(date "+%d/%m/%Y %H:%M:%S")",${PROCESS},${UUID},ERROR_missing" >> ${SCRATCH_DIR}/bwa_logs/bwa_error_messages.csv
		exit 1	
	fi
	  
else
	echo $(date "+%d/%m/%Y %H:%M:%S") "Only one read group. Merging is skipped."
	PROCESS="05_merging_skipped"
	
fi

# ---------------------------------------------------------------------------------
# TIME LOG
# ---------------------------------------------------------------------------------
END_TIME=$(date "+%s")
RUN_TIME=$((END_TIME-START_TIME))
echo "${PROCESS},${RUN_TIME}" >> ${SCRATCH_DIR}/bwa_logs/runtimes/${UUID}_runtimes.log

echo "End of script."
