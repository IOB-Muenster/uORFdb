#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node 16
#SBATCH --mem=16GB
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
# Picard CleanSam
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

PICARD="/path/to/picard/picard.jar"

PROCESS="03_clean"
START_TIME=$(date "+%s")

echo $TYPE
echo $UUID
echo $PATIENT $TISSUE
echo $(date "+%d/%m/%Y %H:%M:%S")
echo "--------------------------------------------------------"

# ---------------------------------------------------------------------------------
# ITERATE OVER READ GROUP IDs
# ---------------------------------------------------------------------------------
ml palma/2019a GCC/8.2.0-2.31.1 SAMtools/1.9

READGROUPS=$(samtools view -H ${DOWNLOADED_BAM} | grep "^@RG" | awk '{gsub(/\t/,"\\t");print $0}' | awk '{gsub(" ","\\b");print $0}')
RG_NUMBER=$(samtools view -H ${DOWNLOADED_BAM} | grep "^@RG" | wc -l)
echo $(date "+%d/%m/%Y %H:%M:%S") "Found ${RG_NUMBER} read group(s): ${READGROUPS}"

for READGROUP in ${READGROUPS}
do
	RG_ID=$(echo ${READGROUP} | sed 's/\\t/\t/g ; s/ID://1' | awk '{print $2}')

	if [ -e ${CLOUD_DIR}/${RG_ID}_cleaned.bam ]; then
		echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}: File exists already in Cloud. ${PROCESS} canceled."
		PROCESS="03_clean_subset"
		#exit 0
		
	elif [ -e ${CLOUD_DIR}/${UUID}_qc.bam ]; then
		echo $(date "+%d/%m/%Y %H:%M:%S") "Alignment Pipeline seems to be finished already. QC file exists already in Cloud. ${PROCESS} canceled."
		exit 1
	
	elif [ -e ${WORKING_DIR}/${RG_ID}_cleaned.bam ] || [ -e ${WORKING_DIR}/intermediates/${RG_ID}_cleaned.bam ]; then
		echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}: File exists already in Scratch. ${PROCESS} canceled. Check for errors in previous processes."
		echo $(date "+%d/%m/%Y %H:%M:%S")",${PROCESS},${UUID},${RG_ID},ERROR_already_in_scratch" >> ${SCRATCH_DIR}/bwa_logs/bwa_error_messages.csv
		exit 1
		
	else
		echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}: File does not exist yet. ${PROCESS} starts."
				
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
		# CLEAN WITH PICARD TOOLS
		# ---------------------------------------------------------------------------------	
		ml palma/2020b Java/11.0.2
		
		java -jar ${PICARD} CleanSam \
		-INPUT ${CLOUD_DIR}/${RG_ID}_aligned.bam \
		-OUTPUT ${WORKING_DIR}/${RG_ID}_cleaned.bam || \
			{ echo "Cleaning failed"; exit 1; }
		
		# ---------------------------------------------------------------------------------
		# CHECK IF CLEANING COMPLETED
		# ---------------------------------------------------------------------------------
		ml palma/2019a GCC/8.2.0-2.31.1 SAMtools/1.9
		
		if [ -e ${WORKING_DIR}/${RG_ID}_cleaned.bam ]; then
			VALIDATE=$(samtools quickcheck -v ${WORKING_DIR}/${RG_ID}_cleaned.bam)
			
			if [ "${VALIDATE}" = "" ]; then
				echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}: Validated ${WORKING_DIR}/${RG_ID}_cleaned.bam."
				
				mkdir ${WORKING_DIR}/intermediates
				mv ${WORKING_DIR}/${RG_ID}_cleaned.bam ${WORKING_DIR}/intermediates/
				echo "${TYPE}" > ${SCRATCH_DIR}/bwa_logs/intermediates_info/${UUID}
			
			else
				echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}: Samtools validation message:"
				echo ${VALIDATE}
				exit 1
			fi
		else
			echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}: Error. Cleaned bam file not in Scratch."
			exit 1
		fi
		
	fi
done

echo $(date "+%d/%m/%Y %H:%M:%S") "${PROCESS} completed for all read groups. Checking if uploaded to Cloud."

# ---------------------------------------------------------------------------------
# CHECK RSYNC
# ---------------------------------------------------------------------------------
for READGROUP in ${READGROUPS}
do
	RG_ID=$(echo ${READGROUP} | sed 's/\\t/\t/g ; s/ID://1' | awk '{print $2}')
	
	echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}: Checking if uploaded to Cloud."
	
	while [ -e ${WORKING_DIR}/intermediates/${RG_ID}_cleaned.bam ]; do
		sleep 180
	done
	
	if [ -e ${CLOUD_DIR}/${RG_ID}_cleaned.bam ]; then
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
