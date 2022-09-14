#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node 18
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
# Conversion of BAM to fastq files with biobambam2
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

TARGET_CLOUD_DIR=${CLOUD_DIR}/Alignment_QC/${TYPE}/${UUID}/

DOWNLOADED_BAM=$(ls /path/to/downloaded/BAMs/${TYPE}/${UUID}/*.bam))

cd ${SCRATCH_DIR}
mkdir ${SCRATCH_DIR}/Alignment_QC/${UUID}
WORKING_DIR=${SCRATCH_DIR}/Alignment_QC/${UUID}/

cd ${WORKING_DIR}

SIZE=$(wc -c ${DOWNLOADED_BAM} | awk '{print $1}')
PROCESS="01_bam2fastq"
START_TIME=$(date "+%s")

echo $TYPE
echo $UUID
echo $PATIENT $TISSUE
echo $(date "+%d/%m/%Y %H:%M:%S")
echo "--------------------------------------------------------"

# ---------------------------------------------------------------------------------
# CHECK IF QC ALREADY EXISTS
# ---------------------------------------------------------------------------------
if [ -e ${TARGET_CLOUD_DIR}/${UUID}_qc.bam ]; then
	echo $(date "+%d/%m/%Y %H:%M:%S") "Alignment Pipeline seems to be finished already. All processes are canceled."
	exit 1
fi

# ---------------------------------------------------------------------------------
# CHECK IF ALREADY CONVERTED FOR EACH READ GROUP
# ---------------------------------------------------------------------------------
ml palma/2019a GCC/8.2.0-2.31.1 SAMtools/1.9
RG_IDS=$(samtools view -H ${DOWNLOADED_BAM} | grep "^@RG" | awk '{print $2}' | sed 's/ID://1')

for RG_ID in ${RG_IDS}
do
	if [ -e ${TARGET_CLOUD_DIR}/${RG_ID}_1.fq.gz ] || [ -e ${TARGET_CLOUD_DIR}/${RG_ID}_s.fq.gz ] || [ -e ${TARGET_CLOUD_DIR}/${RG_ID}_o1.fq.gz ]; then
		echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}: File exists already in Cloud. ${PROCESS} canceled."
		exit 0
		
	elif [ -e ${TARGET_CLOUD_DIR}/${UUID}_qc.bam ]; then
		echo $(date "+%d/%m/%Y %H:%M:%S") "Alignment Pipeline seems to be finished already. QC file exists already in Cloud. ${PROCESS} canceled."
		exit 1
	
	elif [ -e ${WORKING_DIR}/${RG_ID}_1.fq.gz ] || [ -e ${WORKING_DIR}/intermediates/${RG_ID}_1.fq.gz ]; then
		echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}: File exists already in Scratch. ${PROCESS} canceled. Check for errors in previous processes."
		echo $(date "+%d/%m/%Y %H:%M:%S")",${PROCESS},${UUID},${RG_ID},ERROR_already_in_scratch" >> ${SCRATCH_DIR}/bwa_logs/bwa_error_messages.csv
		
	else
		echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}: File does not exist yet."
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
# CONVERT TO FASTQ (AUTOMATICALLY FOR EACH READ GROUP)
# ---------------------------------------------------------------------------------
ml palma/2020b GCCcore/10.2.0 biobambam2/2.0.180

# important! cd to to the directory where intermediate file should be stored!
cd ${WORKING_DIR}

bamtofastq \
collate=1 \
exclude=QCFAIL,SECONDARY,SUPPLEMENTARY \
filename=${DOWNLOADED_BAM} \
gz=1 \
inputformat=bam \
level=5 \
outputdir=${WORKING_DIR} \
outputperreadgroup=1 \
outputperreadgroupsuffixF=_1.fq.gz \
outputperreadgroupsuffixF2=_2.fq.gz \
outputperreadgroupsuffixO=_o1.fq.gz \
outputperreadgroupsuffixO2=_o2.fq.gz \
outputperreadgroupsuffixS=_s.fq.gz \
tryoq=1

# ---------------------------------------------------------------------------------
# CHECK IF CONVERSION COMPLETED FOR EACH READ GROUP ID
# ---------------------------------------------------------------------------------
for RG_ID in ${RG_IDS}
do
	if [ -e ${WORKING_DIR}/${RG_ID}_1.fq.gz ] && [ -e ${WORKING_DIR}/${RG_ID}_2.fq.gz ]; then
		echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}: BAM file successfully converted to 1 and 2."
		
		mkdir ${WORKING_DIR}/intermediates
		mv ${WORKING_DIR}/${RG_ID}_*.fq.gz ${WORKING_DIR}/intermediates/
		echo "${TYPE}" > ${SCRATCH_DIR}/bwa_logs/intermediates_info/${UUID}
		
	elif [ -e ${WORKING_DIR}/${RG_ID}_o1.fq.gz ] && [ -e ${WORKING_DIR}/${RG_ID}_o2.fq.gz ]; then
		echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}: BAM file successfully converted to o1 and o2."
		
		mkdir ${WORKING_DIR}/intermediates
		mv ${WORKING_DIR}/${RG_ID}_*.fq.gz ${WORKING_DIR}/intermediates/
		echo "${TYPE}" > ${SCRATCH_DIR}/bwa_logs/intermediates_info/${UUID}
		
	elif [ -e ${WORKING_DIR}/${RG_ID}_s.fq.gz ]; then
		echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}: BAM file successfully converted to s."
		
		mkdir ${WORKING_DIR}/intermediates
		mv ${WORKING_DIR}/${RG_ID}_*.fq.gz ${WORKING_DIR}/intermediates/
		echo "${TYPE}" > ${SCRATCH_DIR}/bwa_logs/intermediates_info/${UUID}
		
	else
		echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}: Error. Fastq file not in Scratch."
		exit 1
		
	fi
done

# ---------------------------------------------------------------------------------
# CHECK RSYNC
# ---------------------------------------------------------------------------------
for RG_ID in ${RG_IDS}
do
	
	echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}: Checking if uploaded to Cloud."
	
	FILE=${WORKING_DIR}/intermediates/${RG_ID}
	CLOUD_FILE=${TARGET_CLOUD_DIR}/${RG_ID}
	
	while [ -e ${FILE}_1.fq.gz ] || [ -e ${FILE}_2.fq.gz ] || [ -e ${FILE}_o1.fq.gz ] || [ -e ${FILE}_o2.fq.gz ] || [ -e ${FILE}_s.fq.gz ]; do
		sleep 180
	done
	
	if [ -e ${CLOUD_FILE}_1.fq.gz ] || [ -e ${CLOUD_FILE}_2.fq.gz ] || [ -e ${CLOUD_FILE}_o1.fq.gz ] || [ -e ${CLOUD_FILE}_o2.fq.gz ] || [ -e ${CLOUD_FILE}_s.fq.gz ]; then
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
echo "TYPE,${TYPE}" >> ${SCRATCH_DIR}/bwa_logs/runtimes/${UUID}_runtimes.log
echo "UUID,${UUID}" >> ${SCRATCH_DIR}/bwa_logs/runtimes/${UUID}_runtimes.log
echo "PATIENT,${PATIENT}" >> ${SCRATCH_DIR}/bwa_logs/runtimes/${UUID}_runtimes.log
echo "TISSUE,${TISSUE}" >> ${SCRATCH_DIR}/bwa_logs/runtimes/${UUID}_runtimes.log
echo "SIZE,${SIZE}" >> ${SCRATCH_DIR}/bwa_logs/runtimes/${UUID}_runtimes.log
echo "${PROCESS},${RUN_TIME}" >> ${SCRATCH_DIR}/bwa_logs/runtimes/${UUID}_runtimes.log

echo "End of script."
