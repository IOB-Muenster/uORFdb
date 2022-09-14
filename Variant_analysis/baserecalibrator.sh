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
# GATK baserecalibrator for Variant Calling
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
CLOUD_DIR="/path/to/cloud/dir/${TYPE}/${UUID}/"

DOWNLOADED_BAM=$(ls /path/to/downloaded/BAMs/${TYPE}/${UUID}/*.bam)

cd ${SCRATCH_DIR}

mkdir ${SCRATCH_DIR}/Alignment_QC/${UUID}
WORKING_DIR=${SCRATCH_DIR}/Alignment_QC/${UUID}/
cd ${WORKING_DIR}

REF=/path/to/ref/genome/GCF_000001405.39_GRCh38.p13_genomic.fna
KNOWN=/path/to/known/vars/GCF_000001405.39.gz
MARKED_BAM=${CLOUD_DIR}/${UUID}_marked.bam

BQSR=${WORKING_DIR}/${UUID}.bqsr.grp

PROCESS="07_baserecalibrator"
START_TIME=$(date "+%s")

echo $TYPE
echo $UUID
echo $PATIENT $TISSUE
echo $(date "+%d/%m/%Y %H:%M:%S")
echo "--------------------------------------------------------"

# ---------------------------------------------------------------------------------
# INDEX AND DICT REFERENCE GENOME
# ---------------------------------------------------------------------------------
if [ ! -e ${REF}.fai ]; then
	echo $(date "+%d/%m/%Y %H:%M:%S") "Reference genome not indexed yet. Samtools faidx started."
	ml palma/2019a GCC/8.2.0-2.31.1 SAMtools/1.9 
	samtools faidx ${REF}
	echo $(date "+%d/%m/%Y %H:%M:%S") "Samtools faidx finished."
else
	echo $(date "+%d/%m/%Y %H:%M:%S") "Reference genome already indexed."
fi

if [ ! -e ${SCRATCH_DIR}/Reference_genome/ncbi-genomes-2021-07-15/GCF_000001405.39_GRCh38.p13_genomic.dict ]; then
	echo $(date "+%d/%m/%Y %H:%M:%S") "Reference genome not indexed yet. Samtools dict started."
	ml palma/2019a GCC/8.2.0-2.31.1 SAMtools/1.9 
	samtools dict ${REF} -o ${SCRATCH_DIR}/Reference_genome/ncbi-genomes-2021-07-15/GCF_000001405.39_GRCh38.p13_genomic.dict
	echo $(date "+%d/%m/%Y %H:%M:%S") "Samtools dict finished."
else
	echo $(date "+%d/%m/%Y %H:%M:%S") "Reference genome dict already exists."
fi

# ---------------------------------------------------------------------------------
# CHECK IF BQSR ALREADY EXISTS
# ---------------------------------------------------------------------------------
if [ -e ${CLOUD_DIR}/${UUID}.bqsr.grp ]; then
	echo $(date "+%d/%m/%Y %H:%M:%S") "File exists already in Cloud. ${PROCESS} canceled."
	exit 0
	
elif [ -e ${CLOUD_DIR}/${UUID}_qc.bam ]; then
	echo $(date "+%d/%m/%Y %H:%M:%S") "Alignment Pipeline seems to be finished already. QC file exists already in Cloud. ${PROCESS} canceled."
	exit 1

elif [ -e ${BQSR} ] || [ -e ${WORKING_DIR}/intermediates/${UUID}.bqsr.grp ]; then
	echo $(date "+%d/%m/%Y %H:%M:%S") "File exists already in Scratch. ${PROCESS} canceled. Check for errors in previous processes."
	echo $(date "+%d/%m/%Y %H:%M:%S")",${PROCESS},${UUID},ERROR_already_in_scratch" >> ${SCRATCH_DIR}/bwa_logs/bwa_error_messages.csv
	exit 1
	
else
	echo $(date "+%d/%m/%Y %H:%M:%S") "File does not exist yet."
fi

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
# BASERECALIBRATOR
# ---------------------------------------------------------------------------------
ml palma/2019a GCCcore/8.2.0 GATK/4.1.4.1-Java-11

gatk BaseRecalibrator \
-RF MappingQualityNotZeroReadFilter \
-R ${REF} \
-I ${MARKED_BAM} \
--known-sites ${KNOWN} \
-O ${BQSR} || \
    { echo "BaseRecalibrator failed"; exit 1; }
 
# ---------------------------------------------------------------------------------
# CHECK IF SUCCESSFUL
# ---------------------------------------------------------------------------------
if [ -e ${BQSR} ]; then
	echo $(date "+%d/%m/%Y %H:%M:%S") "BQSR was created."
	
	#enable rsync and deletion of input files
	mkdir ${WORKING_DIR}/intermediates
	mv ${BQSR} ${WORKING_DIR}/intermediates/
	echo "${TYPE}" > ${SCRATCH_DIR}/bwa_logs/intermediates_info/${UUID}
	
else
	echo $(date "+%d/%m/%Y %H:%M:%S") "${RG_ID}: Error. File not in Scratch."
	exit 1
fi

# ---------------------------------------------------------------------------------
# CHECK RSYNC
# ---------------------------------------------------------------------------------
echo $(date "+%d/%m/%Y %H:%M:%S") "Checking if uploaded to Cloud."

while [ -e ${WORKING_DIR}/intermediates/${UUID}.bqsr.grp ]; do
	sleep 180
done

if [ -e ${CLOUD_DIR}/${UUID}.bqsr.grp ]; then
	echo $(date "+%d/%m/%Y %H:%M:%S") "Successfully uploaded to Cloud."
else
	echo $(date "+%d/%m/%Y %H:%M:%S") "Missing. Check for errors."
	echo $(date "+%d/%m/%Y %H:%M:%S")",${PROCESS},${UUID},ERROR_missing" >> ${SCRATCH_DIR}/bwa_logs/bwa_error_messages.csv
	exit 1		
fi
    
# ---------------------------------------------------------------------------------
# TIME LOG
# ---------------------------------------------------------------------------------
END_TIME=$(date "+%s")
RUN_TIME=$((END_TIME-START_TIME))
echo "${PROCESS},${RUN_TIME}" >> ${SCRATCH_DIR}/bwa_logs/runtimes/${UUID}_runtimes.log

echo "End of script."












