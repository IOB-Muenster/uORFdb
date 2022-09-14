#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node 16
#SBATCH --mem=48G
#SBATCH --time=10:00:00
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
# GATK GetPileUpSummaries for Variant Calling
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
CLOUD_DIR="/path/to/cloud/dir/"

WORKING_DIR=${SCRATCH_DIR}/Variant_Calling/${PATIENT}/
mkdir ${WORKING_DIR}

QC_BAM=${CLOUD_DIR}/Alignment_QC/${TYPE}/${UUID}/${UUID}_qc.bam

PILEUPS=${UUID}_pileups.table

EXAC=${SCRATCH_DIR}/Reference_genome/exac/small_exac_common_3_renamed.vcf.gz

echo -e "${TYPE}\n${UUID}\n${PATIENT}\n${TISSUE}\n${ID}\n"$(date "+%d/%m/%Y %H:%M:%S")
echo "--------------------------------------------------------"

PROCESS="12_getpileups"
START_TIME=$(date "+%s")

# ---------------------------------------------------------------------------------
# CHECK IF ALREADY PROCESSED
# ---------------------------------------------------------------------------------
if [ -e ${CLOUD_DIR}/Variant_Calling/${TYPE}/${PATIENT}/${PILEUPS} ]; then
	echo $(date "+%d/%m/%Y %H:%M:%S") "File exists already in Cloud. ${PROCESS} canceled."
	exit 0

elif [ -e ${WORKING_DIR}/${PILEUPS} ] || [ -e ${WORKING_DIR}/intermediates/${PILEUPS} ]; then
	echo $(date "+%d/%m/%Y %H:%M:%S") "File exists already in Scratch. ${PROCESS} canceled. Check for errors in previous processes."
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
# GETPILEUPSUMMARIES
# ---------------------------------------------------------------------------------
ml palma/2021a GCCcore/10.3.0 GATK/4.2.3.0-Java-1.8

gatk --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' GetPileupSummaries \
--tmp-dir ${SCRATCH_DIR}/Variant_Calling/tmp_dir/ \
-I ${QC_BAM} \
-V ${EXAC} \
-L ${EXAC} \
-O ${WORKING_DIR}/${PILEUPS} || \
	{ echo "GetPileSummaries failed"; exit 1; }

# ---------------------------------------------------------------------------------
# CHECK IF PROCESSED
# ---------------------------------------------------------------------------------
if [ -e ${WORKING_DIR}/${PILEUPS} ]; then
	
	echo $(date "+%d/%m/%Y %H:%M:%S") "Successfully created output."
	
	mkdir ${WORKING_DIR}/intermediates/
	mv ${WORKING_DIR}/${PILEUPS} ${WORKING_DIR}/intermediates/
	echo "${TYPE}" > ${SCRATCH_DIR}/vcalling_logs/intermediates_info/${PATIENT}

else
	echo $(date "+%d/%m/%Y %H:%M:%S") "Error. File not in Scratch."
	exit 1
fi	

# ---------------------------------------------------------------------------------
# CHECK RSYNC
# ---------------------------------------------------------------------------------
echo $(date "+%d/%m/%Y %H:%M:%S") "Checking if uploaded to Cloud."
	
while [ -e ${WORKING_DIR}/intermediates/${PILEUPS} ]; do
	sleep 180
done

if [ -e ${CLOUD_DIR}/Variant_Calling/${TYPE}/${PATIENT}/${PILEUPS} ]; then
	echo $(date "+%d/%m/%Y %H:%M:%S") "Successfully uploaded to Cloud."
fi

# ---------------------------------------------------------------------------------
# TIME LOG
# ---------------------------------------------------------------------------------
TIME_LOG=${SCRATCH_DIR}/vcalling_logs/runtimes/${UUID}_runtimes.log
END_TIME=$(date "+%s")
RUN_TIME=$((END_TIME-START_TIME))
echo "${PROCESS},${RUN_TIME}" >> ${TIME_LOG}

echo $(date "+%d/%m/%Y %H:%M:%S") "End of script."
