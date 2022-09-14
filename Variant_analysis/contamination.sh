#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node 16
#SBATCH --mem=48GB
#SBATCH --time=30:00:00
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
# GATK CalculateContamination for Variant Calling
#
# ==========================================================================================

# ---------------------------------------------------------------------------------
# VARIABLES
# ---------------------------------------------------------------------------------
while getopts i:t:p:n: flag
do
    case "${flag}" in
        i) T_UUID=${OPTARG};;
        t) TYPE=${OPTARG};;
        p) PATIENT=${OPTARG};;
        n) N_UUID=${OPTARG}
    esac
done

TISSUE="TN"

T_ID=$(echo ${T_UUID} | cut -d'-' -f 1)
N_ID=$(echo ${N_UUID} | cut -d'-' -f 1)
ID="${PATIENT}_${T_ID}_${N_ID}"

SCRATCH_DIR="/path/to/working/dir/"
CLOUD_DIR="/path/to/cloud/dir/"

WORKING_DIR=${SCRATCH_DIR}/Variant_Calling/${PATIENT}/
mkdir ${WORKING_DIR}

TUMOR_PILEUPS=${CLOUD_DIR}/Variant_Calling/${TYPE}/${PATIENT}/${T_UUID}_pileups.table
NORMAL_PILEUPS=${CLOUD_DIR}/Variant_Calling/${TYPE}/${PATIENT}/${N_UUID}_pileups.table

CONTAMINATION=${ID}_contamination.table

TUMOR_SEGMENTATION=${ID}_segments.table

echo -e "${TYPE}\n${UUID}\n${PATIENT}\n${TISSUE}\n${ID}\n"$(date "+%d/%m/%Y %H:%M:%S")
echo "--------------------------------------------------------"

PROCESS="13_contamination"
START_TIME=$(date "+%s")

# ---------------------------------------------------------------------------------
# CHECK IF ALREADY PROCESSED
# ---------------------------------------------------------------------------------
if [ -e ${CLOUD_DIR}/Variant_Calling/${TYPE}/${PATIENT}/${CONTAMINATION} ]; then
	echo $(date "+%d/%m/%Y %H:%M:%S") "File exists already in Cloud. ${PROCESS} canceled."
	exit 0

elif [ -e ${WORKING_DIR}/${CONTAMINATION} ] || [ -e ${WORKING_DIR}/intermediates/${CONTAMINATION} ]; then
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
# CALCULATECONTAMINATION
# ---------------------------------------------------------------------------------
ml palma/2021a GCCcore/10.3.0 GATK/4.2.3.0-Java-1.8

gatk --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' CalculateContamination \
-I ${TUMOR_PILEUPS} \
-matched ${NORMAL_PILEUPS} \
-O ${WORKING_DIR}/${CONTAMINATION} \
-tumor-segmentation ${WORKING_DIR}/${TUMOR_SEGMENTATION} || \
	{ echo "CalculateContamination failed"; exit 1; }

# ---------------------------------------------------------------------------------
# CHECK IF PROCESSED
# ---------------------------------------------------------------------------------
if [ -e ${WORKING_DIR}/${CONTAMINATION} ]; then
	
	echo $(date "+%d/%m/%Y %H:%M:%S") "Successfully calculated contamination."
	
	mkdir ${WORKING_DIR}/intermediates/
	mv ${WORKING_DIR}/*_pileups.table ${WORKING_DIR}/intermediates/
	mv ${WORKING_DIR}/${CONTAMINATION} ${WORKING_DIR}/intermediates/
	mv ${WORKING_DIR}/${TUMOR_SEGMENTATION} ${WORKING_DIR}/intermediates/
	echo "${TYPE}" > ${SCRATCH_DIR}/vcalling_logs/intermediates_info/${PATIENT}

else
	echo $(date "+%d/%m/%Y %H:%M:%S") "Error. Contamination File not in Scratch."
	exit 1
fi	

# ---------------------------------------------------------------------------------
# CHECK RSYNC
# ---------------------------------------------------------------------------------
echo $(date "+%d/%m/%Y %H:%M:%S") "Checking if uploaded to Cloud."
	
while [ -e ${WORKING_DIR}/intermediates/${CONTAMINATION} ]; do
	sleep 180
done

if [ -e ${CLOUD_DIR}/Variant_Calling/${TYPE}/${PATIENT}/${CONTAMINATION} ]; then
	echo $(date "+%d/%m/%Y %H:%M:%S") "Successfully uploaded to Cloud."
fi

# ---------------------------------------------------------------------------------
# TIME LOG
# ---------------------------------------------------------------------------------
TIME_LOG=${SCRATCH_DIR}/vcalling_logs/runtimes/${ID}_runtimes.log
END_TIME=$(date "+%s")
RUN_TIME=$((END_TIME-START_TIME))
echo "${PROCESS},${RUN_TIME}" >> ${TIME_LOG}

echo $(date "+%d/%m/%Y %H:%M:%S") "End of script."


