#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node 16
#SBATCH --mem=16GB
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
# Picard CrossCheckFingerprint for Variant Calling
#
# tisse is always tumor, normal UUID (reference UUID) info is given under flag "n"
# iterate over tumor files and make a pipeline for each tumor file, normal files are involved in this step
# patient info is always given in output files as well, so that matching is easier later
#
# https://gatk.broadinstitute.org/hc/en-us/articles/360037057832-CrosscheckFingerprints-Picard-
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

T_ID=$(echo ${T_UUID} | cut -d'-' -f 1)
N_ID=$(echo ${N_UUID} | cut -d'-' -f 1)
ID="${PATIENT}_${T_ID}_${N_ID}"

SCRATCH_DIR="/path/to/working/dir/"
CLOUD_DIR="/path/to/cloud/dir/"

mkdir ${SCRATCH_DIR}/Variant_Calling/${PATIENT}/
WORKING_DIR=${SCRATCH_DIR}/Variant_Calling/${PATIENT}/

TUMOR_BAM=${CLOUD_DIR}/Alignment_QC/${TYPE}/${T_UUID}/${T_UUID}_qc.bam
NORMAL_BAM=${CLOUD_DIR}/Alignment_QC/${TYPE}/${N_UUID}/${N_UUID}_qc.bam

HAP=${SCRATCH_DIR}/Reference_genome/hapmap/hapmap_bamheader_LO.txt

OUT=${ID}.crosscheck_metrics

PICARD="/path/to/picard/picard.jar"

echo -e "${TYPE}\n${UUID}\n${PATIENT}\n${TISSUE}\n${ID}\n"$(date "+%d/%m/%Y %H:%M:%S")
echo "--------------------------------------------------------"

T_SIZE=$(wc -c ${TUMOR_BAM} | awk '{print $1}')
N_SIZE=$(wc -c ${NORMAL_BAM} | awk '{print $1}')
PROCESS="10_fingerprint"
START_TIME=$(date "+%s")

# ---------------------------------------------------------------------------------
# CHECK IF ALREADY PROCESSED
# ---------------------------------------------------------------------------------
if [ -e ${CLOUD_DIR}/Variant_Calling/${TYPE}/${PATIENT}/${OUT} ]; then

	echo "Exists already in Cloud. Process canceled."
	exit 0
	
elif [ -e ${WORKING_DIR}/${OUT} ] || [ -e ${WORKING_DIR}/intermediates/${OUT} ]; then

	echo "Exists already in Scratch. Process canceled. Please check previous processes."
	exit 1
	
else

	echo "Not processed yet. CrosscheckFingerprints starts."
	
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
	# CROSSCHECKFINGERPRINTS
	# ---------------------------------------------------------------------------------
	ml palma/2020b Java/11.0.2
	
	java -jar ${PICARD} CrosscheckFingerprints \
	-INPUT ${TUMOR_BAM} \
	-INPUT ${NORMAL_BAM} \
	-HAPLOTYPE_MAP ${HAP} \
	-LOD_THRESHOLD -5 \
	-EXPECT_ALL_GROUPS_TO_MATCH true \
	-OUTPUT ${WORKING_DIR}/${OUT} || \
		{ echo "CrosscheckFingerprints failed"; exit 1; }
	
	mkdir ${WORKING_DIR}/intermediates
	mv ${WORKING_DIR}/${OUT} ${WORKING_DIR}/intermediates/
	echo "${TYPE}" > ${SCRATCH_DIR}/vcalling_logs/intermediates_info/${PATIENT}

fi

# ---------------------------------------------------------------------------------
# CHECK RSYNC
# ---------------------------------------------------------------------------------
echo $(date "+%d/%m/%Y %H:%M:%S") "Checking if uploaded to Cloud."
	
while [ -e ${WORKING_DIR}/intermediates/${OUT} ]; do
	sleep 180
done

if [ -e ${CLOUD_DIR}/Variant_Calling/${TYPE}/${PATIENT}/${OUT} ]; then
	echo $(date "+%d/%m/%Y %H:%M:%S") "Successfully uploaded to Cloud."
fi

# ---------------------------------------------------------------------------------
# TIME LOG
# ---------------------------------------------------------------------------------
TIME_LOG=${SCRATCH_DIR}/vcalling_logs/runtimes/${ID}_runtimes.log
END_TIME=$(date "+%s")
RUN_TIME=$((END_TIME-START_TIME))
echo "TYPE,${TYPE}" >> ${TIME_LOG}
echo "PATIENT,${PATIENT}" >> ${TIME_LOG}
echo "T_UUID,${T_UUID}" >> ${TIME_LOG}
echo "N_UUID,${N_UUID}" >> ${TIME_LOG}
echo "T_SIZE,${SIZE}" >> ${TIME_LOG}
echo "N_SIZE,${N_SIZE}" >> ${TIME_LOG}
echo "${PROCESS},${RUN_TIME}" >> ${TIME_LOG}

echo $(date "+%d/%m/%Y %H:%M:%S") "End of script."
