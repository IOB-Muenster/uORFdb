#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node 24
#SBATCH --mem=48GB
#SBATCH --time=168:00:00
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
# GATK Mutect2 for Variant Calling
#
# tissue is always tumor, normal UUID (reference UUID) info is given under flag "n"
# iterate over tumor files and make a pipeline for each tumor file, normal files are involved in this step
# patient info is always given in output files as well, so that matching is easier later
#
# https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2
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

TISSUE="TN"

SCRATCH_DIR="/path/to/working/dir/"
CLOUD_DIR="/path/to/cloud/dir/"

mkdir ${SCRATCH_DIR}/Variant_Calling/${PATIENT}
WORKING_DIR=${SCRATCH_DIR}/Variant_Calling/${PATIENT}/

REF=${SCRATCH_DIR}/Reference_genome/ncbi-genomes-2021-07-15/GCF_000001405.39_GRCh38.p13_genomic.fna
PON=${SCRATCH_DIR}/Reference_genome/PON/pon_LO.vcf.gz
GNOMAD=${SCRATCH_DIR}/Reference_genome/gnomad/gnomad_LO.vcf.gz

DOWNLOADED_BAM=${CLOUD_DIR}/GDC_data/${TYPE}/${T_UUID}/*bam
TUMOR_BAM=${CLOUD_DIR}/Alignment_QC/${TYPE}/${T_UUID}/${T_UUID}_qc.bam
NORMAL_BAM=${CLOUD_DIR}/Alignment_QC/${TYPE}/${N_UUID}/${N_UUID}_qc.bam
FINGERPRINT=${CLOUD_DIR}/Variant_Calling/${TYPE}/${PATIENT}/${ID}.crosscheck_metrics

VARIANTS=${ID}_variants.vcf
F1R2=${ID}_f1r2.tar.gz

echo -e "${TYPE}\n${UUID}\n${PATIENT}\n${TISSUE}\n${ID}\n"$(date "+%d/%m/%Y %H:%M:%S")
echo "--------------------------------------------------------"

PROCESS="11_mutect2"
START_TIME=$(date "+%s")

# ---------------------------------------------------------------------------------
# CHECK FINGERPRINT
# ---------------------------------------------------------------------------------
MISMATCH=$(grep "UNEXPECTED_MISMATCH" ${FINGERPRINT})

if [ "${MISMATCH}" = "" ]; then
	echo $(date "+%d/%m/%Y %H:%M:%S") "Fingerprints match between tumor and normal file: (${FINGERPRINT})."	
else
	echo $(date "+%d/%m/%Y %H:%M:%S") "Fingerprints mismatch between tumor and normal file: (${FINGERPRINT}). Check input files!"
	echo $(date "+%d/%m/%Y %H:%M:%S")",${PROCESS},${T_UUID},${N_UUID},ERROR_mismatch" >> ${SCRATCH_DIR}/vcalling_logs/vcalling_error_messages.csv
	exit 1
fi

# ---------------------------------------------------------------------------------
# DETERMINE NORMAL SAMPLE NAME
# ---------------------------------------------------------------------------------
ml palma/2019a GCC/8.2.0-2.31.1 SAMtools/1.9

NORMAL_SAMPLE=$(samtools view -H ${NORMAL_BAM} | grep "@RG" | grep -ow '[S][M][:][\_.A-Z0-9-]*' | sed 's/SM://1' | head -n1 | cut -d " " -f1)

# ---------------------------------------------------------------------------------
# CHECK IF ALREADY PROCESSED
# ---------------------------------------------------------------------------------
if [ -e ${CLOUD_DIR}/Variant_Calling/${TYPE}/${PATIENT}/${VARIANTS} ]; then
	echo $(date "+%d/%m/%Y %H:%M:%S") "File exists already in Cloud. ${PROCESS} canceled."
	exit 0

elif [ -e ${WORKING_DIR}/${VARIANTS} ] || [ -e ${WORKING_DIR}/intermediates/${VARIANTS} ]; then
	echo $(date "+%d/%m/%Y %H:%M:%S") "File exists already in Scratch. ${PROCESS} canceled. Check for errors in previous processes."
	
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
# MUTECT2
# ---------------------------------------------------------------------------------
echo "starting mutect"
echo "-R ${REF}"
echo "-I ${TUMOR_BAM}"
echo "-I ${NORMAL_BAM}"
echo "-O ${WORKING_DIR}/${VARIANTS}"
echo "--normal-sample ${NORMAL_SAMPLE}"
echo "--f1r2-tar-gz ${WORKING_DIR}/${F1R2}"
echo "--germline-resource ${GNOMAD}"
echo "--panel-of-normals ${PON}"

ml palma/2021a GCCcore/10.3.0 GATK/4.2.3.0-Java-1.8

gatk --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' Mutect2 \
-I ${TUMOR_BAM} \
-I ${NORMAL_BAM} \
-O ${WORKING_DIR}/${VARIANTS} \
-R ${REF} \
--normal-sample ${NORMAL_SAMPLE} \
--f1r2-tar-gz ${WORKING_DIR}/${F1R2} \
--germline-resource ${GNOMAD} \
--panel-of-normals ${PON} || \
	{ echo "Mutect2 failed"; exit 1; }

# ---------------------------------------------------------------------------------
# CHECK IF PROCESS COMPLETED
# ---------------------------------------------------------------------------------
if [ -e ${WORKING_DIR}/${VARIANTS} ] && [ -e ${WORKING_DIR}/${VARIANTS}.stats ] && [ -e ${WORKING_DIR}/${VARIANTS}.idx ]; then
	echo $(date "+%d/%m/%Y %H:%M:%S") "Successfully created Mutect2 output files."
	
	mkdir ${WORKING_DIR}/intermediates
	mv ${WORKING_DIR}/${VARIANTS} ${WORKING_DIR}/intermediates/
	mv ${WORKING_DIR}/${VARIANTS}.stats ${WORKING_DIR}/intermediates/
	mv ${WORKING_DIR}/${VARIANTS}.idx ${WORKING_DIR}/intermediates/
	mv ${WORKING_DIR}/${F1R2} ${WORKING_DIR}/intermediates/
	echo "${TYPE}" > ${SCRATCH_DIR}/vcalling_logs/intermediates_info/${PATIENT}
	
else
	echo $(date "+%d/%m/%Y %H:%M:%S") "Error. Fastq file not in Scratch."
	exit 1
fi

# ---------------------------------------------------------------------------------
# CHECK RSYNC
# ---------------------------------------------------------------------------------
echo $(date "+%d/%m/%Y %H:%M:%S") "Checking if uploaded to Cloud."
	
while [ -e ${WORKING_DIR}/intermediates/${VARIANTS} ]; do
	sleep 180
done

if [ -e ${CLOUD_DIR}/Variant_Calling/${TYPE}/${PATIENT}/${VARIANTS} ]; then
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


