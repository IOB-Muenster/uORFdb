#!/bin/bash

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
# Extract sample names of downloaded BAMs
#
# ==========================================================================================
while getopts l: flag
do
    case "${flag}" in
        l) LOG=${OPTARG}
    esac
done

ml palma/2019a GCC/8.2.0-2.31.1 SAMtools/1.9

echo -e "UUID\tSAMPLE"

echo -e "UUID\tSAMPLE" > ${LOG}

for TYPE in "TCGA-BRCA" "TCGA-COAD" "TCGA-LAML" "TCGA-LUAD" "TCGA-PRAD" "TCGA-SKCM"; do

	echo ${TYPE}

	CLOUD_DIR="/path/to/cloud/dir/${TYPE}/"
	
	cd ${CLOUD_DIR}
	
	UUIDS=$(ls )
	
	for UUID in ${UUIDS}; do
	
		SAMPLE=$(samtools view -H ${CLOUD_DIR}/${UUID}/${UUID}_qc.bam | grep "@RG" | grep -ow '[S][M][:][\_.A-Z0-9-]*' | sed 's/SM://1' | head -n1 | cut -d " " -f1)
		
		echo -e "${UUID}\t${SAMPLE}"
		
		echo -e "${UUID}\t${SAMPLE}" >> ${LOG}
	
	done
	
done
