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
# Filter VCFs to PASS variants and uORF regions and annotate subsequently
# to ClinVar and dbSNP databases
# Additional annotation to positions of ClinVar and dbSNP polymorphisms
#
# ==========================================================================================


while getopts t:r:c:d:e:f: flag
do
    case "${flag}" in
        t) TYPE=${OPTARG};;
		r) REGIONS=${OPTARG};;
		c) CLINVAR_POS=${OPTARG};;
		d) DBSNP_POS=${OPTARG};;
		e) CLINVAR_EXACT=${OPTARG};;
		f) DBSNP_EXACT=${OPTARG}
    esac
done

echo "${TYPE}"

ml palma/2020b GCC/10.2.0 BCFtools/1.11

CLOUD_DIR="/path/to/cloud/dir/${TYPE}/"
WORKING_DIR="/path/to/working/dir/${TYPE}/"
LOG=${WORKING_DIR}/${TYPE}_filter.log

mkdir ${WORKING_DIR}
mkdir ${WORKING_DIR}/original
mkdir ${WORKING_DIR}/filtered
mkdir ${WORKING_DIR}/annotated

echo -e "PATIENT\tVCF\tALL_VAR\tALL_VAR_PASS\tUORF_EXON_PASS_VAR" > ${LOG}

cd ${CLOUD_DIR}/
PATIENTS=$(ls )

for PATIENT in ${PATIENTS}; do

		echo ${PATIENT}

        cd ${CLOUD_DIR}/${PATIENT}/
        VCFS=$(ls *normalized*.vcf)

        for VCF in ${VCFS}; do
        
        		echo ${VCF}
        
                NEW_VCF=$(echo ${VCF} | sed 's/normalized_variants.vcf/filtered_variants.vcf.gz/')
                VCF_ID=$(echo ${NEW_VCF} | sed 's/_filtered_variants.vcf.gz//') 

				if [ ! -e ${WORKING_DIR}/annotated/${VCF_ID}_annotated.vcf ]; then

					# ---------------------------------------------------------------------------------
					# FILTERING
					# ---------------------------------------------------------------------------------
	        		echo "Filtering to uORFexon PASS variants ..."
	        
	                if [ ! -e ${WORKING_DIR}/${VCF}.gz ]; then
	        
	                        cp ${CLOUD_DIR}/${PATIENT}/${VCF} ${WORKING_DIR}/original/
	                        cd ${WORKING_DIR}/original
	                        bgzip -f ${WORKING_DIR}/original/${VCF}
	                        tabix -p vcf ${WORKING_DIR}/original/${VCF}.gz
	                        
	                        
	                        bcftools filter \
	                        -R ${REGIONS} \
	                        -i "%FILTER='PASS'" \
	                        -O z \
	                        -o ${WORKING_DIR}/filtered/${NEW_VCF} \
	                        ${WORKING_DIR}/original/${VCF}.gz
	                        
	                        tabix -p vcf ${WORKING_DIR}/filtered/${NEW_VCF}
	                
	                fi
	                
	                ALL_VAR=$(less ${CLOUD_DIR}/${PATIENT}/${VCF} | grep -v "#" | wc -l)
	                ALL_VAR_PASS=$(less ${CLOUD_DIR}/${PATIENT}/${VCF} | grep -v "#" | grep "PASS" | wc -l)
					UORF_EXON_PASS_VAR=$(less ${WORKING_DIR}/filtered/${NEW_VCF} | grep -v "#" | wc -l)
					
	                echo -e "${PATIENT}\t${VCF}\t${ALL_VAR}\t${ALL_VAR_PASS}\t${UORF_EXON_PASS_VAR}\t" >> ${LOG}
	                echo -e "${PATIENT}\t${VCF}\t${ALL_VAR}\t${ALL_VAR_PASS}\t${UORF_EXON_PASS_VAR}\t"
	                
	                # ---------------------------------------------------------------------------------
					# ANNOTATION
					# ---------------------------------------------------------------------------------
	                echo ${NEW_VCF}
	
					VCF_INPUT="${WORKING_DIR}/filtered/${NEW_VCF}"
					VCF_OUT1="${WORKING_DIR}/annotated/${VCF_ID}_annotated_intermediate1.vcf.gz"
					VCF_OUT2="${WORKING_DIR}/annotated/${VCF_ID}_annotated_intermediate2.vcf.gz"
					VCF_OUT3="${WORKING_DIR}/annotated/${VCF_ID}_annotated_intermediate3.vcf.gz"
					VCF_OUT4="${WORKING_DIR}/annotated/${VCF_ID}_annotated.vcf.gz"
				
					# ------------------------------------
					# ANNOTATE POSITIONS WITH BCF TOOLS
					# ------------------------------------
					echo "Annotating ClinVar POS ..."
					
					bcftools annotate \
					-c CHROM,FROM,TO,=ID \
					-a ${CLINVAR_POS} \
					-o ${VCF_OUT1} \
					-O z \
					${VCF_INPUT}
					
					tabix -p vcf ${VCF_OUT1}
					
					echo "Annotating dbSNP POS ..."
					
					bcftools annotate \
					-c CHROM,FROM,TO,=ID \
					-a ${DBSNP_POS} \
					-o ${VCF_OUT2} \
					-O z \
					${VCF_OUT1}
					
					tabix -p vcf ${VCF_OUT2}
					rm ${VCF_OUT1}
					rm ${VCF_OUT1}.tbi
					
					# ------------------------------------
					# ANNOTATE EXACT WITH BCF TOOLS
					# ------------------------------------
					echo "Annotating ClinVar EXACT ..."
					
					bcftools annotate \
					-c CHROM,POS,=ID \
					-a ${CLINVAR_EXACT} \
					-o ${VCF_OUT3} \
					-O z \
					${VCF_OUT2}
					
					tabix -p vcf ${VCF_OUT3}
					rm ${VCF_OUT2}
					rm ${VCF_OUT2}.tbi
				
					echo "Annotating dbSNP EXACT ..."
					
					bcftools annotate \
					-c CHROM,POS,=ID,INFO \
					-a ${DBSNP_EXACT} \
					-o ${VCF_OUT4} \
					-O z \
					${VCF_OUT3}
					
					tabix -p vcf ${VCF_OUT4}
					gunzip ${VCF_OUT4}
					rm ${VCF_OUT3}
					rm ${VCF_OUT3}.tbi
					
					echo "Finished."
					
				else
				
					echo "Already annotated."
					ALL_VAR=$(less ${CLOUD_DIR}/${PATIENT}/${VCF} | grep -v "#" | wc -l)
	                ALL_VAR_PASS=$(less ${CLOUD_DIR}/${PATIENT}/${VCF} | grep -v "#" | grep "PASS" | wc -l)
					UORF_EXON_PASS_VAR=$(less ${WORKING_DIR}/filtered/${NEW_VCF} | grep -v "#" | wc -l)
					echo -e "${PATIENT}\t${VCF}\t${ALL_VAR}\t${ALL_VAR_PASS}\t${UORF_EXON_PASS_VAR}\t" >> ${LOG}
	                echo -e "${PATIENT}\t${VCF}\t${ALL_VAR}\t${ALL_VAR_PASS}\t${UORF_EXON_PASS_VAR}\t"
					
				fi
                 
        done

done

