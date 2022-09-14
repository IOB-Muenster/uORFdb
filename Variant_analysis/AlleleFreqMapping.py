#!/usr/bin/env python3

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
# Mapping of allele frequencies to variants, merging identical variants to one entry
# (patient data is merged to lists separated by ";"
#
# ==========================================================================================

# ---------------------------------------------------------------------------------
# PACKAGES
# ---------------------------------------------------------------------------------
import pandas as pd
import argparse
from pandarallel import pandarallel
from map_AF import *
from natsort import natsorted
import numpy as np
from tqdm import tqdm

# ---------------------------------------------------------------------------------
# DEFINE INPUT / OUTPUT
# ---------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("-d","--directory",help="working dir")
parser.add_argument("--dbsnp",help="directory that contains dbSNP TSVs with allele counts")
args = parser.parse_args()

directory = args.directory
dbSNP_dir = args.dbsnp

variants_file = directory+"VariantsWithCoverage.tsv"

# ---------------------------------------------------------------------------------
# INFO
# ---------------------------------------------------------------------------------
print("+-----------------------------------+")
print("|       dbSNP AF Mapping 1.0        |")
print("+-----------------------------------+\n")
print("working directory: "+directory)
print("dbSNP_data: "+dbSNP_dir)
pandarallel.initialize(progress_bar=True, nb_workers=8)#, use_memory_fs=False
print("")

# ---------------------------------------------------------------------------------
# LOADING VARIANTS
# ---------------------------------------------------------------------------------
print("Loading variants data ... ")
selected_variants = pd.read_csv(variants_file, delimiter="\t")

# ---------------------------------------------------------------------------------
# PROCESSING EACH CHR: ADDING DBSNP DATA AND MERGING
# ---------------------------------------------------------------------------------
print("Adding dbSNP counts and merging variants ... ")

chromosomes = natsorted(selected_variants["CHROM"].unique())

for chromosome in chromosomes:
    
    variant_subset = selected_variants[(selected_variants.CHROM == chromosome)]
    print("")
    print(chromosome+"\tNumber of variants: "+str(len(variant_subset)))

    # ---------------------------------------------------------------------------------
    # ADDING dbSNP NUMBERS FROM JSON
    # ---------------------------------------------------------------------------------
    contig_alias = {"NC_000001.11" : "chr1","NC_000002.12" : "chr2","NC_000003.12" : "chr3","NC_000004.12" : "chr4","NC_000005.10" : "chr5",
                    "NC_000006.12" : "chr6","NC_000007.14" : "chr7","NC_000008.11" : "chr8","NC_000009.12" : "chr9","NC_000010.11" : "chr10",
                    "NC_000011.10" : "chr11","NC_000012.12" : "chr12","NC_000013.11" : "chr13","NC_000014.9" : "chr14","NC_000015.10" : "chr15",
                    "NC_000016.10" : "chr16","NC_000017.11" : "chr17","NC_000018.10" : "chr18","NC_000019.10" : "chr19","NC_000020.11" : "chr20",
                    "NC_000021.9" : "chr21","NC_000022.11" : "chr22","NC_000023.11" : "chrX","NC_000024.10" : "chrY"}
    
    print("    Loading dbSNP data ... ")
    dbSNP_data = dbSNP_dir+contig_alias[chromosome]+"_uORFexonv2_AF.tsv"
    dbSNP_data = pd.read_csv(dbSNP_data, delimiter = "\t")
    dbSNP_data.fillna("N", inplace=True)
    dbSNP_data = dbSNP_data.to_dict("index")
    dbSNP_dict = {}
    for entry in dbSNP_data:
        identifier = str(dbSNP_data[entry]["ID"])+"_"+dbSNP_data[entry]["REF"]+"_"+dbSNP_data[entry]["ALT"]
        dbSNP_dict[identifier] = dbSNP_data[entry]
    
    print("    Adding dbSNP data ... ")    
    variant_subset = variant_subset.parallel_apply(add_dbSNP_numbers,dbSNP_dict=dbSNP_dict,axis=1)
    print("")

    # ---------------------------------------------------------------------------------
    # MERGING VARIANTS
    # ---------------------------------------------------------------------------------
    identifier_columns = ["CHROM","POS","ALT","REF"]

    merge_columns = ["N_GT","T_GT","N_REF_AlleleDepth","T_REF_AlleleDepth",
                     "N_ALT_AlleleDepth","T_ALT_AlleleDepth","N_ReadDepth","T_ReadDepth",
                     "N_SomaticMutationFrequency","T_SomaticMutationFrequency","Proportion_T_SMF/N_SMF","CUSTOM_FILTER",
                     "BAM_ReadDepth_N", "BAM_ReadDepth_T",
                     "cancer_type","patient","T_sample_type","N_UUID","T_UUID","gender"]
    
    sum_columns = ["TCGA_ALT_AlleleCount","TCGA-BRCA_ALT_AlleleCount", "TCGA-COAD_ALT_AlleleCount", "TCGA-LAML_ALT_AlleleCount",
                   "TCGA-LUAD_ALT_AlleleCount", "TCGA-PRAD_ALT_AlleleCount", "TCGA-SKCM_ALT_AlleleCount"]
    
    cancer_types = ["TCGA-BRCA","TCGA-COAD","TCGA-LAML","TCGA-LUAD","TCGA-PRAD","TCGA-SKCM"]
    
    print("    Merging variants and adding cancer type counts ... ") 
    merged_variants = merge_variants(variant_dataframe=variant_subset, identifier_columns=identifier_columns, merge_columns=merge_columns, sum_columns=sum_columns, cancer_types=cancer_types)
    merged_variants = merged_variants.parallel_apply(convert_list_to_text, dataframe=merged_variants, sep=";", axis = 1)
    print("")
    print("    Merged: "+str(len(variant_subset))+" -> "+str(len(merged_variants)))
    
    print("    Adding variant allele frequencies ...")
    merged_variants = merged_variants.parallel_apply(addCancerAF, axis = 1)
    print("")
    
    # ---------------------------------------------------------------------------------
    # OUTPUT
    # ---------------------------------------------------------------------------------
    sorted_columns = [  # ID columns
                        "CHROM","POS","FROM_0BHO","TO_0BHO","REF","ALT",
                        # same values
                        "dbSNP_ID_POS","ClinVar_ID_POS","dbSNP_ID_EXACT","ClinVar_ID_EXACT",
                        "ExAC_fromVCF","GnomAD_fromVCF","TOPMED_fromVCF","ID.index",
                        "ExAC.1_ref/alt/total","GnomAD.3_ref/alt/total","TOPMED.2_ref/alt/total",
                        "TOPMED.3_ref/alt/total","match_info","match_info2",
                        # merge columns -> lists
                        "N_GT","T_GT","N_REF_AlleleDepth","T_REF_AlleleDepth",
                        "N_ALT_AlleleDepth","T_ALT_AlleleDepth","N_ReadDepth","T_ReadDepth",
                        "N_SomaticMutationFrequency","T_SomaticMutationFrequency","Proportion_T_SMF/N_SMF","CUSTOM_FILTER",
                        "BAM_ReadDepth_N", "BAM_ReadDepth_T","VCF_PatCount",
                        # same values
                        "ReadDepth_mean_N/T", "ReadDepth_median_N/T", "ReadDepth_std_N/T",
                        # add columns -> sum
                        "VCF_PatCount",
                        "TCGA_ALT_AlleleCount","TCGA-BRCA_ALT_AlleleCount", "TCGA-COAD_ALT_AlleleCount", "TCGA-LAML_ALT_AlleleCount",
                        "TCGA-LUAD_ALT_AlleleCount", "TCGA-PRAD_ALT_AlleleCount", "TCGA-SKCM_ALT_AlleleCount",
                        # same values
                        "TCGA_AlleleCount","TCGA-BRCA_AlleleCount", "TCGA-COAD_AlleleCount", "TCGA-LAML_AlleleCount",
                        "TCGA-LUAD_AlleleCount", "TCGA-PRAD_AlleleCount", "TCGA-SKCM_AlleleCount",
                        "TCGA_AF", "TCGA-BRCA_AF", "TCGA-COAD_AF", "TCGA-LAML_AF", "TCGA-LUAD_AF", "TCGA-PRAD_AF", "TCGA-SKCM_AF",
                        # merge columns -> lists
                        "cancer_type","patient","T_sample_type","N_UUID","T_UUID","gender"]
    merged_variants = merged_variants.reindex(columns=sorted_columns)

    print("    Writing to output ... ")
    merged_variants.to_csv(directory+"/"+chromosome+"_VariantsWithAF.tsv", sep = "\t", index = False)

print("Done.")

