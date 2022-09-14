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
# Reformat variants from multiple VCFs and assign metadata
#
# ==========================================================================================

# ---------------------------------------------------------------------------------
# PACKAGES
# ---------------------------------------------------------------------------------
import pandas as pd
import argparse
from pandarallel import pandarallel
from preprocess_var import *
from natsort import natsorted
import numpy as np
from tqdm import tqdm
tqdm.pandas()

# ---------------------------------------------------------------------------------
# DEFINE INPUT / OUTPUT
# ---------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("-i","--input",help="input directory with vcf files")
parser.add_argument("-o","--output",help="output directory")
parser.add_argument("--meta",help="CSV with metadata")
parser.add_argument("--vcflist",help="list of vcfs to analyze")
parser.add_argument("--samplenames",help="tsv with sample names and UUIDs")
args = parser.parse_args()

vcf_list = args.vcflist
directory = args.input
output_dir = args.output
metadata_file = args.meta
samplenames = args.samplenames

selection_file = output_dir+"/selection.tsv"
variants_file = output_dir+"/variants.tsv"
bed_file = output_dir+"/VarCoverage_positions.bed"

# ---------------------------------------------------------------------------------
# INFO
# ---------------------------------------------------------------------------------
print("+-----------------------------------+")
print("|     VARIANT PREPROCESSING 1.0     |")
print("+-----------------------------------+\n")
print("vcflist: "+vcf_list)
print("input directory: "+directory)
print("output directory: "+output_dir)
print("TCGA metadata: "+metadata_file)
print("Sample names: "+samplenames)
pandarallel.initialize(progress_bar=True, nb_workers=8)
print("")

# ---------------------------------------------------------------------------------
# CONVERT INPUT
# ---------------------------------------------------------------------------------
print("Loading metadata ... ")
metadata = (pd.read_csv(metadata_file)).to_dict("index")
metadata_dict = {}
for i in metadata:
    identifier = metadata[i]["patient_id"]+"-"+metadata[i]["UUID"][0:8]
    metadata_dict[identifier] = metadata[i]

print("Loading VCF list ... ")
vcfs = []
with open (vcf_list, "r") as file:
    for line in file:
        vcfs.append(line[0:-1])
        
print("Loading sample names ...")
sample_dict = {}
with open (samplenames, "r") as file:
    for line in file:
        uuid = line[:-1].split("\t")[0]
        sample = line[:-1].split("\t")[1]
        sample_dict[uuid] = sample

# ---------------------------------------------------------------------------------
# SELECT SAMPLES
# ---------------------------------------------------------------------------------
print("Collecting VCFs ... ")
dataframes = []
for vcf in vcfs:
    dataframe = cut_vcfs(vcf=vcf, input_directory=directory, output_directory=output_dir, metadata_dict=metadata_dict, sample_dict=sample_dict)
    dataframes.append(dataframe)
vcf_df = pd.concat(dataframes)

print("Selecting VCFs ... ") 
vcf_df = vcf_df.progress_apply(select_samples, dataframe = vcf_df, axis = 1)

sorted_columns = ["cancer_type","patient","T_UUID","N_UUID","tumor_sample","normal_sample","T_sample_type","N_sample_type","file_ids","file","tumor_bam_size","gender",
                  "variant_count","variant_index","selection_status"] #added "gender" 17.08.22
vcf_df = vcf_df.reindex(columns=sorted_columns)
print("Writing selection log output ... ")
vcf_df.to_csv(selection_file, sep="\t", index = False)

# ---------------------------------------------------------------------------------
# ADDING VARIANTS DATA
# ---------------------------------------------------------------------------------
selected_vcfs = vcf_df[(vcf_df.selection_status == "")]
print("Selected "+str(len(selected_vcfs))+" out of "+str(len(vcf_df))+" VCF files.")

selected_vcfs = selected_vcfs.explode("variant_index").reset_index(drop=True)

vcf_files = selected_vcfs["file"].unique()

variants_dataframes = []
for file in vcf_files: 
    variants = pd.read_csv(file, sep = "\t")
    variants["file"] = file
    variants_dataframes.append(variants)
variants = pd.concat(variants_dataframes)

print("Adding variant data ... ")

selected_variants = selected_vcfs.parallel_apply(add_variant_data, variants_dataframe=variants, axis = 1)
print("")
print("Found "+str(len(selected_variants))+" variants.")

sorted_columns = ["CHROM","POS","FROM_0BHO","TO_0BHO","REF","ALT","dbSNP_ID_POS","ClinVar_ID_POS","dbSNP_ID_EXACT","ClinVar_ID_EXACT",
                  "ExAC_fromVCF","GnomAD_fromVCF","TOPMED_fromVCF","N_GT","T_GT",
                  "N_REF_AlleleDepth","T_REF_AlleleDepth","N_ALT_AlleleDepth","T_ALT_AlleleDepth","N_ReadDepth","T_ReadDepth",
                  "N_SomaticMutationFrequency","T_SomaticMutationFrequency","Proportion_T_SMF/N_SMF","CUSTOM_FILTER",
                  "cancer_type","patient","N_sample_type","T_sample_type","N_UUID","T_UUID","tumor_sample","normal_sample","file","gender"
                  ]

selected_variants = selected_variants.reindex(columns=sorted_columns)

print("Writing variants output ... ")
selected_variants.to_csv(variants_file, sep = "\t", index = False)

# ---------------------------------------------------------------------------------
# MAKE POS BED FILE FOR COVERAGE ANALYSIS
# ---------------------------------------------------------------------------------
position_data = selected_variants.reindex(columns=["CHROM","FROM_0BHO","TO_0BHO"])
position_data = position_data.sort_values(by=["CHROM","FROM_0BHO","TO_0BHO"])
position_data = position_data.drop_duplicates()
position_data.to_csv(bed_file, sep = "\t", index = False, header = False)     
        