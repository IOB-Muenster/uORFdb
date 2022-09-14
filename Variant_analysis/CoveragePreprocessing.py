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
# Reformat read depth data for downstream analysis
#
# ==========================================================================================

# ---------------------------------------------------------------------------------
# PACKAGES
# ---------------------------------------------------------------------------------
import pandas as pd
from natsort import natsorted
import argparse
from progress.bar import ChargingBar

# ---------------------------------------------------------------------------------
# DEFINE INPUT / OUTPUT
# ---------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("-d","--dir", help="working directory")
parser.add_argument("-i","--cov_dir", help="directory with coverage data")
parser.add_argument("-c","--coverage", help="coverage list")
parser.add_argument("-v","--variants", help="files with variants")
args = parser.parse_args()

cov_dir = args.cov_dir
directory = args.dir
cov_list = args.coverage
variants = args.variants

# ---------------------------------------------------------------------------------
# INPUT
# ---------------------------------------------------------------------------------
variants = pd.read_csv(variants, delimiter = "\t")

files = []
file_count = 0 
with open (cov_list, "r") as file:
    for line in file:
        files.append(line[:-1])
        file_count += 1

metadata = variants.reindex(columns=["patient","N_UUID","T_UUID"]).drop_duplicates()
patient_dict = metadata.set_index("N_UUID").to_dict("index")

chromosomes = natsorted(variants["CHROM"].unique())

# ---------------------------------------------------------------------------------
# READING COVERAGE FROM FILES
# ---------------------------------------------------------------------------------
print("Reading coverage data from "+str(file_count)+" files ...")

coverage_data = []

with ChargingBar(max=file_count) as bar:

    for tissue in ["N","T"]:
    
        for file in files:
            if "_"+tissue+"_" in file:
            
                coverage = cov_dir+file
                coverage = pd.read_csv(coverage, delimiter = "\t")
                coverage = coverage.astype({"#CHROM":"str","POS":"str"})
                coverage["CHR_POS"] = coverage["#CHROM"]+"-"+coverage["POS"]
                coverage = coverage.set_index("CHR_POS")
                coverage.drop(["#CHROM","POS"], axis=1, inplace=True)
                coverage = coverage.transpose()
                coverage.rename(columns=coverage.iloc[0]).drop(coverage.index[0])
                coverage["tissue"] = tissue
                coverage_data.append(coverage)
                
            bar.next()

print("Merging ...")
merged_coverage = pd.concat(coverage_data)

# ---------------------------------------------------------------------------------
# ADDING STATS
# ---------------------------------------------------------------------------------
print("Calculating mean, median, stdev ...")

def add_stats(dataframe, tissue):
    
    stats = dataframe[dataframe.tissue == tissue]
    stats = stats.drop(columns=["tissue"]).transpose()
    stats["mean"] = dataframe[merged_coverage.tissue == tissue].mean(numeric_only=True)
    stats["median"] = dataframe[merged_coverage.tissue == tissue].median(numeric_only=True)
    stats["std"] = dataframe[merged_coverage.tissue == tissue].std(numeric_only=True)
    
    return stats

N_stats = add_stats(merged_coverage, tissue="N")
T_stats = add_stats(merged_coverage, tissue="T")

# ---------------------------------------------------------------------------------
# COLLECTING DATA FOR EACH CHROMOSOME, MATCHING TUMOR AND NORMAL
# ---------------------------------------------------------------------------------
coverage_meta_list = []

for chromosome in chromosomes:
    
    print(chromosome)
    print("Reformatting to dict ...")
    
    N_stats_dict = N_stats[N_stats.index.str.contains(chromosome)].to_dict("index")
    T_stats_dict = T_stats[T_stats.index.str.contains(chromosome)].to_dict("index")

    new_dict = {}

    for var in N_stats_dict:
        
        new_dict[var] = {}
        new_dict[var]["TCGA-BRCA_PatCount"] = 0
        new_dict[var]["TCGA-COAD_PatCount"] = 0
        new_dict[var]["TCGA-LAML_PatCount"] = 0
        new_dict[var]["TCGA-LUAD_PatCount"] = 0
        new_dict[var]["TCGA-PRAD_PatCount"] = 0
        new_dict[var]["TCGA-SKCM_PatCount"] = 0
        new_dict[var]["TCGA_PatCount"] = 0
        
        for file in N_stats_dict[var]:
            
            if "/" in file:
                TYPE = file.split("/")[-3]
                N_UUID = file.split("/")[-2]
                
                if N_UUID in patient_dict:
                    T_UUID = patient_dict[N_UUID]["T_UUID"]
                    patient = patient_dict[N_UUID]["patient"]
                    ID = patient+"_"+N_UUID+"_"+T_UUID
                    T_file = "/path/to/cloud/dir/"+TYPE+"/"+T_UUID+"/"+T_UUID+"_qc.bam"
                    
                    if file in N_stats_dict[var] and T_file in T_stats_dict[var]:  
                        N_COV = int(N_stats_dict[var][file])
                        T_COV = int(T_stats_dict[var][T_file])
                
                        new_dict[var][ID] = str(N_COV) + "," + str(T_COV)
                        
                        if N_COV >= 10 and T_COV >= 10:
                            new_dict[var][TYPE+"_PatCount"] += 1
                            new_dict[var]["TCGA_PatCount"] += 1

        new_dict[var]["mean_N/T"] = str(N_stats_dict[var]["mean"]) + "," + str(T_stats_dict[var]["mean"])
        new_dict[var]["median_N/T"] = str(N_stats_dict[var]["median"]) + "," + str(T_stats_dict[var]["median"])
        new_dict[var]["std_N/T"] = str(N_stats_dict[var]["std"]) + "," + str(T_stats_dict[var]["std"])
    
    print("Writing to output ...")
    
    coverage_data = pd.DataFrame.from_dict(new_dict, "index")

    meta_columns = ["TCGA_PatCount",
                    "TCGA-BRCA_PatCount", "TCGA-COAD_PatCount", "TCGA-LAML_PatCount",
                    "TCGA-LUAD_PatCount", "TCGA-PRAD_PatCount", "TCGA-SKCM_PatCount",
                    "mean_N/T", "median_N/T", "std_N/T"]
    
    # exact values are stored in output file for each chromosome
    coverage_values = coverage_data.drop(columns=meta_columns)
    coverage_values.to_csv(directory+chromosome+"_ReadDepth.tsv", sep="\t")

    # metadata is stored in list to merge later to big file
    coverage_meta = coverage_data.reindex(columns=meta_columns)
    coverage_meta_list.append(coverage_meta)

print("Writing stats to output ...")
coverage_meta = pd.concat(coverage_meta_list)
coverage_meta.to_csv(directory+"ReadDepthMeta.tsv", sep="\t")
        
print("Done.")
        

