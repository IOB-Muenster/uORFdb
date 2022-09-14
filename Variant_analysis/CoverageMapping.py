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
# Obtained read depths (Samtools) is mapped to variants according to POS and CHROM
#
# ==========================================================================================

# ---------------------------------------------------------------------------------
# PACKAGES
# ---------------------------------------------------------------------------------
import pandas as pd
from natsort import natsorted
import time
import argparse
from tqdm import tqdm
tqdm.pandas()
from pandarallel import pandarallel
pandarallel.initialize(progress_bar=True)

# ---------------------------------------------------------------------------------
# DEFINE INPUT / OUTPUT
# ---------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("-d","--dir", help="working directory")
parser.add_argument("-i","--cov_dir", help="directory with preprocessed coverage data")
parser.add_argument("-v","--variants", help="files with variants")
args = parser.parse_args()

cov_dir = args.cov_dir
directory = args.dir
variants = args.variants

# ---------------------------------------------------------------------------------
# DEFINE INPUT / OUTPUT
# ---------------------------------------------------------------------------------
variants = pd.read_csv(variants, delimiter = "\t")

intermediate_file = directory+"/VariantsWithCoverage_intermediate.tsv"
outfile = directory+"/VariantsWithCoverage.tsv"

coverage_stats = cov_dir+"/ReadDepthMeta.tsv"
coverage_stats = pd.read_csv(coverage_stats, delimiter = "\t", index_col=0)
coverage_stats = coverage_stats.to_dict("index")

# ---------------------------------------------------------------------------------
# MAP PATIENT COUNTS AND STATS
# ---------------------------------------------------------------------------------
def map_coverage_meta(var:pd.Series, coverage_stats):
    
    count_columns = ["TCGA-BRCA_PatCount", "TCGA-COAD_PatCount", "TCGA-LAML_PatCount",
                      "TCGA-LUAD_PatCount", "TCGA-PRAD_PatCount", "TCGA-SKCM_PatCount"]
    stat_columns = ["mean_N/T", "median_N/T", "std_N/T"]

    # select best dataset (the one with the highest patient count, if same take first position)
    
    selected_TCGA_PatCount = 0

    for POS in range(var["POS"],var["TO_0BHO"]+1):
        
        CHR_POS = var["CHROM"]+"-"+str(POS)

        TCGA_PatCount = 0
        for cancer_type in count_columns:
            TCGA_PatCount += int(coverage_stats[CHR_POS][cancer_type])
            
        if TCGA_PatCount > selected_TCGA_PatCount:
            selected_TCGA_PatCount = TCGA_PatCount
            selected_CHR_POS = CHR_POS
     
    # calculate to allele counts       
    for column in stat_columns:
        var["ReadDepth_"+column] = coverage_stats[selected_CHR_POS][column]
    for column in count_columns:
        var[column[:9]+"_AlleleCount"] = int(coverage_stats[selected_CHR_POS][column])*2
    
    var["TCGA_AlleleCount"] = selected_TCGA_PatCount*2
            
    return var

variants_with_patcount = variants.parallel_apply(map_coverage_meta, coverage_stats=coverage_stats, axis=1)

variants_with_patcount.to_csv(intermediate_file, sep="\t", index=False)


# ---------------------------------------------------------------------------------
# MAP EXACT READ DEPTH TO RESPECTIVE VARIANTS
# ---------------------------------------------------------------------------------
def map_ReadDepth(var:pd.Series, coverage_data):
    
    cancer_types = ["TCGA-BRCA", "TCGA-COAD", "TCGA-LAML",
                    "TCGA-LUAD", "TCGA-PRAD", "TCGA-SKCM"]
    
    var["BAM_ReadDepth_N"] = []
    var["BAM_ReadDepth_T"] = []
    
    ID = var["patient"]+"_"+var["N_UUID"]+"_"+var["T_UUID"]

    for POS in range(var["POS"],var["TO_0BHO"]+1):

        CHR_POS = var["CHROM"]+"-"+str(POS)
        
        if CHR_POS in coverage_data and ID in coverage_data[CHR_POS]:

            read_depth = coverage_data[CHR_POS][ID].split(",")

            var["BAM_ReadDepth_N"].append(int(read_depth[0]))
            var["BAM_ReadDepth_T"].append(int(read_depth[1]))
    
    if var["BAM_ReadDepth_N"] != [] and var["BAM_ReadDepth_T"] != []:
        if max(var["BAM_ReadDepth_N"]) >= 10 and max(var["BAM_ReadDepth_T"]) >= 10:
            var["TCGA_ALT_AlleleCount"] = var["T_GT"].count("1")
        else:
            var["TCGA_ALT_AlleleCount"] = 0
    else:
        var["TCGA_ALT_AlleleCount"] = "."
        print("no RD found")

    if var["BAM_ReadDepth_N"] == []:
        var["BAM_ReadDepth_N"] = "."  
    else:
        var["BAM_ReadDepth_N"] = max(var["BAM_ReadDepth_N"])
        
    if var["BAM_ReadDepth_T"] == []:
        var["BAM_ReadDepth_T"] = "."
    else:
        var["BAM_ReadDepth_T"] = max(var["BAM_ReadDepth_T"])
        
    for cancer_type in cancer_types:
        if cancer_type == var["cancer_type"]:
            var[cancer_type+"_ALT_AlleleCount"] = var["TCGA_ALT_AlleleCount"]
        else:
            var[cancer_type+"_ALT_AlleleCount"] = 0
       
    return var

chromosomes = natsorted(variants["CHROM"].unique())
all_variants_with_coverage = []

for chromosome in chromosomes:
    
    print("")
    print(chromosome)
    
    coverage_data = cov_dir+"/"+chromosome+"_ReadDepth.tsv"
    coverage_data = pd.read_csv(coverage_data, delimiter = "\t", index_col=0)
    coverage_data = coverage_data.to_dict("index")
    
    variants_with_coverage = variants_with_patcount[(variants_with_patcount["CHROM"] == chromosome)].progress_apply(map_ReadDepth, coverage_data=coverage_data, axis=1)
    
    all_variants_with_coverage.append(variants_with_coverage)
    
all_variants_with_coverage = pd.concat(all_variants_with_coverage)
 
all_variants_with_coverage = all_variants_with_coverage.reindex(
    columns=["CHROM", "POS", "FROM_0BHO", "TO_0BHO", "REF", "ALT", "dbSNP_ID_POS", "ClinVar_ID_POS", "dbSNP_ID_EXACT", "ClinVar_ID_EXACT",
    "ExAC_fromVCF", "GnomAD_fromVCF", "TOPMED_fromVCF", "N_GT", "T_GT",
    "N_REF_AlleleDepth", "T_REF_AlleleDepth", "N_ALT_AlleleDepth", "T_ALT_AlleleDepth", "N_ReadDepth", "T_ReadDepth",
    "N_SomaticMutationFrequency", "T_SomaticMutationFrequency", "Proportion_T_SMF/N_SMF", "CUSTOM_FILTER",
    "BAM_ReadDepth_N", "BAM_ReadDepth_T",
    "TCGA_ALT_AlleleCount",
    "TCGA-BRCA_ALT_AlleleCount", "TCGA-COAD_ALT_AlleleCount", "TCGA-LAML_ALT_AlleleCount", "TCGA-LUAD_ALT_AlleleCount", "TCGA-PRAD_ALT_AlleleCount", "TCGA-SKCM_ALT_AlleleCount",
    "TCGA_AlleleCount","TCGA-BRCA_AlleleCount", "TCGA-COAD_AlleleCount", "TCGA-LAML_AlleleCount", "TCGA-LUAD_AlleleCount", "TCGA-PRAD_AlleleCount", "TCGA-SKCM_AlleleCount",
    "ReadDepth_mean_N/T", "ReadDepth_median_N/T", "ReadDepth_std_N/T",
    "cancer_type", "patient", "N_sample_type", "T_sample_type", "N_UUID", "T_UUID", "tumor_sample", "normal_sample", "file", "gender"])

all_variants_with_coverage.to_csv(outfile, sep="\t", index=False)
    