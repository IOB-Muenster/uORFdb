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
# uORF-related variant analysis, uORF effect prediction
#
# ==========================================================================================

# ---------------------------------------------------------------------------------
# PACKAGES
# ---------------------------------------------------------------------------------
import pandas as pd
from natsort import natsorted
from pyfaidx import Fasta
from add_uORFdata import *
from analyze_uORFvar import *
from pandarallel import pandarallel
import argparse
import time
import numpy as np
from tqdm import tqdm
tqdm.pandas()

# ---------------------------------------------------------------------------------
# DEFINE INPUT / OUTPUT
# ---------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("-o","--output", help="output directory")
parser.add_argument("-i","--input", help="directory merged variants (produced by 03_AlleleFreqMapping.py)")
parser.add_argument("--uorf", help="tsv with uORF data")
parser.add_argument("-f", "--fasta", help="genomic fasta file")
parser.add_argument("--reverse", action="store_true", help="reverse chromosome order for testing (smallest first for less runtime)")
parser.add_argument("--chrom", help="selected chromosome to run", default="NONE")
args = parser.parse_args()

output_dir = args.output
input_dir = args.input
uORF_data_file = args.uorf
genomic_fasta = args.fasta
selected_chr = args.chrom
reverse = args.reverse

# ---------------------------------------------------------------------------------
# INFO
# ---------------------------------------------------------------------------------
print("+-----------------------------------+")
print("|     uORF VARIANT ANALYSIS 1.0     |")
print("+-----------------------------------+\n")
print("input directory: "+input_dir)
print("output directory: "+output_dir)
print("uORF data: "+uORF_data_file)
print("Reference genome: "+genomic_fasta)
pandarallel.initialize(progress_bar=True, nb_workers=8)
print("")

# ---------------------------------------------------------------------------------
# Reformat input
# ---------------------------------------------------------------------------------
print("Loading uORF data ...")
uORF_data = pd.read_csv(uORF_data_file, delimiter = "\t")
uORF_data.fillna("NA", inplace=True)

print("Loading sequence data ...")
genome_sequences = Fasta(genomic_fasta)

chrom_alias = {"1" : "NC_000001.11","2" : "NC_000002.12","3" : "NC_000003.12","4" : "NC_000004.12","5" : "NC_000005.10","6" : "NC_000006.12",
               "7" : "NC_000007.14","8" : "NC_000008.11","9" : "NC_000009.12","10" : "NC_000010.11","11" : "NC_000011.10","12" : "NC_000012.12",
               "13" : "NC_000013.11","14" : "NC_000014.9","15" : "NC_000015.10","16" : "NC_000016.10","17" : "NC_000017.11","18" : "NC_000018.10",
               "19" : "NC_000019.10","20" : "NC_000020.11","21" : "NC_000021.9","22" : "NC_000022.11","23" : "NC_000023.11","24" : "NC_000024.10"}

if selected_chr != "NONE":
    chromosomes = []
    for chrom in selected_chr.split(","):
        chromosomes.append(chrom_alias[chrom])
else:
    chromosomes = natsorted(uORF_data["chrom_id"].unique(), reverse=reverse)

print("Selected chromosomes: ")
print(chromosomes)
print("")

# ---------------------------------------------------------------------------------
# ITERATE OVER CHROMOSOMES
# ---------------------------------------------------------------------------------
for chromosome in chromosomes:

    print("")
    print(chromosome)
    start_time = time.time()

    print("    Loading variant data ...")
    uORF_dict = uORF_data[uORF_data["chrom_id"] == chromosome].set_index("uORF_uniqueID").to_dict("index")
    variants = input_dir+"/"+chromosome+"_VariantsWithAF.tsv"
    variants = pd.read_csv(variants, delimiter = "\t")
    
    # -----------------------------------------------------------------------------
    # ADD uORF IDs
    # -----------------------------------------------------------------------------
    print("    Adding uORF IDs ...")
    chrom_uORF_data = uORF_data[uORF_data["chrom_id"] == chromosome]
    variants = variants.parallel_apply(add_uORF_IDs, uORF_data=chrom_uORF_data, axis = 1)
    
    # -----------------------------------------------------------------------------
    # SELECTION: REMOVE VARIANTS THAT ARE PARTLY IN INTRON
    # -----------------------------------------------------------------------------
    intron_variants = variants[(variants["variant_location"] == "partly_not_in_uORF")]
    print("    Selection: "+str(len(intron_variants))+" variants were removed because they were not completely in uORFexon.")
    
    # -----------------------------------------------------------------------------
    # EXPAND ROWS OF SELECTED VARIANTS
    # -----------------------------------------------------------------------------
    variants = variants[(variants["variant_location"] != "partly_not_in_uORF")]
    
    print("    Expanding rows for each uORF-variant combination ...")
    variants = variants.explode("uORF_uniqueID").reset_index(drop=True)
    
    # -----------------------------------------------------------------------------
    # ADD uORF DATA AND ANALYSE uORF VARIANTS
    # -----------------------------------------------------------------------------
    print("    Adding uORF data ...")
    variants = variants.progress_apply(add_uORF_data, uORF_data=uORF_dict, genome_sequences=genome_sequences, axis = 1)

    print("    Analyze uORF substitutions ...")
    substitutions = variants[(variants.mutation_type == "point_mutation") | (variants.mutation_type == "multiple_bp_substitution")]
    substitutions = substitutions.progress_apply(analyze_substitutions, axis = 1)
    
    print("    Adding uORF indel data ...")
    indels = variants[(variants.mutation_type != "point_mutation") & (variants.mutation_type != "multiple_bp_substitution")]
    indels = indels.progress_apply(analyze_indels, axis = 1)
    
    sorted_columns = ["CHROM", "chrom", "HGVS", "POS", "REF", "ALT", "codon_variant", "dbSNP_ID_POS", "ClinVar_ID_POS", "dbSNP_ID_EXACT", "ClinVar_ID_EXACT",
                      "ExAC_fromVCF","GnomAD_fromVCF","TOPMED_fromVCF",
                      "ExAC.1_ref/alt/total","GnomAD.3_ref/alt/total", "TOPMED.2_ref/alt/total", "TOPMED.3_ref/alt/total",
                      "N_GT", "T_GT", "N_REF_AlleleDepth", "T_REF_AlleleDepth",
                      "N_ALT_AlleleDepth", "T_ALT_AlleleDepth", "N_ReadDepth", "T_ReadDepth", "N_SomaticMutationFrequency", "T_SomaticMutationFrequency",
                      "Proportion_T_SMF/N_SMF","CUSTOM_FILTER",
                      "BAM_ReadDepth_N", "BAM_ReadDepth_T",
                      
                      "ReadDepth_mean_N/T", "ReadDepth_median_N/T", "ReadDepth_std_N/T",
                      "VCF_PatCount",
                      "TCGA_ALT_AlleleCount","TCGA-BRCA_ALT_AlleleCount", "TCGA-COAD_ALT_AlleleCount", "TCGA-LAML_ALT_AlleleCount",
                      "TCGA-LUAD_ALT_AlleleCount", "TCGA-PRAD_ALT_AlleleCount", "TCGA-SKCM_ALT_AlleleCount",
                      "TCGA_AlleleCount","TCGA-BRCA_AlleleCount", "TCGA-COAD_AlleleCount", "TCGA-LAML_AlleleCount",
                      "TCGA-LUAD_AlleleCount", "TCGA-PRAD_AlleleCount", "TCGA-SKCM_AlleleCount",
                      "TCGA_AF","TCGA-BRCA_AF", "TCGA-COAD_AF", "TCGA-LAML_AF",
                      "TCGA-LUAD_AF", "TCGA-PRAD_AF", "TCGA-SKCM_AF",
                      
                      "cancer_type", "patient", "T_sample_type", "N_UUID", "T_UUID",
                      "strand", "gene_symbol",
                      "id", "tv_count_for_gene", "tv_count_with_same_uORF", "TV_with_uORF/TV_total", "GRCh38_REF",
                      "comparison_GRCh38/REF/ALT", "mutation_type", "cds_info", "uORFdb_ID", "variant-uorf-id", "variant_location_mRNA", "uORF_start_effect",
                      "uORF_stop_effect", "uORF_kozak_effect", "uORF_seq_effect", "REF_uORF", "ALT_uORF", "REF_start_codon", "ALT_start_codon", "REF_stop_codon",
                      "ALT_stop_codon", "REF_kozak_context", "ALT_kozak_context", "REF_kozak_strength", "ALT_kozak_strength", "REF_stop_located_in","ALT_stop_located_in",
                      "REF_uORF_length", "ALT_uORF_length", "REF_CDS_distance", "ALT_CDS_distance",
                      "DEBUG_add_uORF_data","DEBUG_analyze_substitutions","#fasta_sequences",
                      "CDS_in_other_transcripts", "REF_uORF_type", "uORF_type_effect", "gender"] #added gender on 17.08.22
    
    substitutions = substitutions.reindex(columns=sorted_columns)
    indels = indels.reindex(columns=sorted_columns)

    print("    Writing debug data to output ...")
    debug_data = substitutions[(substitutions.DEBUG_add_uORF_data != "") | (substitutions.DEBUG_analyze_substitutions != "")]
    debug_data.to_csv(output_dir+"/Intermediates/"+chromosome+"_DebugData.tsv", sep = "\t", index = False)
    
    print("    Writing discarded variants to output ...")
    discarded_variants = pd.concat([indels,intron_variants])
    indels.to_csv(output_dir+"/Intermediates/"+chromosome+"_DiscardedVariants.tsv", sep = "\t", index = False)
    
    print("    Writing substitutions to output ...")
    substitutions.drop(columns=["DEBUG_add_uORF_data","DEBUG_analyze_substitutions","#fasta_sequences"]).to_csv(output_dir+"/"+chromosome+"_AnalyzedSubstitutions.tsv", sep = "\t", index = False)

    print("\n"+str(round(((time.time() - start_time)/ 60),2))+" min")

print("Done")
        