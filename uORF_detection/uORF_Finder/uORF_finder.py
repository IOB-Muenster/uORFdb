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
# uORF finder tool, BASIC_FUNCTIONS and PANDA_FUNCTIONS are necessary
#
# ==========================================================================================

# ---------------------------------------------------------------------------------
# PACKAGES AND FUNCTIONS
# ---------------------------------------------------------------------------------
import argparse
import pandas as pd
import re
import time
from PANDA_FUNCTIONS import *
from BASIC_FUNCTIONS import *
import time
from pandarallel import pandarallel
from natsort import natsorted

# ---------------------------------------------------------------------------------
# DEFINE INPUT / OUTPUT
# ---------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--refseq",help="path to refseq file, required")
parser.add_argument("--output", help="path to output directory, required")
parser.add_argument("--fasta", help="path to fasta file that contains seq of all chromosomes, required")

parser.add_argument("-v", "--verbose", action="store_true", help="verbose output, default = False")
parser.add_argument("-p", "--preprocessing", action="store_true", help="only run preprocessing, default = False")
parser.add_argument("-d", "--uorfdetection", action="store_true", help="only run uORF detection, default = False")
parser.add_argument("-c", "--chrom_output", action="store_true", help="enable one output file per chromosome, recommended!, default = False")
parser.add_argument("--complement", action="store_true", help="only run chromosomes which were not processed yet (--chrom_output is set to True as well), default = False")
parser.add_argument("--duplicate", action="store_true", help="includes duplicate transcript ids on different chromosomes, default = False")
parser.add_argument("--test21", action="store_true", help="test run chr21 (does not work for all species, recommended for Homo sapiens)")

args = parser.parse_args()

RefSeq_file = args.refseq
genomic_fasta = args.fasta
output_dir = args.output
duplicate = args.duplicate
verbose = args.verbose
chrom_output = args.chrom_output
preprocessing = args.preprocessing
uorfdetection = args.uorfdetection
complement = args.complement
test21 = args.test21

big_output = output_dir+"/uORF_data_v2.tsv"
transcripts_tsv = output_dir+"/transcripts.tsv"

# ---------------------------------------------------------------------------------
# SETTINGS
# ---------------------------------------------------------------------------------
print("")
print("uORF FINDER version 2.2\n")
print("SETTINGS")

if duplicate == False:
    print("duplicate regions: excluded")
elif duplicate == True:
    print("duplicate regions: included")
    
if verbose == True:
    print("output: verbose")
elif verbose == False:
    print("ouput: default")
    
if preprocessing == False and uorfdetection == False:
    preprocessing = True
    uorfdetection = True

if preprocessing == True:
    print("preprocessing: enabled")
elif preprocessing == False:
    print("preprocessing: disabled")
    
if uorfdetection == True:
    print("uorfdetection: enabled")
elif uorfdetection == False:
    print("uorfdetection: disabled")
    
if complement == True:
    print("complement: enabled")
    chrom_output = True
elif complement == False:
    print("complement: disabled")
    
if chrom_output == True:
    print("output file: additional tsv for each chromosome")
elif chrom_output == False:
    print("ouput file: one tsv for all chromosomes")
    
    
pandarallel.initialize(progress_bar=True)
print("")

# ---------------------------------------------------------------------------------
# PREPROCESSING
# ---------------------------------------------------------------------------------
if preprocessing == True:
    print("PREPROCESSING")
    # ---------------------------------------------------------------------------------
    # CONVERT INPUT
    # ---------------------------------------------------------------------------------
    column_names = ["bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames"]
    transcripts = pd.read_csv(RefSeq_file, delimiter = "\t", names = column_names, header = None, keep_default_na=False)
    transcripts.fillna("nan", inplace=True)
    transcripts["index"] = transcripts.index
    
    chromosomes = headers_to_IDs(genomic_fasta)

    print("")
    print("Extracted chromosomes:")
    for chromosome in chromosomes:
        print(chromosome+" -> "+chromosomes[chromosome])
    print("")
    
    # ---------------------------------------------------------------------------------
    # PREPROCESSING
    # ---------------------------------------------------------------------------------
    print("Filtering transcripts ... ")
    start_time = time.time()
    if test21 == True:
        transcripts = transcripts[(transcripts.chrom == "chr21")]
    transcripts = transcripts.parallel_apply(filter_RefSeq, chromosomes = chromosomes, axis = 1)
    print("\n"+str(round(((time.time() - start_time)/ 60),2))+" min")

    # mark duplicates
    transcripts["duplicate"] = "."
    duplications = transcripts.duplicated("duplicate_candidate", keep = False)
    duplicated_transcripts = transcripts.copy()[(duplications == True) & (transcripts.duplicate_candidate != ".")]
    duplicated_transcripts["duplicate"] = "duplicate"
    not_duplicated_transcripts = transcripts.loc[(duplications != True) | (transcripts.duplicate_candidate == ".")]
    transcripts = pd.concat([duplicated_transcripts,not_duplicated_transcripts])
    transcripts = transcripts.sort_values(by=["index"])
    transcripts = transcripts.drop(columns=["duplicate_candidate", "index"])
    
    # transcripts with missing 5'UTR are counted for transcript count for each gene too!
    print("Adding transcript counts ... ")
    start_time = time.time()
    transcripts = transcripts.parallel_apply(add_tv_counts, dataframe=transcripts, axis = 1)
    print("\n"+str(round(((time.time() - start_time)/ 60),2))+" min")

    print("Extracting sequences ... ")
    start_time = time.time()
    transcripts = transcripts.apply(extract_sequences, genomic_fasta=genomic_fasta, axis = 1)
    print("\n"+str(round(((time.time() - start_time)/ 60),2))+" min")
    
    sorted_columns = ["bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount",
                      "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames",
                      "fasta_header", "RtxStart", "RtxEnd", "RcdsStart", "RcdsEnd",
                      "RexonStarts", "RexonEnds", "R5utrStart", "R5utrEnd", "R5utrLength",
                      "CDS_kozak_context", "CDS_kozak_strength", "tv_count_for_gene", "mRNA_seq", "selection_status","duplicate"]
    transcripts = transcripts.reindex(columns=sorted_columns)

    # ---------------------------------------------------------------------------------
    # WRITE TO OUTPUT
    # ---------------------------------------------------------------------------------
    print("Writing to output ... ", end = '')
    start_time = time.time()
    transcripts.to_csv(transcripts_tsv, sep="\t", index = False)
    print(str(round(((time.time() - start_time)/ 60),2))+" min")
    print("\n")

# ---------------------------------------------------------------------------------
# UORF DETECTION
# ---------------------------------------------------------------------------------
if uorfdetection == True:
    print("UORF DETECTION")
    
    # ---------------------------------------------------------------------------------
    # INPUT SECTION
    # ---------------------------------------------------------------------------------
    transcripts = pd.read_csv(transcripts_tsv, delimiter = "\t")
    
    selected_transcripts = transcripts[(transcripts.selection_status == "SELECTED")]
    selected_transcripts = selected_transcripts.drop(columns=["bin", "score", "cdsStartStat", "cdsEndStat", "exonFrames"])

    chromosomes = natsorted(selected_transcripts["chrom"].unique())
    
    start_triplets = ["ATG", "AAG", "ATT", "TTG", "AGG", "ATC", "GTG", "ACG", "CTG", "ATA"]
    stop_triplets = ["TAG", "TGA", "TAA"]
    
    # ---------------------------------------------------------------------------------
    # MAKE DICT WITH CDSSTARTS
    # ---------------------------------------------------------------------------------
    plus_transcripts_dict = {}
    minus_transcripts_dict = {}

    transcripts = transcripts[(transcripts.selection_status == "SELECTED") | (transcripts.selection_status == "no_5'UTR,")]
    
    for index, row in transcripts.iterrows():
        if row["strand"] == "+":
            if row["cdsStart"] not in plus_transcripts_dict:
                plus_transcripts_dict[row["cdsStart"]] = [row["name"]+":"+row["name2"]]
            else:
                plus_transcripts_dict[row["cdsStart"]].append(row["name"]+":"+row["name2"])
            
        elif row["strand"] == "-":
            if row["cdsEnd"] not in minus_transcripts_dict:
                minus_transcripts_dict[row["cdsEnd"]] = [row["name"]+":"+row["name2"]]
            else:
                minus_transcripts_dict[row["cdsEnd"]].append(row["name"]+":"+row["name2"])

    # ---------------------------------------------------------------------------------
    # uORF DETECTION
    # ---------------------------------------------------------------------------------
    uORF_data_all_chromosomes = []
    
    whole_start_time = time.time()
    
    for chrom in chromosomes:
        
        output_file = output_dir+"/"+chrom+"_uORF_data.tsv"
        
        if complement == True:
            try:
                with open(output_file, "r") as file:
                    print(chrom+" "+output_file+" already exists. Chromosome is skipped.")
                    print(chrom+" Appending data ... ")
                    uORF_data = pd.read_csv(output_file, delimiter = "\t")
                    uORF_data_all_chromosomes.append(uORF_data)
                    continue
            except:
                print(chrom+" Starting ....")
        
        print(chrom+" Detecting uORFs ... ")
        start_time = time.time()
        uORF_data = selected_transcripts[(selected_transcripts.chrom == chrom)].apply(detect_start_pos, axis = 1)

        # DEBUG
        if len(selected_transcripts[(selected_transcripts.chrom == chrom)]) == 0:
            print("No transcripts selected. Proceeding with next chromosome.")
            break

        uORF_data = uORF_data[(uORF_data.start_pos_mRNA != ".")].explode("start_pos_mRNA").reset_index(drop=True)
        uORF_data = uORF_data.parallel_apply(detect_uORFs, minus_transcripts_dict=minus_transcripts_dict, plus_transcripts_dict=plus_transcripts_dict, axis = 1)
        print("\n"+str(round(((time.time() - start_time)/ 60),2))+" min")
        
        uORF_data = uORF_data.drop(columns=["fasta_header","name2","name","selection_status"])
        
        sorted_columns = [ "uORF_uniqueID", "uORFdb_ID", "uORF_ID",
            "chrom_id","chrom", "id", "strand", "gene_symbol", "start", "stop", "start_codon", "stop_codon", "start_pos_mRNA", "end_pos_mRNA",
            "uorf_length_incl_introns", "uorf_length_excl_introns", "cds_start_dist", "cds_start_dist_no_introns", "stop_located_in", "frame",
            "uORF_kozak_context", "uORF_kozak_fragment", "uORF_kozak_strength", "CDS_kozak_context", "CDS_kozak_strength", "seq_no_introns", "translated_uORF",
            "tv_count_for_gene", "tv_count_with_same_uORF", "TV_with_uORF/TV_total", "uORF_type", "CDS_in_other_transcripts", "CDS_stop_info", "duplicate",
            "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "uORF_exonStarts", "uORF_exonEnds", 
            "RtxStart", "RtxEnd", "RcdsStart", "RcdsEnd", "RexonStarts", "RexonEnds", "R5utrStart","R5utrEnd", "R5utrLength", "mRNA_seq",
            "start_incl_kozak","stop_incl_kozak","kozak_status"]
        uORF_data = uORF_data.reindex(columns=sorted_columns)
        
        print(chrom+" Adding transcript counts ... ")
        start_time = time.time()
        uORF_data = uORF_data.parallel_apply(add_transcript_and_gene_data, dataframe = uORF_data, axis = 1)
        print("\n"+str(round(((time.time() - start_time)/ 60),2))+" min")

        if verbose == False:
            uORF_data = uORF_data.drop(columns=["uORF_ID","txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
                "RtxStart", "RtxEnd", "RcdsStart", "RcdsEnd", "RexonStarts", "RexonEnds", "R5utrStart","R5utrEnd", "R5utrLength", "mRNA_seq",
                "start_incl_kozak","stop_incl_kozak","kozak_status", "uORF_kozak_fragment", "uORF_exonStarts", "uORF_exonEnds"])
            
        if duplicate == False:
            uORF_data = uORF_data[uORF_data["duplicate"] != "duplicate"]
            uORF_data = uORF_data.drop(columns=["duplicate"])
        
        print(chrom+" Appending data ... ")
        uORF_data_all_chromosomes.append(uORF_data)
        
        if chrom_output == True:
            print(chrom+" Writing to output ... ")
            uORF_data.to_csv(output_file, sep="\t", index = False)

        print("\n")
        
    uORF_data_all_chromosomes_df = pd.concat(uORF_data_all_chromosomes)
    uORF_data_all_chromosomes_df.to_csv(big_output, sep="\t", index = False)
        
    print("whole time: "+str(round(((time.time() - whole_start_time)/ 60),2))+" min")

