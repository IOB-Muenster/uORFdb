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
# get uORF exon regions from uORF finder output
#
# convert coordinates from half-open, 0-based to fully-closed, 1-based
# https://genome-blog.soe.ucsc.edu/blog/2016/12/12/the-ucsc-genome-browser-coordinate-counting-systems/
#
# ==========================================================================================

import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--uORF_data",help="input file (uORF finder output table)")
parser.add_argument("--uORF_regions_file",help="output file (table with relevant uORF coordinates)")
args = parser.parse_args()
uORF_data = args.uORF_data
uORF_regions_file = args.uORF_regions_file

uORF_dataframe = pd.read_csv(uORF_data, delimiter = "\t")

def string_to_list(string):
    my_list = list(string.split(","))
    my_list.pop()
    my_list = list(map(int, my_list))
    return(my_list)

file_text = ""
chrom = ""

for index, row in uORF_dataframe[["chrom_id","uORF_exonStarts","uORF_exonEnds"]].iterrows():
    
    if chrom != row["chrom_id"]:
        print(row["chrom_id"])
    
    chrom = row["chrom_id"]
    
    #uORF_exon coordinates contain kozak context!
    uORF_exonStarts = string_to_list(row["uORF_exonStarts"])
    uORF_exonEnds = string_to_list(row["uORF_exonEnds"])
    
    for i in range(len(uORF_exonStarts)):
        start = uORF_exonStarts[i]+1
        stop = uORF_exonEnds[i]

        file_text += chrom+"\t"+str(start)+"\t"+str(stop)+"\n"
        
with open (uORF_regions_file, "w") as file:
    file.write(file_text)
    
print("Done.")
    