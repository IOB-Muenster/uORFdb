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
# get relevant IDs from dbSNP and ClinVar VCFs for selected uORF exon regions
#
# ==========================================================================================

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input",help="input file")
parser.add_argument("--output",help="output file")
args = parser.parse_args()
vcf = args.input
regions_file = args.output

text = "#CHROM\tFROM\tTO\tID\n"

with open (vcf, "r") as file:
    for line in file:
        if "##" not in line and "#CHROM" not in line:
            line = line.split("\t")
            
            CHROM = line[0]
            POS = int(line[1])
            ID = "POS-"+line[2]
            REF = line[3]
            
            FROM = str(POS)
            TO = str(POS+(len(REF)-1))
            
            text += CHROM+"\t"+FROM+"\t"+TO+"\t"+ID+"\n"
            
with open (regions_file, "w") as file:
    file.write(text)
    
    
    