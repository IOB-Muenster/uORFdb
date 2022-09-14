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
# Functions used for uORF finder
#
# ==========================================================================================

# ---------------------------------------------------------------------------------
# FUNCTIONS
# ---------------------------------------------------------------------------------
import pandas as pd
from Bio.Seq import Seq
import hashlib
from Bio import SeqIO

def extract_chrID(fasta_header):
    import re
    try:
        chrID = ((re.search('(?:chromosome )([0-9 a-z A-Z _]*)', fasta_header)).group(0)).replace("omosome ","")
    except:
        chrID = "chr_ID_not_found"

    print(fasta_header+" -> "+chrID)
    return(chrID)

def headers_to_IDs(chrFile):

    chromosomes = {}

    for record in SeqIO.parse(chrFile,"fasta"):
        header = record.description
    
        if "scaffold" not in header and "chromosome" in header and "NW" not in header:

            chromosomes[extract_chrID(header)] = record.id
        elif "scaffold" not in header and "chr" in header and "NW" not in header and "M" not in header and "random" not in header and "Un" not in header and "alt" not in header:

            chromosomes[record.id] = record.id

    return chromosomes

def string_to_list(string):
    my_list = list(string.split(","))
    my_list.pop()
    my_list = list(map(int, my_list))
    return(my_list)

def list_to_string(my_list):
    string = ""
    for i in my_list:
        string += (str(i)+",")
    return(string)

def calculate_relative_coordinates(exonStarts,exonEnds,strand,coordinate):
    
    #input: exonStart and exonEnd coordinates + one coordinate that should be calculated as relative coordinate in this mRNA, respect of minus strand
    exons_sum = 0
    
    for i in range(0, len(exonStarts)):

        if coordinate in range(exonStarts[i],exonEnds[i]+1):#coordinate is in exon
            relative_coordinate = (exons_sum + (coordinate - exonStarts[i]))
            exons_sum += (exonEnds[i] - exonStarts[i])
        
        elif i != (len(exonStarts)-1) and coordinate in range(exonEnds[i],exonStarts[i+1]+1):
            exons_sum += (exonEnds[i] - exonStarts[i])
            relative_coordinate = exons_sum
            
        else:
            exons_sum += (exonEnds[i] - exonStarts[i])
    
    #reverse for minus strand
    if strand == "-":
        relative_coordinate = exons_sum - relative_coordinate
    
    return(relative_coordinate)


def frequencyCount(string, substr): #only needed for function "get_codon_positions"
    count = 0
    pos = 0
    while(True):
        pos = string.find(substr, pos)
        if pos > -1:
            count = count + 1
            pos += 1
        else:
            break
    return count


def get_codon_positions(codons,seq,start=0,end=None):
    
    if end == None:#if no argument for end, set to length of the entered sequence
        end = len(seq)
    
    codon_positions = []
    
    for codon in codons:
        
        codon_counts = frequencyCount(seq,codon) #counts how often specific codon-triplet is in sequence, if nothing is found, it does not run for loop
        
        #start is usually 0
        pos = start - 1

        for i in range(0, codon_counts):
            pos = seq.find(codon,pos+1,end)
            codon_positions.append(pos)
    
    codon_positions.sort()
            
    return codon_positions

def calculate_absolute_coordinates(exonStarts,exonEnds,RexonStarts,RexonEnds,strand,relative_coordinate):
    
    for i in range(0,len(RexonStarts)):
        
        if relative_coordinate in range(RexonStarts[i],RexonEnds[i]+1):
            
            if strand == "-":
                absolute_coordinate = exonEnds[-(i+1)] - (relative_coordinate - RexonStarts[i])
                
            elif strand == "+":
                absolute_coordinate = exonStarts[i] + (relative_coordinate - RexonStarts[i])
        
    return(absolute_coordinate)

def stop_location(stop_pos,R5utrStart,R5utrEnd,RcdsStart,RcdsEnd,RtxEnd):
    
    if (stop_pos + 2) in range(R5utrStart,R5utrEnd):
        stop_located_in = "5'UTR"
        
    elif (stop_pos + 2) in range(RcdsStart,RcdsEnd):
        
        if stop_pos in range(R5utrStart,R5utrEnd):
            stop_located_in = "5'UTR_CDS_overlap"
        
        else:
            stop_located_in = "CDS"
        
    elif (stop_pos + 2) in range(RcdsEnd,RtxEnd):
        
        if stop_pos in range(RcdsStart,RcdsEnd):
            stop_located_in = "CDS_3'UTR_overlap"
            
        elif (stop_pos + 2) == (RtxEnd - 1):
            stop_located_in = "transcript_end"
            
        else:
            stop_located_in = "3'UTR"
    
    return(stop_located_in)


def determine_kozak_strength(kozak_context):
    if (kozak_context[3] == "G" or kozak_context[3] == "A") and (kozak_context[9] == "G"):
        kozak_strength = "strong"
    elif (kozak_context[3] == "G" or kozak_context[3] == "A") and (kozak_context[9] != "G"): 
        kozak_strength = "adequate"
    elif kozak_context[3] != "G" and kozak_context[3] != "A" and (kozak_context[9] == "G"): 
        kozak_strength = "adequate"
    elif kozak_context[3] != "G" and kozak_context[3] != "A" and kozak_context[9] != "G":
        kozak_strength = "weak"

    return kozak_strength
