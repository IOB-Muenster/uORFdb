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
# Functions used in Variant Analysis scripts
#
# ==========================================================================================


import pandas as pd
from Bio.Seq import Seq

def remove_internal_stop(REF_triplets,ALT_triplets,seq):
    
    stop_codons = ["TAG", "TGA", "TAA"]
    
    REF_triplets_stop_removed = []
    ALT_triplets_stop_removed = []
    
    for triplet_number in range(0,len(REF_triplets)):
        
        if REF_triplets[triplet_number] not in stop_codons or triplet_number == (len(REF_triplets)-1):
            REF_triplets_stop_removed.append(REF_triplets[triplet_number])
            ALT_triplets_stop_removed.append(ALT_triplets[triplet_number])
        else:
            REF_triplets_stop_removed.append("---")
            ALT_triplets_stop_removed.append("---")
                
    if seq == "REF":
        return REF_triplets_stop_removed
    elif seq == "ALT":
        return ALT_triplets_stop_removed
    
    
def translate_to_list(SEQ,side):
    
    SEQ_triplets = []
    modulo = len(SEQ) % 3
    
    if side == "right":
        for i in range(0,len(SEQ)-modulo):
            if i % 3 == 0:
                triplet = SEQ[i:i+3]
                SEQ_triplets.append(triplet)
        
    elif side == "left":
        for i in range(modulo,len(SEQ)):
            if (i+(3-modulo)) % 3 == 0:
                triplet = SEQ[i:i+3]
                SEQ_triplets.append(triplet)
    
    return SEQ_triplets


def detect_triplet_positions(triplets,codons):
    
    start_codons = ["ATG", "AAG", "ATT", "TTG", "AGG", "ATC", "GTG", "ACG", "CTG", "ATA"]          
    stop_codons = ["TAG", "TGA", "TAA"]
    
    if codons == "start_codons":
        codon_list = start_codons
    elif codons == "stop_codons":
        codon_list = stop_codons

    codon_positions=[]
    for i in range(0,len(triplets)):
        if triplets[i] in codon_list:
            codon_positions.append(i)
            
    return codon_positions


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


def stop_location(stop_pos,R5utrStart,R5utrEnd,RcdsStart,RcdsEnd,RtxEnd):
    
    if (stop_pos + 2) in range(R5utrStart,R5utrEnd):
        stop_located_in = "5'UTR"
        
    elif (stop_pos + 2) in range(RcdsStart,RcdsEnd):
        
        if stop_pos in range(R5utrStart,R5utrEnd):
            stop_located_in = "5'UTR CDS overlap"
        
        else:
            stop_located_in = "CDS"
        
    elif (stop_pos + 2) in range(RcdsEnd,RtxEnd):
        
        if stop_pos in range(RcdsStart,RcdsEnd):
            stop_located_in = "CDS 3'UTR overlap"
            
        elif (stop_pos + 2) == (RtxEnd - 1):
            stop_located_in = "transcript end"
            
        else:
            stop_located_in = "3'UTR"
    
    return(stop_located_in)


def analyze_substitutions(var:pd.Series):

    var["DEBUG_analyze_substitutions"] = ""
    
    #-----------------------------------------------
    # BASIC DATA
    #-----------------------------------------------
    REF_uORF = var["REF_uORF"]
    ALT_uORFregion = var["ALT_uORFregion"]
    ALT_SEQregion = var["ALT_SEQregion"]

    REF_uORF_triplets = translate_to_list(REF_uORF,side="right")
    ALT_uORF_triplets = translate_to_list(ALT_uORFregion,side="right")
    
    #-----------------------------------------------
    # CODON MUTATION (Lara, 01.06.22)
    #-----------------------------------------------
    codon_variant = ""
    
    for i in range(0,len(REF_uORF_triplets)):
        
        if REF_uORF_triplets[i] != ALT_uORF_triplets[i]:
            codon_variant += REF_uORF_triplets[i]+"->"+ALT_uORF_triplets[i]+","
    
    if codon_variant == "":
        var["codon_variant"] = "."
    else:
        var["codon_variant"] = codon_variant 

    #-----------------------------------------------
    # START CODON
    #-----------------------------------------------
    var["ALT_uORF_start_pos_mRNA"] = "."
    
    start_codons = ["ATG", "AAG", "ATT", "TTG", "AGG", "ATC", "GTG", "ACG", "CTG", "ATA"]
    
    REF_start = REF_uORF_triplets[0]
    ALT_start = ALT_uORF_triplets[0]
    
    REF_start_positions = detect_triplet_positions(REF_uORF_triplets,codons="start_codons")
    ALT_start_positions = detect_triplet_positions(ALT_uORF_triplets,codons="start_codons")
    
    if REF_start == ALT_start:
        # start codon remains intact
        var["ALT_start_codon"] = "."
        var["uORF_start_effect"] = "."
        
    else:
        #check if kozak can complete the start codon
        if REF_uORF in ALT_SEQregion:
            var["DEBUG_analyze_substitutions"] += "REF_uORF is still in ALT_SEQ,"

        # determine type of start mutation
        if REF_start == "ATG":
            ref_start_type = "uAUG"
        elif REF_start != "ATG" and REF_start in start_codons:
            ref_start_type = "aTIS"
        else:
            print(REF_start)
            print(var["REF_kozak_context"])
            print(var["uORF_uniqueID"])
            
        if ALT_start == "ATG":
            alt_start_type = "uAUG"
            var["ALT_start_codon"] = ALT_start
        elif ALT_start != "ATG" and ALT_start in start_codons:
            alt_start_type = "aTIS"
            var["ALT_start_codon"] = ALT_start
        
        elif ALT_start not in start_codons:
            
            if len(ALT_start_positions) == len(REF_start_positions):
                #this means that the start codon only changed
                if ALT_start_positions[0] != REF_start_positions[0]:
                    
                    var["DEBUG_analyze_substitutions"] += "start codon shift,"
                    ALT_start = ALT_uORF_triplets[ALT_start_positions[0]]
                    if ALT_start != "ATG" and ALT_start in start_codons:
                        alt_start_type = "aTISpos"
                    elif ALT_start == "ATG":
                        alt_start_type = "uAUGpos"
                    
                    ALT_uORF = ""
                    for triplet_number in range(ALT_start_positions[0],len(ALT_uORF_triplets)):
                        ALT_uORF += ALT_uORF_triplets[triplet_number]
                    var["ALT_uORF_start_pos_mRNA"] = var["uORF_end_pos_mRNA"]-len(ALT_uORF)
                    var["ALT_start_codon"] = ALT_start
                    
            else:
                alt_start_type = "uSTARTloss"
                var["ALT_start_codon"] = "loss"
                ALT_uORF = "loss"
            
        var["uORF_start_effect"] = ref_start_type+"->"+alt_start_type
    
    #-----------------------------------------------
    # KOZAK CONTEXT
    #-----------------------------------------------
    REF_kozak_fragment = var["REF_kozak_fragment"].replace(".", "")
    ALT_kozak_fragmentregion = var["ALT_kozak_fragmentregion"].replace(".", "")

    if var["ALT_start_codon"] == "loss":
        var["ALT_kozak_context"] = "."
        var["ALT_kozak_strength"] = "."
        var["uORF_kozak_effect"] = "."
 
    else:
        # respect start codon shift
        if "pos" in var["uORF_start_effect"]:
            if var["ALT_uORF_start_pos_mRNA"] < 6:
                var["ALT_kozak_context"] = (6-var["ALT_uORF_start_pos_mRNA"])*"-" + var["VAR_mRNA_SEQ"][var["ALT_uORF_start_pos_mRNA"]:var["ALT_uORF_start_pos_mRNA"]+4]

            else:
                var["ALT_kozak_context"] = var["VAR_mRNA_SEQ"][(var["ALT_uORF_start_pos_mRNA"]-6):(var["ALT_uORF_start_pos_mRNA"]+4)]
        else:
            var["ALT_kozak_context"] = (6-len(ALT_kozak_fragmentregion))*"-"+ALT_kozak_fragmentregion+ALT_uORFregion[:4]

        if var["ALT_kozak_context"][3] != "-":
            var["ALT_kozak_strength"] = determine_kozak_strength(var["ALT_kozak_context"])
 
        else:
            var["ALT_kozak_strength"] = "."

        if (var["ALT_kozak_context"][:6]+var["ALT_kozak_context"][:10]) != (var["REF_kozak_context"][:6]+var["REF_kozak_context"][:10]) and var["ALT_kozak_strength"] == var["REF_kozak_strength"]:
            var["uORF_kozak_effect"] = "altered sequence"
         
        elif (var["ALT_kozak_context"][:6]+var["ALT_kozak_context"][:10]) != (var["REF_kozak_context"][:6]+var["REF_kozak_context"][:10]) and var["ALT_kozak_strength"] != var["REF_kozak_strength"]:
            var["uORF_kozak_effect"] =  var["REF_kozak_strength"]+"->"+var["ALT_kozak_strength"] 

        elif (var["ALT_kozak_context"][:6]+var["ALT_kozak_context"][:10]) == (var["REF_kozak_context"][:6]+var["REF_kozak_context"][:10]):
            var["uORF_kozak_effect"] = "."

    ALT_kozak_context_int = var["ALT_kozak_context"]
    if var["REF_kozak_strength"] == var["ALT_kozak_strength"]:
        var["ALT_kozak_strength"] = "."
    if var["REF_kozak_context"] == var["ALT_kozak_context"]:
        var["ALT_kozak_context"] = "."
        
    #-----------------------------------------------
    # STOP CODON
    #-----------------------------------------------
    if var["ALT_start_codon"] == "loss":
        var["uORF_stop_effect"] = "."
        var["ALT_stop_codon"] = "."
        ALT_uORF_end_pos_mRNA = "."
        
    else:
        REF_stop_positions = detect_triplet_positions(REF_uORF_triplets,codons="stop_codons")
        ALT_stop_positions = detect_triplet_positions(ALT_uORF_triplets,codons="stop_codons")
        
        if len(REF_stop_positions) > 1:
            var["DEBUG_analyze_substitutions"] += "stop codon readthrough,"

            REF_uORF_triplets_rm = remove_internal_stop(REF_triplets=REF_uORF_triplets,ALT_triplets=ALT_uORF_triplets,seq="REF")
            ALT_uORF_triplets_rm = remove_internal_stop(REF_triplets=REF_uORF_triplets,ALT_triplets=ALT_uORF_triplets,seq="ALT")
            REF_uORF_triplets = REF_uORF_triplets_rm
            ALT_uORF_triplets = ALT_uORF_triplets_rm
            if REF_uORF_triplets == ALT_uORF_triplets:
                var["DEBUG_analyze_substitutions"] += "stop codon readthrough variant,"
            REF_stop_positions = detect_triplet_positions(REF_uORF_triplets,codons="stop_codons")
            ALT_stop_positions = detect_triplet_positions(ALT_uORF_triplets,codons="stop_codons")

        if len(ALT_stop_positions) != 0 and len(REF_stop_positions) != 0:
            REF_stop = REF_uORF_triplets[REF_stop_positions[0]]
            ALT_stop = ALT_uORF_triplets[ALT_stop_positions[0]]
        
        # original uORF has no stop codon -> end is end of mRNA
        if len(REF_stop_positions) == 0:
            REF_stop = "NA"
            if var["REF_stop_located_in"] != "no_matching_stop_codon":
                var["DEBUG_analyze_substitutions"] += "No stop in REF_uORF found but not identified as 'no matching stop codon' by uORF finder,"
            
            if len(ALT_stop_positions) != 0:
                # case 1: ALT uORF has stop
                # new premature stop codon (previously no stop codon in REF uORF)
                var["uORF_stop_effect"] = "NA->uSTOPinSeq"
                var["ALT_stop_codon"] = ALT_uORF_triplets[ALT_stop_positions[0]]
                
                #changed 15.06.22
                ALT_uORF = ""
                for triplet_number in range(0,ALT_stop_positions[0]+1):
                    ALT_uORF += ALT_uORF_triplets[triplet_number]
                ALT_uORF_end_pos_mRNA = var["uORF_start_pos_mRNA"]+len(ALT_uORF)

                
            elif len(ALT_stop_positions) == 0:
                # case 2: ALT uORF has also no stop -> no effect
                var["uORF_stop_effect"] = "."
                var["ALT_stop_codon"] = "."
                ALT_uORF = ALT_uORFregion
                ALT_uORF_end_pos_mRNA = "."
        
        # original uORF has stop codon
        elif len(REF_stop_positions) != 0:
            # take last stop triplet in case of internal stop codons
            REF_stop = REF_uORF_triplets[REF_stop_positions[-1]]
            
            # case 1: there is at least one stop in ALT uORF
            if len(ALT_stop_positions) != 0:
                ALT_stop = ALT_uORF_triplets[ALT_stop_positions[0]]
                # premature stop codon because first stop is different
                if REF_stop_positions[0] != ALT_stop_positions[0]:
                    var["uORF_stop_effect"] = "uSTOP->uSTOPinSeq"
                    var["ALT_stop_codon"] = ALT_uORF_triplets[ALT_stop_positions[0]]
                    ALT_uORF = ""
                    for triplet_number in range(0,ALT_stop_positions[0]+1):
                        ALT_uORF += ALT_uORF_triplets[triplet_number]
                    ALT_uORF_end_pos_mRNA = var["uORF_start_pos_mRNA"]+len(ALT_uORF)
                    
                # different stop codon but at same position
                elif REF_stop_positions[0] == ALT_stop_positions[0] and REF_stop != ALT_stop:
                    var["uORF_stop_effect"] = "uSTOP->uSTOP"
                    var["ALT_stop_codon"] = ALT_uORF_triplets[ALT_stop_positions[0]]
                    ALT_uORF = ALT_uORFregion
                    ALT_uORF_end_pos_mRNA = "."
                    
                # neither position nor codon of stop changed
                elif REF_stop_positions[0] == ALT_stop_positions[0] and REF_stop == ALT_stop:
                    var["uORF_stop_effect"] = "."
                    var["ALT_stop_codon"] = "."
                    if "pos" not in var["uORF_start_effect"]:
                        ALT_uORF = ALT_uORFregion
                    ALT_uORF_end_pos_mRNA = "."
            
            # case 2: there is no stop codon in ALT uORF -> stop loss, search downstream for new stop codon    
            elif len(ALT_stop_positions) == 0:
                # extract downstream mRNA sequence to uORF sequence and format to triplets
                downstream_seq = var["mRNA_SEQ"][int(var["uORF_end_pos_mRNA"]):]
                downstream_triplets = translate_to_list(downstream_seq,side="right")
                ALT_stop_positions = detect_triplet_positions(downstream_triplets,codons="stop_codons")
                
                if len(ALT_stop_positions) == 0:
                    var["uORF_stop_effect"] = "uSTOP->loss"
                    var["ALT_stop_codon"] = "NA"
                    ALT_uORF = ALT_uORFregion+downstream_seq
                    ALT_uORF_end_pos_mRNA = var["uORF_end_pos_mRNA"]+len(downstream_seq)
                    
                else:
                    var["uORF_stop_effect"] = "uSTOP->uSTOPdownstream"
                    var["ALT_stop_codon"] = downstream_triplets[ALT_stop_positions[0]]
                    # append new sequence to uORF
                    appended_seq = ""
                    for triplet_number in range(0,ALT_stop_positions[0]+1):
                        appended_seq += downstream_triplets[triplet_number]
                    ALT_uORF = ALT_uORFregion+appended_seq
                    ALT_uORF_end_pos_mRNA = var["uORF_end_pos_mRNA"]+len(appended_seq)
                    
                    if len(ALT_uORF) != (ALT_uORF_end_pos_mRNA - var["uORF_start_pos_mRNA"]):
                        var["DEBUG_analyze_substitutions"] += "len(ALT_uORF) != (ALT_uORF_end_pos_mRNA - uORF_start_pos_mRNA)"
                        
    #-----------------------------------------------
    # uORF SEQUENCE and LENGTH
    #-----------------------------------------------
    if var["ALT_start_codon"] == "loss":
        var["uORF_seq_effect"] = "."
        ALT_uORF_length = "."
        
    else:
        if ALT_uORF == REF_uORF:
            var["uORF_seq_effect"] = "."
            ALT_uORF_length = "."
        elif len(ALT_uORF) < len(REF_uORF):
            var["uORF_seq_effect"] = "shorter sequence"
            ALT_uORF_length = len(ALT_uORF)
        elif len(ALT_uORF) > len(REF_uORF):
            var["uORF_seq_effect"] = "longer sequence"
            ALT_uORF_length = len(ALT_uORF)
        elif len(ALT_uORF) == len(REF_uORF) and ALT_uORF != REF_uORF:
            var["uORF_seq_effect"] = "altered sequence"
            ALT_uORF_length = "."

    #-----------------------------------------------
    # CDS distance
    #-----------------------------------------------
    if ALT_uORF_length == "." or ALT_uORF_end_pos_mRNA == ".":
        ALT_CDS_distance = "."
    else:
        ALT_CDS_distance = var["RcdsStart"]-ALT_uORF_end_pos_mRNA
        
    #-----------------------------------------------
    # uORF type effect, added on 28.06.2022
    #-----------------------------------------------
    if ALT_CDS_distance != "." and ALT_CDS_distance < 0 and var["REF_CDS_distance"] >= 0:
        var["uORF_type_effect"] = var["REF_uORF_type"]+"->oORF"
    elif ALT_CDS_distance != "." and ALT_CDS_distance >= 0 and var["REF_CDS_distance"] < 0:
        var["uORF_type_effect"] = var["REF_uORF_type"]+"->uORF"
    else:
        var["uORF_type_effect"] = "."
        
    #-----------------------------------------------
    # stop located in
    #-----------------------------------------------
    if ALT_uORF_end_pos_mRNA == ".":
        ALT_stop_located_in = "."
    elif ALT_uORF_end_pos_mRNA == "uORF loss":
        ALT_stop_located_in = "uORF loss"
    elif var["ALT_stop_codon"] == "NA":
        ALT_stop_located_in = "no matching stop codon"
    else:
        ALT_stop_located_in = stop_location(stop_pos=ALT_uORF_end_pos_mRNA,R5utrStart=var["R5utrStart"],R5utrEnd=var["R5utrEnd"],
                                            RcdsStart=var["RcdsStart"],RcdsEnd=var["RcdsEnd"],RtxEnd=var["RtxEnd"])
        
    #-----------------------------------------------
    # OUTPUT
    #----------------------------------------------- 
    if ALT_uORF == REF_uORF:
        var["ALT_uORF"] = "."
    else:
        var["ALT_uORF"] = ALT_uORF
        
    var["ALT_CDS_distance"] = ALT_CDS_distance
    var["ALT_uORF_length"] = ALT_uORF_length
    var["ALT_uORF_end_pos_mRNA"] = ALT_CDS_distance
    var["ALT_stop_located_in"] = ALT_stop_located_in
    if var["REF_stop_located_in"] == var["ALT_stop_located_in"]:
        var["ALT_stop_located_in"] = "."
    
    if var["uORF_stop_effect"] == ".":
        var["ALT_stop_codon"] = "."
        
    ###DEBUGGING
    var["#fasta_sequences"] = ">REF_SEQ\n"
    var["#fasta_sequences"] += var["REF_kozak_fragment"]+"NNN"+var["REF_uORF"]+"\n"
    var["#fasta_sequences"] += ">VAR_SEQ\n"
    var["#fasta_sequences"] += var["ALT_kozak_fragmentregion"]+"NNN"+var["ALT_uORFregion"]+"\n"
    var["#fasta_sequences"] += ">ALT_SEQ\n"
    var["#fasta_sequences"] += ALT_kozak_context_int[:6]+"NNN"+ALT_uORF+"\n"
    
    #-----------------------------------------------
    # HGVS https://varnomen.hgvs.org/bg-material/simple/
    #----------------------------------------------- 
    if var["mutation_type"] == "point_mutation":
        var["HGVS"] = var["CHROM"]+":g."+str(var["POS"])+var["REF"]+">"+var["ALT"]
    elif var["mutation_type"] == "multiple_bp_substitution":
        #https://www.ncbi.nlm.nih.gov/snp/rs121913227#hgvs_tab
        if str(Seq(var["ALT"]).reverse_complement()) == var["REF"]:
            var["HGVS"] = var["CHROM"]+":g."+str(var["POS"])+"_"+str(var["POS"]+(len(var["ALT"])-1))+"inv"
        else:
            var["HGVS"] = var["CHROM"]+":g."+str(var["POS"])+"_"+str(var["POS"]+(len(var["ALT"])-1))+"delins"+var["ALT"]

    return var


def analyze_indels(var:pd.Series):

    var["uORF_start_effect"] = "<>"
    var["uORF_stop_effect"] = "<>"
    var["uORF_kozak_effect"] = "<>"
    var["uORF_seq_effect"] = "<>"
    var["ALT_start_codon"] = "<>"
    var["ALT_stop_codon"] = "<>"
    var["ALT_kozak_context"] = "<>"
    var["ALT_kozak_strength"] = "<>"
    var["ALT_stop_located_in"] = "<>"
    var["ALT_uORF_length"] = "<>"
    var["ALT_CDS_distance"] = "<>"
    var["ALT_uORF"] = "<>"
    var["HGVS"] = "<>"

    var["#fasta_sequences"] = ">REF_SEQ\n"
    var["#fasta_sequences"] += var["REF_kozak_fragment"]+"NNN"+var["REF_uORF"]+"\n"
    var["#fasta_sequences"] += ">VAR_SEQ\n"
    var["#fasta_sequences"] += var["ALT_kozak_fragmentregion"]+"NNN"+var["ALT_uORFregion"]+"\n"
    
    return var


def merge_uORF_variants(variant_dataframe, identifier_columns, merge_columns):
    
    variants = variant_dataframe.to_dict("index")
    
    unique_variants = {}
    
    for index in variants:
        
        ID = ""
        for column in identifier_columns:
            if ID != "":
                ID += "_"
            ID += str(variants[index][column])
        
        if ID in unique_variants:
            for column in variants[index]:
                if column in merge_columns:
                    unique_variants[ID][column].append(variants[index][column])
                    
        elif ID not in unique_variants:
            unique_variants[ID] = {}
            for column in variants[index]:
                if column in merge_columns:
                    unique_variants[ID][column] = [(variants[index][column])]
                else:
                    unique_variants[ID][column] = (variants[index][column])

    merged_variants = pd.DataFrame.from_dict(unique_variants, 'index')
    
    return merged_variants

