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
# Functions used for Variant Analysis scripts
#
# ==========================================================================================

import pandas as pd
from Bio.Seq import Seq

def string_to_list(string):
    my_list = list(string.split(","))
    my_list.pop()
    my_list = list(map(int, my_list))
    return(my_list)

def calculate_relative_pos(exonStarts,exonEnds,strand,pos):
    
    relative_position_plus = 0
    position_dict_plus = {}
    
    for i in range(len(exonStarts)):
        
        for j in range(exonStarts[i],exonEnds[i]):

            position_dict_plus[j] = relative_position_plus
            relative_position_plus += 1
        
    if pos not in position_dict_plus:
        if pos in range(exonStarts[0],exonEnds[-1]):
            relative_position = "intron"
        elif pos < exonStarts[0] or pos >= exonEnds[-1]:
            relative_position = "out_of_range"
        else:
            relative_position = "out_of_range"
            
    else:
        if strand == "-":
            relative_position = (len(position_dict_plus)-1) - position_dict_plus[pos]
        elif strand == "+": 
            relative_position = position_dict_plus[pos]
    
    return(relative_position)


def add_uORF_IDs(row:pd.Series, uORF_data):
    
    row["DEBUG_add_uORF_IDs"] = ""
    row["variant_location"] = "."
    row["uORF_uniqueID"] = []
    out_of_range_IDs = []
    intron_IDs = []
    
    CHROM = row["CHROM"]
    FROM = int(row["FROM_0BHO"])
    TO = int(row["TO_0BHO"])

    matching_uORFs = uORF_data[(uORF_data["chrom_id"] == CHROM) & 
                               ((uORF_data["start_incl_kozak"] <= FROM) & (uORF_data["stop_incl_kozak"] > FROM) |
                                (uORF_data["start_incl_kozak"] < TO) & (uORF_data["stop_incl_kozak"] >= TO))
                               ]

    if len(matching_uORFs) == 0:
        row["DEBUG_add_uORF_IDs"] += "no uORFs match this variant,"
    
    for index, uORF in matching_uORFs.iterrows():

        uORF["uORF_exonStarts"]
        
        uORF_exonStarts = string_to_list(uORF["uORF_exonStarts"])
        uORF_exonEnds = string_to_list(uORF["uORF_exonEnds"])
    
        # checking if in intron or out of range, strand is not relevant yet
        FROM_in_uORF = calculate_relative_pos(exonStarts=uORF_exonStarts, exonEnds=uORF_exonEnds, strand="+", pos=FROM)
        TO_in_uORF = calculate_relative_pos(exonStarts=uORF_exonStarts, exonEnds=uORF_exonEnds, strand="+", pos=TO-1)
    
        if FROM_in_uORF != "intron" and TO_in_uORF != "intron" and FROM_in_uORF != "out_of_range" and TO_in_uORF != "out_of_range":
            row["uORF_uniqueID"].append(uORF["uORF_uniqueID"])
        elif (FROM_in_uORF == "out_of_range" and TO_in_uORF != "out_of_range") or (FROM_in_uORF != "out_of_range" and TO_in_uORF == "out_of_range"):
            out_of_range_IDs.append(uORF["uORF_uniqueID"])
        elif (FROM_in_uORF != "intron" and TO_in_uORF == "intron") or (FROM_in_uORF == "intron" and TO_in_uORF != "intron"):
            intron_IDs.append(uORF["uORF_uniqueID"])
        elif FROM_in_uORF == "intron" and TO_in_uORF == "intron":
            intron_IDs.append(uORF["uORF_uniqueID"])

    if len(row["uORF_uniqueID"]) == 0 and len(intron_IDs) != 0:
        row["DEBUG_add_uORF_IDs"] += "partly_in_intron,"
        row["variant_location"] = "partly_not_in_uORF"
        row["uORF_uniqueID"] = intron_IDs
        
    elif len(row["uORF_uniqueID"]) == 0 and len(out_of_range_IDs) != 0:
        row["DEBUG_add_uORF_IDs"] += "partly_out_of_range"
        row["variant_location"] = "partly_not_in_uORF"
        row["uORF_uniqueID"] = out_of_range_IDs
        
    elif len(row["uORF_uniqueID"]) == 0 and len(out_of_range_IDs) == 0 and len(intron_IDs) == 0:
        row["DEBUG_add_uORF_IDs"] += "all matching uORFs are completely out of range"
        row["variant_location"] = "partly_not_in_uORF"
        row["uORF_uniqueID"] = out_of_range_IDs
        print(out_of_range_IDs)
        print("out of range!")
        print(CHROM)
        print(FROM)
        print(TO)
         
    return row


def mutation_type(ALT,REF):
        
    if len(REF) > len(ALT):
        mutation_type = "deletion"
    elif len(REF) < len(ALT):
        mutation_type = "insertion"
    elif len(REF) == len(ALT):
        if len(REF) > 1:
            mutation_type = "multiple_bp_substitution"
        if len(REF) == 1:
            mutation_type = "point_mutation"
            
    return mutation_type


def extract_exonic_seq(genomic_seq,strand,exonStarts,exonEnds,chrom):
    
    from Bio.Seq import Seq
    
    fasta_header_dict = {}
    for header in genomic_seq.keys():
        fasta_header_dict[header.split(" ")[0]] = header
    
    exonic_seq = ""
    for x in range(len(exonStarts)):
        exonic_seq += str(genomic_seq[fasta_header_dict[chrom]][exonStarts[x]:exonEnds[x]])
    exonic_seq = exonic_seq.upper()
    
    if strand == "-":
        exonic_seq = (Seq(exonic_seq)).reverse_complement()
        
    return str(exonic_seq)


def location_of_variant(cdsStart,cdsEnd,exonStarts,exonEnds,FROM,TO,sep,strand):
    
    if strand == "+":
        regions = {"TLS" : range(exonStarts[0],cdsStart),
                   "CDS" : range(cdsStart,cdsEnd),
                   "3'UTR" : range(cdsEnd,exonEnds[-1])}
    elif strand == "-":
        regions = {"TLS" : range(cdsEnd,exonEnds[-1]),
                   "CDS" : range(cdsStart,cdsEnd),
                   "3'UTR" : range(exonStarts[0],cdsStart)}
    
    affected_regions = ""
    for region in regions:
        if FROM in regions[region] or TO-1 in regions[region]:
            affected_regions += region+sep
            
    if affected_regions == "TLS"+sep+"3'UTR"+sep:
        affected_regions = "TLS"+sep+"CDS"+sep+"3'UTR"+sep
    
    return affected_regions


def correct_uORFdbID(ID):
    if ID.count(".") == 3 and ID.count("_") == 1:
        id_list = ID.split(".")
        new_id = id_list[0]+"."+id_list[1]+"_"+id_list[2]+"."+id_list[3]
    elif id.count(".") == 2 and id.count("_") == 2:
        new_id = ID
    return new_id


def add_uORF_data(var:pd.Series, uORF_data, genome_sequences):
    
    var["DEBUG_add_uORF_data"] = ""
    
    try:
    
        uORF = uORF_data[var["uORF_uniqueID"]]
        
    except:
        print(var)
        print(var["uORF_uniqueID"])
    
    if uORF["CDS_stop_info"] == "stop_codon_readthrough":
        var["DEBUG_add_uORF_data"] += "stop_codon_readthrough,"
    
    # copy columns, added uORF type and CDS in other transcripts on 28.06.2022
    copy_columns = ["uORFdb_ID","chrom","strand","gene_symbol","id","tv_count_for_gene","tv_count_with_same_uORF","TV_with_uORF/TV_total","kozak_status",
                    "uORF_type","RcdsStart","RcdsEnd","RtxEnd","R5utrStart","R5utrEnd","CDS_in_other_transcripts"]

    for column in copy_columns:
        var[column] = uORF[column]
    
    # only necessary if uORFdbID in wrong format: "NM_001242784.3.CTG.20" -> "NM_001242784.3_CTG.20"
    var["uORFdb_ID"] = correct_uORFdbID(uORF["uORFdb_ID"])

    var["mutation_type"] = mutation_type(ALT=var["ALT"],REF=var["REF"])
    
    var["REF_kozak_context"] = uORF["uORF_kozak_context"]
    var["REF_kozak_fragment"] = uORF["uORF_kozak_fragment"]
    var["REF_kozak_strength"] = uORF["uORF_kozak_strength"]
    var["REF_CDS_distance"] = uORF["cds_start_dist_no_introns"]
    var["REF_uORF_length"] = uORF["uorf_length_excl_introns"]
    var["REF_stop_located_in"] = uORF["stop_located_in"]
    if var["REF_kozak_fragment"] == "NA":  
        var["REF_kozak_fragment"] = ""
        
    var["REF_uORF_type"] = uORF["uORF_type"]
    
    #------------------------------------------------------------
    # VARIABLES
    #------------------------------------------------------------
    uORF_exonStarts = string_to_list(uORF["uORF_exonStarts"])
    uORF_exonEnds = string_to_list(uORF["uORF_exonEnds"])
    exonStarts = string_to_list(uORF["exonStarts"])
    exonEnds = string_to_list(uORF["exonEnds"])
    strand = uORF["strand"]
    FROM = var["FROM_0BHO"]
    TO = var["TO_0BHO"]
    var["uORF_start_pos_mRNA"] = uORF["start_pos_mRNA"]
    var["uORF_end_pos_mRNA"] = uORF["end_pos_mRNA"]

    #------------------------------------------------------------
    # get correct uORF start and stop COORDINATES and REF and ALT
    #------------------------------------------------------------
    if strand == "+":
        FROM_in_uORF = calculate_relative_pos(exonStarts=uORF_exonStarts, exonEnds=uORF_exonEnds, strand=strand, pos=FROM)
        TO_in_uORF = calculate_relative_pos(exonStarts=uORF_exonStarts, exonEnds=uORF_exonEnds, strand=strand, pos=TO-1)+1
        FROM_in_mRNA = calculate_relative_pos(exonStarts=exonStarts, exonEnds=exonEnds, strand=strand, pos=FROM)
        TO_in_mRNA = calculate_relative_pos(exonStarts=exonStarts, exonEnds=exonEnds, strand=strand, pos=TO-1)+1
        REF = var["REF"]
        ALT = var["ALT"]
        
    elif strand == "-":
        FROM_in_uORF = calculate_relative_pos(exonStarts=uORF_exonStarts, exonEnds=uORF_exonEnds, strand=strand, pos=TO-1)
        TO_in_uORF = calculate_relative_pos(exonStarts=uORF_exonStarts, exonEnds=uORF_exonEnds, strand=strand, pos=FROM)+1
        FROM_in_mRNA = calculate_relative_pos(exonStarts=exonStarts, exonEnds=exonEnds, strand=strand, pos=TO-1)
        TO_in_mRNA = calculate_relative_pos(exonStarts=exonStarts, exonEnds=exonEnds, strand=strand, pos=FROM)+1
        REF = str(Seq(var["REF"]).reverse_complement())
        ALT = str(Seq(var["ALT"]).reverse_complement())
    
    #-----------------------------------------------
    # get sequence, SEQ refers to GRCh38 sequence!!!
    #-----------------------------------------------
    var["mRNA_SEQ"] = extract_exonic_seq(genomic_seq=genome_sequences,
                                         strand=strand,
                                         exonStarts=exonStarts,
                                         exonEnds=exonEnds,
                                         chrom=uORF["chrom_id"])
    
    var["VAR_mRNA_SEQ"] = var["mRNA_SEQ"][:FROM_in_mRNA] + ALT + var["mRNA_SEQ"][TO_in_mRNA:]

    SEQ = extract_exonic_seq(genomic_seq=genome_sequences,
                             strand=strand,
                             exonStarts=uORF_exonStarts,
                             exonEnds=uORF_exonEnds,
                             chrom=uORF["chrom_id"])
    
    #DEBUG: double-check sequences
    if uORF["uORF_kozak_fragment"] == "NA":
        uORF_kozak_fragment_debug = ""
    else:
        uORF_kozak_fragment_debug = uORF["uORF_kozak_fragment"]
    if SEQ != uORF_kozak_fragment_debug+uORF["seq_no_introns"]:
        var["DEBUG_add_uORF_data"] += "uORF sequence could not be restored with uORF coordinates and genomic fasta,\n"
    if uORF["mRNA_seq"] != var["mRNA_SEQ"]:
        var["DEBUG_add_uORF_data"] += "mRNA_seq from uORF finder and extracted mRNA by this function do not match,"
        
    #-----------------------------
    # compare REF and GRCh38
    #-----------------------------
    try:
        GRCh38_REF = SEQ[FROM_in_uORF:TO_in_uORF]
    except:
        print(uORF)
        print(SEQ)
        print(FROM_in_uORF)
        print(TO_in_uORF)

    if REF != GRCh38_REF:
        additional_germline = "additional_germline"
        var["DEBUG_add_uORF_data"] += "REF does not match GRCh38 but may be germline mutation,"
        
        if GRCh38_REF == ALT:
            additional_germline = "back-mutation"

    elif REF == GRCh38_REF:
        additional_germline = "only_somatic"
        
    #-----------------------------
    # build REF and ALT sequences
    #-----------------------------
    start_pos = len(var["REF_kozak_fragment"])
    
    # sequences incl kozak context, no matter if kozak is complete (fragment is extracted)
    REF_SEQ = SEQ[:FROM_in_uORF] + REF + SEQ[TO_in_uORF:]
    ALT_SEQ = REF_SEQ[:FROM_in_uORF] + ALT + REF_SEQ[TO_in_uORF:]

    REF_uORF = REF_SEQ[start_pos:]
    REF_kozak = REF_SEQ[:start_pos]
    
    if start_pos == 0:
        REF_kozak = ""
    
    # variant is in kozak context
    if FROM_in_uORF < start_pos and TO_in_uORF <= start_pos:
        ALT_kozak = REF_SEQ[:FROM_in_uORF] + ALT + REF_SEQ[TO_in_uORF:start_pos]
        ALT_uORF = REF_uORF
    
    # variant is in uORF sequence
    elif FROM_in_uORF >= start_pos and TO_in_uORF > start_pos:
        ALT_kozak = REF_kozak
        ALT_uORF = REF_SEQ[start_pos:FROM_in_uORF] + ALT + REF_SEQ[TO_in_uORF:] 
    
    # variant is in kozak and uORF
    elif FROM_in_uORF < start_pos and TO_in_uORF > start_pos:
        if var["mutation_type"] == "deletion":
            # respect + and - because - is right-aligned
            if strand == "+":
                ALT_kozak = REF_SEQ[:FROM_in_uORF+1]
                ALT_uORF = REF_SEQ[TO_in_uORF:]
                    
            elif strand == "-":
                ALT_kozak = REF_SEQ[:FROM_in_uORF]
                ALT_uORF = REF_SEQ[TO_in_uORF-1:]
                    
        elif var["mutation_type"] == "multiple_bp_substitution":
            ALT_kozak = ALT_SEQ[:start_pos]
            ALT_uORF = ALT_SEQ[start_pos:]
            
        else:
            var["DEBUG_add_uORF_data"] += "mutation affects both kozak and uORF but mutation type is not deletion or multiple_bp_substitution,"

    if ALT_kozak.replace(".", "")+ALT_uORF != ALT_SEQ:
        var["DEBUG_add_uORF_data"] += "ALT_kozak+ALT_uORF != ALT_SEQ,"

    # check which mRNA regions are affected
    var["variant_location_mRNA"] = location_of_variant(cdsStart=uORF["cdsStart"],cdsEnd=uORF["cdsEnd"],exonStarts=exonStarts,exonEnds=exonEnds,FROM=FROM,TO=TO,sep=",",strand=var["strand"])
    
    if "CDS" in var["variant_location_mRNA"]:
        var["cds_info"] = "affected"
    else:
        var["cds_info"] = "."
        
    if var["cds_info"] == "affected" and uORF["uORF_type"] == "non-overlapping":
        var["DEBUG_add_uORF_data"] += "Variant affects CDS but uORF is non-overlapping,"
    
    #-----------------------------
    # OUTPUT
    #-----------------------------
    if strand == "-":
        var["GRCh38_REF"] = Seq(GRCh38_REF).reverse_complement()
    elif strand == "+":
        var["GRCh38_REF"] = GRCh38_REF
        
    var["comparison_GRCh38/REF/ALT"] = additional_germline
    
    var["REF_start_codon"] = REF_uORF[:3]
    var["REF_stop_codon"] = REF_SEQ[-3:]

    var["FROM_in_uORF"] = FROM_in_uORF
    var["TO_in_uORF"] = TO_in_uORF
    var["FROM_in_mRNA"] = FROM_in_mRNA
    var["TO_in_mRNA"] = TO_in_mRNA
    
    var["REF_SEQ"] = REF_SEQ
    var["ALT_SEQregion"] = ALT_SEQ
    
    var["REF_uORF"] = REF_uORF
    var["ALT_uORFregion"] = ALT_uORF#region where uORF was
    var["REF_kozak_fragment"] = REF_kozak
    var["ALT_kozak_fragmentregion"] = ALT_kozak#region where kozak was

    var["variant-uorf-id"] = var["CHROM"]+"_"+str(var["POS"])+"_"+var["REF"]+"_"+var["ALT"]+"-"+str(var["uORF_uniqueID"])
    
    return var


