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
# Panda functions used for uORF finder
#
# ==========================================================================================

# ---------------------------------------------------------------------------------
# FUNCTIONS
# ---------------------------------------------------------------------------------
import pandas as pd
from Bio.Seq import Seq
import hashlib
from BASIC_FUNCTIONS import *
import math

def add_tv_counts(row: pd.Series, dataframe) -> pd.Series:
    
    if row["selection_status"] == "SELECTED" or row["selection_status"] == "no_5'UTR,":
        tv_count = dataframe[(dataframe.selection_status == "SELECTED") | (dataframe.selection_status == "no_5'UTR,")]["name2"].values.tolist().count(row["name2"])
        row["tv_count_for_gene"] = tv_count
    else:
        row["tv_count_for_gene"] = "."
    
    return row

def filter_RefSeq(row: pd.Series, chromosomes) -> pd.Series:
    
    # add some data
    if row["chrom"] in chromosomes:
        row["fasta_header"] = chromosomes[row["chrom"]]
    else:
        row["fasta_header"] = "?"

    row["selection_status"] = ""
    
    # preselect transcripts
    if "NR" in row["name"] or "YP" in row["name"]:
        row["selection_status"] += "invalid_id,"
    
    # changed on 19.04.22:
    if "alt" in row["chrom"] or "fix" in row["chrom"] or "random" in row["chrom"] or "chrUn" in row["chrom"]:
        row["selection_status"] += "invalid_chr,"
        
    if type(row["name2"]) != str:
        if math.isnan(row["name2"]) == True:
            print("NAN still found")
            row["selection_status"] += "invalid_gene,"
    
    if row["selection_status"] == "":
        row["duplicate_candidate"] = "duplicate_candidate"+row["name"]
    else:
        row["duplicate_candidate"] = "."
        
    if row["cdsStart"] == row["cdsEnd"]:
        row["selection_status"] += "no_CDS,"

    exonStarts = string_to_list(row["exonStarts"])
    exonEnds = string_to_list(row["exonEnds"])
                
    RtxStart = calculate_relative_coordinates(exonStarts, exonEnds, row["strand"], row["txStart"])
    RtxEnd = calculate_relative_coordinates(exonStarts, exonEnds, row["strand"], row["txEnd"])
            
    RcdsStart = calculate_relative_coordinates(exonStarts, exonEnds, row["strand"], row["cdsStart"])
    RcdsEnd = calculate_relative_coordinates(exonStarts, exonEnds, row["strand"], row["cdsEnd"])
            
    RexonStarts = []
    for exonStart in exonStarts:
        RexonStart = calculate_relative_coordinates(exonStarts, exonEnds, row["strand"], exonStart)
        RexonStarts.append(RexonStart)
    RexonStarts.sort()
    RexonStarts = list_to_string(RexonStarts)
            
    RexonEnds = []
    for exonEnd in exonEnds:
        RexonEnd = calculate_relative_coordinates(exonStarts, exonEnds, row["strand"], exonEnd)
        RexonEnds.append(RexonEnd)
    RexonEnds.sort()
    RexonEnds = list_to_string(RexonEnds)
            
    if row["strand"] == "-":#swap start/end for minus strand
        RtxStart, RtxEnd = RtxEnd, RtxStart
        RcdsEnd, RcdsStart = RcdsStart, RcdsEnd
        RexonEnds, RexonStarts = RexonStarts, RexonEnds

    R5utrStart = RtxStart
    R5utrEnd = RcdsStart
    R5utrLength = R5utrEnd - R5utrStart

    row["RtxStart"] = RtxStart
    row["RtxEnd"] = RtxEnd
    row["RcdsStart"] = RcdsStart
    row["RcdsEnd"] = RcdsEnd
    row["RexonStarts"] = RexonStarts
    row["RexonEnds"] = RexonEnds
    row["R5utrStart"] = R5utrStart
    row["R5utrEnd"] = R5utrEnd
    row["R5utrLength"] = R5utrLength

    if R5utrLength == 0:
        row["selection_status"] += "no_5'UTR,"
        
    if row["selection_status"] == "":
        row["selection_status"] = "SELECTED"

    return row

def extract_sequences(row: pd.Series, genomic_fasta) -> pd.Series:

    from pyfaidx import Fasta

    if "invalid" not in row["selection_status"] and "no_CDS" not in row["selection_status"] and row["fasta_header"] != "?":
        
        genomic_seq = Fasta(genomic_fasta)
        
        #format strings from RefSeq to iterable lists
        exonStarts = string_to_list(row["exonStarts"])
        exonEnds = string_to_list(row["exonEnds"])
        
        mRNA_seq = ""
        
        fasta_header = row["fasta_header"]
        
        #calculate transcript seq without introns
        for i in range(0, len(exonStarts)):
            if genomic_seq[fasta_header][exonStarts[i]:exonEnds[i]] == None:
                mRNA_seq = None
                print(row["name"])
                print("sequence extraction not possible")
                row["selection_status"] += "sequence extraction not possible,"
                #break
            else:
                mRNA_seq += str(genomic_seq[fasta_header][exonStarts[i]:exonEnds[i]])
        
        if row["strand"] == "-":
            mRNA_seq = (Seq(mRNA_seq)).reverse_complement()
        
        #check if extracted correctly
        if mRNA_seq == None or mRNA_seq == "":
            row["selection_status"] = "no_seq,"
            row["mRNA_seq"] = "."
            
        else:
            
            mRNA_seq = str(mRNA_seq).upper()

            RcdsStart = row["RcdsStart"]
            if RcdsStart >= 6:
                CDS_kozak_context = mRNA_seq[RcdsStart-6:RcdsStart+4]
                CDS_kozak_strength = determine_kozak_strength(CDS_kozak_context)
            else:
                CDS_kozak_fragment = mRNA_seq[0:RcdsStart+4]
                CDS_kozak_context = (10-len(CDS_kozak_fragment))*"-"+CDS_kozak_fragment
                if CDS_kozak_context[3] != "-":
                    CDS_kozak_strength = determine_kozak_strength(CDS_kozak_context)
                else:
                    CDS_kozak_strength = "."
            
            row["mRNA_seq"] = mRNA_seq
            row["CDS_kozak_context"] = CDS_kozak_context
            row["CDS_kozak_strength"] = CDS_kozak_strength
    
    else:
        row["mRNA_seq"] = "."
        row["CDS_kozak_context"] = "."
        row["CDS_kozak_strength"] = "."
    
    return row

def detect_start_pos(row : pd.Series):
    
    start_triplets = ["ATG", "AAG", "ATT", "TTG", "AGG", "ATC", "GTG", "ACG", "CTG", "ATA"]
    
    seq = row["mRNA_seq"]

    R5utrEnd = int(row["R5utrEnd"])

    start_positions = get_codon_positions(start_triplets,seq[0:R5utrEnd])

    row["start_pos_mRNA"] = start_positions
    
    if len(start_positions) == 0:
        row["start_pos_mRNA"] = "."
    
    return row

def detect_uORFs(row : pd.Series, minus_transcripts_dict, plus_transcripts_dict):
    
    stop_triplets = ["TAG", "TGA", "TAA"]
    seq = row["mRNA_seq"]
    ID = row["name"]
    start_pos = row["start_pos_mRNA"]
    reading_frame = start_pos % 3
    exonStarts = string_to_list(row["exonStarts"])
    exonEnds = string_to_list(row["exonEnds"])
    RexonStarts = string_to_list(row["RexonStarts"])
    RexonEnds = string_to_list(row["RexonEnds"])
    RcdsStart = row["RcdsStart"]
    RcdsEnd = row["RcdsEnd"]
    R5utrEnd = row["R5utrEnd"]
    stop_positions = get_codon_positions(stop_triplets,seq)
    
    matching_stop_positions = [number for number in stop_positions if number > start_pos and number % 3 == reading_frame]
    
    # DEFINE STOP CODON
    
    #check if matching stop codons are found    
    if len(matching_stop_positions) != 0:
        stop_pos = min(matching_stop_positions)
        stop_codon = seq[stop_pos:stop_pos+3]
    else:
        stop_pos = (row["RtxEnd"]-3)
        stop_codon = "NA"

    # reading frame in relation to CDS (CDS start = reading frame 1), added on 23.03.22
    reading_frame_cds = ((RcdsStart - start_pos) % 3) + 1

    # 5'terminal extensions of CDS (overlapping and in-frame)
    if reading_frame_cds == 1 and stop_pos >= RcdsStart:

        # normal 5'terminal extension, CDSstop == uORFstop
        if stop_pos == (RcdsEnd-3):
            CDS_stop_info = "shared_stop"
        
        # internal stop codon + CDS normal stop -> set uORF stop to CDS stop
        elif stop_pos < (RcdsEnd-3) and seq[(RcdsEnd-3):RcdsEnd] in stop_triplets:
            CDS_stop_info = "stop_codon_readthrough,set_to_shared_stop"
            stop_pos = (RcdsEnd-3)
            stop_codon = seq[stop_pos:stop_pos+3]

        # internal stop codon + CDS strange stop -> search downstream for uORF stop
        elif stop_pos < (RcdsEnd-3) and seq[(RcdsEnd-3):RcdsEnd] not in stop_triplets:
            CDS_stop_info = "stop_codon_readthrough,invalid_CDS_stop_codon"
            # search for next stop codon after CDSend
            for pos in sorted(matching_stop_positions):
                if pos > (RcdsEnd-3):
                    stop_pos = pos
                    stop_codon = seq[stop_pos:stop_pos+3]
                    break
                else:
                    stop_pos = (row["RtxEnd"]-3)
                    stop_codon = "NA"

        # uORF stop was detected after CDS end because CDS has strange stop
        elif stop_pos >= RcdsEnd and seq[(RcdsEnd-3):RcdsEnd] not in stop_triplets:
            CDS_stop_info = "invalid_CDS_stop_codon"

        else:
            CDS_stop_info = "???"

    # all other uORFs
    else:
        CDS_stop_info = "."
    #---------------------------------------------------------------------------------

    #calculate chromosome coordinates
    Astart_coordinate = calculate_absolute_coordinates(exonStarts,exonEnds,RexonStarts,RexonEnds,row["strand"],start_pos)

    if row["strand"] == "-":
        Astop_coordinate = (calculate_absolute_coordinates(exonStarts,exonEnds,RexonStarts,RexonEnds,row["strand"],(stop_pos+2)))-1   
    elif row["strand"] == "+":
        Astop_coordinate = (calculate_absolute_coordinates(exonStarts,exonEnds,RexonStarts,RexonEnds,row["strand"],(stop_pos+2)))+1 
    
    # --------------------------------------------------------------------------------------
    # +/- strand not relevant yet as mRNA coordinates are used
    # --------------------------------------------------------------------------------------
    if stop_codon == "NA":
        stop_located_in = "no_matching_stop_codon"
    else:
        stop_located_in = stop_location(stop_pos,row["R5utrStart"],row["R5utrEnd"],row["RcdsStart"],row["RcdsEnd"],row["RtxEnd"])


    if start_pos >=6:
        kozak_status = "complete"
        Rstart_incl_kozak = start_pos-6
        kozak_context = seq[start_pos-6:start_pos+4]
        uORF_kozak_strength = determine_kozak_strength(kozak_context)
        row["uORF_kozak_fragment"] = seq[start_pos-6:start_pos]
    else:
        uORF_kozak_fragment = seq[0:start_pos+4]
        row["uORF_kozak_fragment"] = seq[0:start_pos]
        Rstart_incl_kozak = 0
        kozak_context = (10-len(uORF_kozak_fragment))*"-"+uORF_kozak_fragment
        if kozak_context[3] != "-":
            kozak_status = "incomplete"
            uORF_kozak_strength = determine_kozak_strength(kozak_context)
        else:
            kozak_status = "."
            uORF_kozak_strength = "."
        
    uorf_length_incl_introns = abs(Astop_coordinate-Astart_coordinate)
    uorf_length_excl_introns = ((stop_pos+3)-start_pos)

    cds_start_dist_no_introns = row["RcdsStart"]-(stop_pos+3)
    #cds_start_dist_no_introns = (stop_pos+3)-row["RcdsStart"]
    end_pos_mRNA = (stop_pos+3)
    seq_no_introns = seq[start_pos:stop_pos+3]

    # --------------------------------------------------------------------------------------
    # respect to +/- strand! reversed if necessary as chromosome coordinates are involved!!!
    # --------------------------------------------------------------------------------------
    #cds start distance with introns
    if row["strand"] == "-":#cdsEnd is cdsStart on minus strand
        Astop_coordinate, Astart_coordinate = Astart_coordinate, Astop_coordinate
        cds_start_dist = Astart_coordinate - row["cdsEnd"]
        #kozak start position
        Astart_incl_kozak = Astart_coordinate
        Astop_incl_kozak = calculate_absolute_coordinates(exonStarts,exonEnds,RexonStarts,RexonEnds,row["strand"],Rstart_incl_kozak)
        #switch start and stop for chromosome coordinates for minus-strand
        
    elif row["strand"] == "+":
        cds_start_dist = row["cdsStart"] - Astop_coordinate
        Astart_incl_kozak = calculate_absolute_coordinates(exonStarts,exonEnds,RexonStarts,RexonEnds,row["strand"],Rstart_incl_kozak)
        Astop_incl_kozak = Astop_coordinate
        
    if row["strand"] == "+":
        if Astart_coordinate in plus_transcripts_dict:
            CDS_in_other_transcripts = ""
            for transcript in plus_transcripts_dict[Astart_coordinate]:
                transcript_id = transcript.split(":")[0]
                gene_symbol = transcript.split(":")[1]
                if gene_symbol == row["name2"]:
                    CDS_in_other_transcripts += transcript_id+","
        else:
            CDS_in_other_transcripts = "."
        
    elif row["strand"] == "-":
        if Astop_coordinate in minus_transcripts_dict:
            CDS_in_other_transcripts = ""
            for transcript in minus_transcripts_dict[Astop_coordinate]:
                transcript_id = transcript.split(":")[0]
                gene_symbol = transcript.split(":")[1]
                if gene_symbol == row["name2"]:
                    CDS_in_other_transcripts += transcript_id+","
        else:
            CDS_in_other_transcripts = "."
    
    # unique identifiers, human-readable and checksum
    seq_checksum = hashlib.sha256((seq_no_introns).encode('utf-8')).hexdigest()
    uORF_uniqueID = row["chrom"]+"_"+ID+"_"+str(Astart_coordinate)+"_"+str(Astop_coordinate)+"_"+seq_checksum
    uORF_ID = row["chrom"]+"_"+str(Astart_coordinate)+"_"+str(Astop_coordinate)+"_"+seq_checksum
    
    modulo = len(seq_no_introns) % 3
    if modulo != 0:
        seq_trimmed = seq_no_introns[:-modulo]
        translated_uORF = Seq(seq_trimmed).translate() + "?"
    else:
        translated_uORF = Seq(seq_no_introns).translate()
    M_translated_uORF = "M" + str(translated_uORF)[1:]
    
    if cds_start_dist_no_introns >= 0:
        uORF_type = "non-overlapping"
    elif cds_start_dist_no_introns < 0 and reading_frame_cds == 1:
        if (Astop_coordinate == row["cdsEnd"] and row["strand"] == "+") or (Astart_coordinate == row["cdsStart"] and row["strand"] == "-"):
            uORF_type = "N-terminal_extension"
        else:
            #print("check shared stop calculation")
            uORF_type = "possible_N-terminal_extension"
    elif cds_start_dist_no_introns < 0 and reading_frame_cds != 1:
        if row["R5utrLength"] < 3:
            uORF_type = "overlapping_start"
        elif row["R5utrLength"] >= 3:
            uORF_type = "overlapping"
        
    start_codon = seq[start_pos:start_pos+3]

    start_positions_for_this_start_codon = get_codon_positions([start_codon],seq[0:R5utrEnd])
    pos_dict = {}
    for i in range(len(start_positions_for_this_start_codon)):
        pos_dict[start_positions_for_this_start_codon[i]] = i+1
    uORFdb_ID = ID+"_"+start_codon+"."+str(pos_dict[start_pos])


    uORF_exonStarts = []
    for exonStart in exonStarts:
        if exonStart > Astart_incl_kozak and exonStart < Astop_incl_kozak:
            uORF_exonStarts.append(exonStart)
    uORF_exonStarts.append(Astart_incl_kozak)
    uORF_exonStarts.sort()
    uORF_exonStarts=list_to_string(uORF_exonStarts)
    
    uORF_exonEnds = []
    for exonEnd in exonEnds:
        if exonEnd > Astart_incl_kozak and exonEnd < Astop_incl_kozak:
            uORF_exonEnds.append(exonEnd)
    uORF_exonEnds.append(Astop_incl_kozak)
    uORF_exonEnds.sort()
    uORF_exonEnds=list_to_string(uORF_exonEnds)

    row["chrom_id"] = row["fasta_header"].split(" ")[0]
    row["start"] = Astart_coordinate
    row["stop"] = Astop_coordinate
    row["id"] = ID
    row["uORF_ID"] = uORF_ID
    row["uORF_uniqueID"] = uORF_uniqueID
    row["uORFdb_ID"] = uORFdb_ID
    row["gene_symbol"] = row["name2"]
    
    row["start_codon"] = start_codon
    row["stop_codon"] = stop_codon
    
    row["end_pos_mRNA"] = end_pos_mRNA
    
    row["kozak_context"] = kozak_context
    
    row["uorf_length_incl_introns"] = uorf_length_incl_introns
    row["uorf_length_excl_introns"] = uorf_length_excl_introns
    
    row["cds_start_dist"] = cds_start_dist
    row["cds_start_dist_no_introns"] = cds_start_dist_no_introns
    
    row["stop_located_in"] = stop_located_in
    
    row["frame"] = reading_frame_cds
    
    row["uORF_kozak_context"] = kozak_context
    row["uORF_kozak_strength"] = uORF_kozak_strength
    
    row["seq_no_introns"] = seq_no_introns
    row["translated_uORF"] = M_translated_uORF #translated_uORF
    
    row["uORF_type"] = uORF_type
    row["CDS_in_other_transcripts"] = CDS_in_other_transcripts
    row["CDS_stop_info"] = CDS_stop_info
    
    row["start_incl_kozak"] = Astart_incl_kozak
    row["stop_incl_kozak"] = Astop_incl_kozak
    row["kozak_status"] = kozak_status
    
    row["uORF_exonStarts"] = uORF_exonStarts
    row["uORF_exonEnds"] = uORF_exonEnds

    return row


def add_transcript_and_gene_data(row: pd.Series, dataframe) -> pd.Series:
    #count only transcripts of same gene!
    tv_count_with_same_uORF = dataframe[(dataframe.gene_symbol == row["gene_symbol"])]["uORF_ID"].values.tolist().count(row["uORF_ID"])
    
    row["tv_count_with_same_uORF"] = tv_count_with_same_uORF
    row["TV_with_uORF/TV_total"] = round((tv_count_with_same_uORF / int(row["tv_count_for_gene"])), 4)
    
    return row
