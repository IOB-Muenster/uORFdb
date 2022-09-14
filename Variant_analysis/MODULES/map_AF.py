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
# Functions used in Variant Analysis workflow
#
# ==========================================================================================

import pandas as pd

def add_dbSNP_numbers(row:pd.Series,dbSNP_dict):

    studies = { "ExAC.1":"ExAC_fromVCF",
                "GnomAD.3":"GnomAD_fromVCF",
                "TOPMED.2":"TOPMED_fromVCF",
                "TOPMED.3":"TOPMED_fromVCF"}
    
    row["ID.index"] = "."
    for study in studies:
        row[study+"_ref/alt/total"] = "."
    row["match_info"] = ""

    identifier = row["dbSNP_ID_EXACT"].split(";")[0]+"_"+row["REF"]+"_"+row["ALT"]
    identifier_NREF = row["dbSNP_ID_EXACT"].split(";")[0]+"_"+"N"+"_"+row["ALT"]
    
    ##DEBUG
    if row["dbSNP_ID_EXACT"] != "." and identifier_NREF in dbSNP_dict:
        print(identifier_NREF)
    
    if row["dbSNP_ID_EXACT"] != "." and identifier in dbSNP_dict:

        IDindex = dbSNP_dict[identifier]["ID.index"]
        row["ID.index"] = IDindex

        for study in studies:
            
            VCF_string = row[studies[study]]
            
            if VCF_string != ".":

                VCF_value = VCF_string.split(",")[int(IDindex.split(".")[1])]
                
                if "/" in str(dbSNP_dict[identifier]["REF_"+study]) and "/" in str(dbSNP_dict[identifier]["ALT_"+study]) and VCF_value != ".":
                    '''
                    COMPARISON OF ALLELE FREQUENCIES FROM JSON AND VCF
                    example data:
                    VCF_string = "1,2.473e-05,.,." (index in ID.index should be correct index of this list)
                    JSON_string = "3/121298" -> AF is calculated from counts
                    JSON calculated AF is compared to value at index from VCF data
                    number of digits after comma is taken over from VCF data (already rounded, JSON data may be more exact)
                    '''
                    ref_count = dbSNP_dict[identifier]["REF_"+study].split("/")[0]
                    alt_count = dbSNP_dict[identifier]["ALT_"+study].split("/")[0]
                    total_count = dbSNP_dict[identifier]["REF_"+study].split("/")[1]
                    row[study+"_ref/alt/total"] = ref_count+","+alt_count+","+total_count
                    
                    VCF_float = format(float(VCF_value),'.10f').rstrip('0').rstrip('.')
                    if float(VCF_float) == 0:
                        VCF_float = str(round(float(VCF_float),0))
                    digits = '.'+str(len(str(VCF_float).split(".")[1]))+'f'  
                        
                    JSON_string = dbSNP_dict[identifier]["ALT_"+study]
                    JSON_float = format(float(int(JSON_string.split("/")[0])/int(JSON_string.split("/")[1])),digits)
                    
                    if VCF_float == JSON_float:
                        row["match_info"] += "MATCH:"+study+":"+str(VCF_float)+","
                    else:
                        row["match_info"] += "MISMATCH:"+study+":"+str(VCF_float)+"!="+str(JSON_float)+","
                        
                else:
                    row[study+"_ref/alt/total"] = "."
            else:
                row[study+"_ref/alt/total"] = "."

    else:
        for study in studies:
            row[study+"_ref/alt/total"] = "."
            
    if row["match_info"] == "":
        row["match_info"] = "."
        
    if "MISMATCH" in row["match_info"]:
        if "MISMATCH:TOPMED.2" in row["match_info"] and "MATCH:TOPMED.3" in row["match_info"] and row["match_info"].count("MISMATCH") == 1:
            row["match_info2"] = "PASS"
            row["TOPMED.2_ref/alt/total"] = "."
        else:
            row["match_info2"] = "MISMATCH"
    else:
        row["match_info2"] = "PASS"
        # in case that AF is same in TOPMED studies -> select newer study
        if "MATCH:TOPMED.2" in row["match_info"] and "MATCH:TOPMED.3" in row["match_info"]:
            row["TOPMED.2_ref/alt/total"] = "."

    return row


def merge_variants(variant_dataframe, identifier_columns, merge_columns, sum_columns, cancer_types):
    
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
                
                elif column in sum_columns:
                    unique_variants[ID][column] += int((variants[index][column]))
                    
        elif ID not in unique_variants:
            unique_variants[ID] = {}
            for column in variants[index]:
                if column in merge_columns:
                    unique_variants[ID][column] = [(variants[index][column])]
                    
                elif column in sum_columns:
                    unique_variants[ID][column] = int((variants[index][column]))
                    
                else:
                    unique_variants[ID][column] = (variants[index][column])

        unique_variants[ID]["VCF_PatCount"] = len(unique_variants[ID]["patient"])
       
    merged_variants = pd.DataFrame.from_dict(unique_variants, 'index')
    
    return merged_variants


def addCancerAF(row: pd.Series):

    groups = ["TCGA","TCGA-BRCA","TCGA-COAD","TCGA-LAML","TCGA-LUAD","TCGA-PRAD","TCGA-SKCM"]

    for group in groups:
        if int(row[group+"_AlleleCount"]) != 0:
            row[group+"_AF"] = round( int(row[group+"_ALT_AlleleCount"]) / int(row[group+"_AlleleCount"]) ,8)
        else:
            row[group+"_AF"] = "."
        
    return row


def convert_list_to_text(row: pd.Series, dataframe, sep):

    for column in dataframe.columns.values:
    
        if isinstance(row[column], list):
            text = ""
            for i in row[column]:
                text += str(i)+sep
            row[column] = text
            
    return row

