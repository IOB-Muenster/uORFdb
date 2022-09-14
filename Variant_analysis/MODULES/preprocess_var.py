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

def cut_vcfs(vcf, input_directory, output_directory, metadata_dict, sample_dict):
    
    vcf_dict = {}
        
    VCF = input_directory+vcf

    VCFcut = output_directory+vcf+".cut"
    vcf_text = ""
    with open (VCF, "r") as file:
        for line in file:
            if "##" not in line:
                vcf_text += line
     
    with open (VCFcut, "w") as file:
        file.write(vcf_text)

    variants = pd.read_csv(VCFcut, delimiter = "\t")
    
    variant_count = len(variants)
    
    patient = vcf[0:12]
    file_ids = vcf[13:30]
    N_identifier = vcf[0:12]+"-"+vcf[22:30]
    cancer_type = metadata_dict[N_identifier]["project_id"]
    T_identifier = vcf[0:12]+"-"+vcf[13:21]
    N_identifier = vcf[0:12]+"-"+vcf[22:30]
    
    N_UUID = metadata_dict[N_identifier]["UUID"]
    T_UUID = metadata_dict[T_identifier]["UUID"]
    
    normal_sample = sample_dict[N_UUID]
    tumor_sample = sample_dict[T_UUID]
    
    selection_status = ""
    
    if "FORMAT" in [variants.columns[-1],variants.columns[-2]]:
        selection_status = "error_in_VCF"
        print(selection_status)
    
    elif sorted([normal_sample,tumor_sample]) != sorted([variants.columns[-1],variants.columns[-2]]):
        selection_status = "sample_not_in_VCF"
        print(selection_status)
    
    variant_index = list(range(0,variant_count))
    
    vcf_dict[T_identifier] = {  "cancer_type":cancer_type,
                                "patient":patient,
                                "file_ids":file_ids,
                                "tumor_sample":tumor_sample,
                                "normal_sample":normal_sample,
                                "T_sample_type":metadata_dict[T_identifier]["sample_type"],
                                "N_sample_type":metadata_dict[N_identifier]["sample_type"],
                                "T_UUID":T_UUID,
                                "N_UUID":N_UUID,
                                "variant_count":variant_count,
                                "selection_status":selection_status,
                                "file":VCFcut,
                                "tumor_bam_size":metadata_dict[T_identifier]["file_size"],
                                "gender": metadata_dict[N_identifier]["gender"], # added on 17.08.22
                                "variant_index":variant_index
                                }
    
    vcf_df = pd.DataFrame.from_dict(vcf_dict, "index")
    
    return vcf_df


def select_samples(sample1: pd.Series, dataframe):
    
    #sample_types = ['Blood Derived Normal', 'Metastatic', 'Primary Blood Derived Cancer - Peripheral Blood', 'Primary Tumor', 'Solid Tissue Normal']
    patient = sample1["patient"]
    T_UUID = sample1["T_UUID"]

    matching_patients = dataframe[(dataframe["patient"] == patient) & (dataframe["T_UUID"] != T_UUID)]
    
    if sample1["variant_count"] == 0:
        sample1["selection_status"] += "discarded_no_variants,"

    for index, sample2 in matching_patients.iterrows():
        
        if "Metastatic" in sample1["T_sample_type"] and "Metastatic" not in sample2["T_sample_type"] and sample2["variant_count"] != 0:
            sample1["selection_status"] += "discarded_metastatic,"
        
        if sample2["variant_count"] > sample1["variant_count"]:

            if "Metastatic" in sample1["T_sample_type"] and "Metastatic" in sample2["T_sample_type"]:
                sample1["selection_status"] += "discarded_lower_count,"

            elif "Metastatic" not in sample1["T_sample_type"] and "Metastatic" not in sample2["T_sample_type"]:
                sample1["selection_status"] += "discarded_lower_count,"
                
        elif sample2["variant_count"] == sample1["variant_count"]:
            if sample2["tumor_bam_size"] > sample1["tumor_bam_size"] and "Metastatic" not in sample2["T_sample_type"]:
                sample1["selection_status"] += "discarded_bam_size,"
    
    return sample1


def get_dbSNP_AF(INFO,studies):
    freq_dict = {}
    for study in studies:
        freq_dict[study] = "."
    for entry in INFO.split(";"):
        if "FREQ" in entry:
            for freq in (entry.replace("FREQ=", "")).split("|"):
                study = freq.split(":")[0]
                AF = freq.split(":")[1]
                if study in studies:
                    freq_dict[study] = AF
    return freq_dict


def add_variant_data(row:pd.Series, variants_dataframe):
    
    file = row["file"]
    normal_sample = row["normal_sample"]
    tumor_sample = row["tumor_sample"]
    patient = row["patient"]
    file_ids = row["file_ids"]
    
    matched = variants_dataframe[(variants_dataframe.index == row["variant_index"]) & (variants_dataframe.file == file)]
    
    for index, line in matched.iterrows():

        row["CHROM"] = line["#CHROM"]
        row["POS"] = line["POS"]
        row["REF"] = line["REF"]
        row["ALT"] = line["ALT"]
        
        IDs = line.ID
        FORMAT = line.FORMAT.split(":")
        NORMAL_DATA = line[normal_sample].split(":")
        TUMOR_DATA = line[tumor_sample].split(":")

        REF = line.REF
        POS = line.POS

        row["FROM_0BHO"] = POS-1
        row["TO_0BHO"] = row["FROM_0BHO"]+len(REF)
        
        IDs = IDs.split(";")
         
        row["ClinVar_ID_POS"] = ""
        row["ClinVar_ID_EXACT"] = ""
        row["dbSNP_ID_POS"] = ""
        row["dbSNP_ID_EXACT"] = ""
        
        Allele_frequencies_dict = get_dbSNP_AF(INFO = line.INFO, studies = ["GnomAD","TOPMED","ExAC"])
        for study in Allele_frequencies_dict:
            row[study+"_fromVCF"] = Allele_frequencies_dict[study]
        
        for ID in IDs:  
            if "rs" not in ID and "POS" in ID and ID != ".":
                row["ClinVar_ID_POS"] += ID[4:]+","
            elif "rs" in ID and "POS" in ID and ID != ".":
                row["dbSNP_ID_POS"] += ID[4:]+","
            
            elif "rs" not in ID and "POS" not in ID and ID != ".":
                if row["ClinVar_ID_EXACT"] != "":
                    print("Multiple ClinVar_ID_EXACT???")# (I added ";" before)
                row["ClinVar_ID_EXACT"] += ID
            elif "rs" in ID and "POS" not in ID and ID != ".":
                if row["dbSNP_ID_EXACT"] != "":
                    print("Multiple dbSNP_ID_EXACT???")
                row["dbSNP_ID_EXACT"] += ID
                
        if row["ClinVar_ID_POS"] == "":
            row["ClinVar_ID_POS"] = "."
        if row["ClinVar_ID_EXACT"] == "":
            row["ClinVar_ID_EXACT"] = "."
        if row["dbSNP_ID_POS"] == "":
            row["dbSNP_ID_POS"] = "."
        if row["dbSNP_ID_EXACT"] == "":
            row["dbSNP_ID_EXACT"] = "."
        
        tumor_dict = {}
        normal_dict = {}
        
        for i in range(0,len(FORMAT)):
            normal_dict[FORMAT[i]] = NORMAL_DATA[i]
            tumor_dict[FORMAT[i]] = TUMOR_DATA[i]
        
        row["N_GT"] = normal_dict["GT"]
        row["T_GT"] = tumor_dict["GT"]
        row["N_SomaticMutationFrequency"] = normal_dict["AF"]
        row["T_SomaticMutationFrequency"] = tumor_dict["AF"] 
        row["Proportion_T_SMF/N_SMF"] = float(tumor_dict["AF"]) / float(normal_dict["AF"])
        row["CUSTOM_FILTER"] = ""
        row["N_REF_AlleleDepth"] = normal_dict["AD"].split(",")[0]
        row["N_ALT_AlleleDepth"] = normal_dict["AD"].split(",")[1]
        row["N_ReadDepth"] = normal_dict["DP"]
        row["T_REF_AlleleDepth"] = tumor_dict["AD"].split(",")[0]
        row["T_ALT_AlleleDepth"] = tumor_dict["AD"].split(",")[1]
        row["T_ReadDepth"] = tumor_dict["DP"]
        row["VCF_normal_sample"] = normal_sample
        row["VCF_tumor_sample"] = tumor_sample
        row["patient"] = patient
        row["file_ids"] = file_ids
        
        if row["Proportion_T_SMF/N_SMF"] < 4:
            row["CUSTOM_FILTER"] += "NTratio < 4,"
        if int(row["N_ReadDepth"]) < 10 or int(row["T_ReadDepth"]) < 10:
            row["CUSTOM_FILTER"] += "ReadDepth < 10,"
        if int(row["T_ALT_AlleleDepth"]) < 3:
            row["CUSTOM_FILTER"] += "AlleleDepth < 3,"
            
        if row["CUSTOM_FILTER"] == "":
            row["CUSTOM_FILTER"] = "PASS"
        else:
            row["CUSTOM_FILTER"] = "FAIL:"+row["CUSTOM_FILTER"]
            
    return row

