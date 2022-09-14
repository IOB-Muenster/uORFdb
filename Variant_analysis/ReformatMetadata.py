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
# Reformat downloaded JSON metadata to TSV with relevant information
#
# ==========================================================================================

import json
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--directory",help="path to directory with JSON for each cancer type")

args = parser.parse_args()

directory = args.directory

types = ["BRCA","COAD","LAML","LUAD","PRAD","SKCM"]

summary_dict = {}

summary_file = directory+"/file_summary.csv"

all_types = directory+"/TCGA_metadata.csv"
all_types_list = []

for cancer_type in types:
    
    print(cancer_type)

    json_file = directory+"/"+cancer_type+"/"+cancer_type+"_metadata.json"
    
    csv_file = directory+"/"+cancer_type+"/"+cancer_type+"_metadata_incl_gender.csv"
    
    with open (json_file, "r") as file:
        object = json.load(file)
    
    dict = {}
    
    for case in object:
        
        id = case["file_id"]
        
        tissue = (case["cases"][0]["samples"][0]["sample_type"])
        
        if "Normal" in tissue:
            tissue_key = "N"
        elif "Cancer" in tissue or "Tumor" in tissue or "Metastatic" in tissue:
            tissue_key = "T"
        else:
            tissue_key = "unknown"
        
        case_dict = {"UUID" : case["file_id"],
                     "tissue_key" : tissue_key,
                     "patient-tissue" : case["cases"][0]["submitter_id"]+"-"+tissue_key,
                     "patient_id" : case["cases"][0]["submitter_id"],
                     "gender" : case["cases"][0]["demographic"]["gender"],
                     "submission_date" : case["metadata_files"][0]["created_datetime"],
                     "sample_type" : case["cases"][0]["samples"][0]["sample_type"],
                     "project_id" : case["cases"][0]["project"]["project_id"],
                     "submitter_id" : case["cases"][0]["submitter_id"],
                     "sample_id" : case["cases"][0]["samples"][0]["sample_id"],
                     "case_id" : case["associated_entities"][0]["case_id"],
                     "entity_id" : case["associated_entities"][0]["entity_id"],
                     "file_id_UUID" : case["file_id"],
                     "file_name" : case["file_name"],
                     "submitter_id_2" : case["submitter_id"],
                     "file_size" : case["file_size"],
                     "sample+patient" : case["cases"][0]["samples"][0]["sample_type"]+"-"+case["cases"][0]["submitter_id"],
                     "aliquot_date" : case["cases"][0]["samples"][0]["portions"][0]["analytes"][0]["aliquots"][0]["updated_datetime"],
                     "date1" : case["cases"][0]["samples"][0]["portions"][0]["analytes"][0]["updated_datetime"],
                     "date2" : case["cases"][0]["samples"][0]["portions"][0]["updated_datetime"],
                     "metadata_date" : case["metadata_files"][0]["created_datetime"]
                     }
        
        dict[id] = case_dict
        
    data = (pd.DataFrame.from_dict(dict, 'index'))
    
    data2=data.sort_values(by=["patient_id","tissue_key","submission_date"])
    
    #mark duplicates
    duplicate = data2.duplicated("patient-tissue")
    data2.loc[duplicate == True, "tissue_key"] += "-duplicate"
    data2.loc[duplicate == True, "patient-tissue"] += "-duplicate"
    data2.loc[data2["tissue_key"] == "T-duplicate", "tissue_key"] = "T"
    
    no_duplicate = data2.duplicated("patient_id", keep = False)
    data2.loc[no_duplicate != True, "tissue_key"] += "-single_file"
    data2.loc[no_duplicate != True, "patient-tissue"] += "-single_file"
    
    data2.to_csv(csv_file, index = False)
    
    all_types_list.append(data2)
    
    #overview of all cancertypes
    summary_dict[cancer_type]={"cancer_type" : cancer_type,
                               "file_number" : len(data2),
                               "patients" : len(data2["patient_id"].unique()),
                               "N" : len(data2.loc[data2["tissue_key"] == "N"]),
                               "T" : len(data2.loc[data2["tissue_key"] == "T"]),
                               "N-duplicates" : len(data2.loc[data2["tissue_key"] == "N-duplicate"]),
                               "T-duplicates" : len(data2.loc[data2["tissue_key"] == "T-duplicate"]),
                               "T-single_file" : len(data2.loc[data2["tissue_key"] == "T-single_file"]),
                               "N-single_file" : len(data2.loc[data2["tissue_key"] == "N-single_file"]),
                               "files_per_patient" : len(data2)/(len(data2["patient_id"].unique())),
                               }
    
summary = (pd.DataFrame.from_dict(summary_dict, 'index'))
summary.to_csv(summary_file, index = False)

all_types_list = pd.concat(all_types_list)

all_types_list.to_csv(all_types, index = False)
    
for i in summary_dict:
    print(summary_dict[i])
