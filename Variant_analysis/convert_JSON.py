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
# convert downloaded dbSNP JSON files to manageable tables for each chromosome 
# https://api.ncbi.nlm.nih.gov/variation/v0/
#
# ==========================================================================================


#----------------------------------------------------------
# PACKAGES
#----------------------------------------------------------
import pandas as pd
import json
import time
import re
import argparse
from progress.bar import ChargingBar

#----------------------------------------------------------
# INPUT SECTION
#----------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--chrom",help="chromosome name")
parser.add_argument("--inp",help="input_dir")
parser.add_argument("--out",help="output_dir")
parser.add_argument("--ids",help="path to id file")
args = parser.parse_args()
chrom = args.chrom
output_dir = args.out
input_dir = args.inp
id_file = args.ids

whole_start_time = time.time()

json_file = input_dir+"/refsnp-"+chrom+".json"
output_file = output_dir+"/"+chrom+"_uORFexonv2_AF_nan.tsv"
nonan_output_file = output_dir+"/"+chrom+"_uORFexonv2_AF.tsv"

print("+----------------------------+")
print("| JSON converter version 1.0 |")
print("+----------------------------+")
print("chromosome: "+chrom)
print("selected IDs: "+id_file)
print("input json: "+json_file)
print("output tsv: "+output_file)

#----------------------------------------------------------
# SOME VARIABLES
#----------------------------------------------------------
alias = {"chr1" : "NC_000001.11", "chr2" : "NC_000002.12", "chr3" : "NC_000003.12", "chr4" : "NC_000004.12", "chr5" : "NC_000005.10", "chr6" : "NC_000006.12",
         "chr7" : "NC_000007.14", "chr8" : "NC_000008.11", "chr9" : "NC_000009.12", "chr10" : "NC_000010.11", "chr11" : "NC_000011.10", "chr12" : "NC_000012.12",
         "chr13" : "NC_000013.11", "chr14" : "NC_000014.9", "chr15" : "NC_000015.10", "chr16" : "NC_000016.10", "chr17" : "NC_000017.11", "chr18" : "NC_000018.10",
         "chr19" : "NC_000019.10", "chr20" : "NC_000020.11", "chr21" : "NC_000021.9", "chr22" : "NC_000022.11", "chrX" : "NC_000023.11", "chrY" : "NC_000024.10"}

chromosome = alias[chrom]

print("Loading input data ... ")

ids_df = pd.read_csv(id_file, delimiter = "\t", names = ["chr","id"], header = None)
ids_df = ids_df[(ids_df.chr == chromosome)]

ids = {}
for index, row in ids_df.iterrows():
    ids[row["id"]] = row["chr"]

counter = 0
uORF_id_counter = 0
discarded = 0

print("Processing JSON lines ... ")

with open (json_file, "r") as file:
    json_line_count = 0
    for line in file:
        json_line_count += 1
        
#----------------------------------------------------------
# iterate over JSON lines
#----------------------------------------------------------
with open (json_file, "r") as file:
    
    data = {}
    
    with ChargingBar(max=json_line_count) as bar:
    
        for line in file:
            
            counter += 1
                
            rs_id = re.search('"refsnp_id":("[0-9]*")', line)[0].split(":")[1][1:-1]
                
            if "rs"+rs_id in ids:
                
                if "TOPMED" not in line and "GnomAD" not in line and "ExAC" not in line:
                    discarded += 1
                    
                uORF_id_counter += 1
            
                json_data = json.loads(line)
                
                rsID = str(json_data["refsnp_id"])
               
                #----------------------------------------------------------
                # determine which one is reference and which is alternative
                #----------------------------------------------------------
                placements_with_allele = json_data["primary_snapshot_data"]["placements_with_allele"]
                
                true_counter = 0
                
                for i in range (0,len(placements_with_allele)):
                    
                    if placements_with_allele[i]["is_ptlp"] == True:
                        
                        true_counter += 1
                        
                        reference_index = []
                        alternative_index = []
                    
                        reference_seq = ""
                        alternative_seq = {}
                        
                        for j in range (0,len(placements_with_allele[i]["alleles"])):
                            
                            DEL = placements_with_allele[i]["alleles"][j]["allele"]["spdi"]["deleted_sequence"]
                            INS = placements_with_allele[i]["alleles"][j]["allele"]["spdi"]["inserted_sequence"]
                            
                            if DEL == INS:
                                reference_index.append(j)
                                if reference_seq == "":
                                    reference_seq = DEL
                                else:
                                    reference_seq += ","+DEL
                                    print("Caution: multiple alleles for reference!")
                                
                            elif DEL != INS:
                                alternative_index.append(j)
                                alternative_seq[j]=INS
                
                if true_counter > 1:
                    print("Caution: multiple blocks contained 'ptlp:true'!")            
                
                #----------------------------------------------------------
                # extract data
                #----------------------------------------------------------
                for index in alternative_index:
                    
                    ID = rsID+"."+str(index)
                    
                    data[ID] = {}
                    
                    data[ID]["ID.index"] = "rs"+str(ID)
                    
                    data[ID]["ID"] = "rs"+str(rsID)
                    
                    data[ID]["index"] = str(index)
                    
                    data[ID]["REF"] = reference_seq
                    data[ID]["ALT"] = alternative_seq[index]

                    ALT_data = json_data["primary_snapshot_data"]["allele_annotations"][index]["frequency"]
                    
                    for i in range (0,len(ALT_data)):
        
                        study_name = ALT_data[i]["study_name"]
                        study_version =  ALT_data[i]["study_version"]
              
                        if study_name == "TOPMED" or study_name == "GnomAD" or study_name == "ExAC":
    
                            data[ID]["ALT"+"_"+study_name+"."+str(study_version)] = str(ALT_data[i]["allele_count"])+"/"+str(ALT_data[i]["total_count"])
                
                    for index in reference_index:
                       
                        REF_data = json_data["primary_snapshot_data"]["allele_annotations"][index]["frequency"]
                               
                        for i in range (0,len(REF_data)):
            
                            study_name = REF_data[i]["study_name"]
                            study_version =  REF_data[i]["study_version"]
                  
                            if study_name == "TOPMED" or study_name == "GnomAD" or study_name == "ExAC":
            
                                if "REF"+"_"+study_name+"."+str(study_version) in data[ID]:
                                    data[ID]["REF"+"_"+study_name+"."+str(study_version)] += str(REF_data[i]["allele_count"])+"/"+str(REF_data[i]["total_count"])
                                else:
                                    data[ID]["REF"+"_"+study_name+"."+str(study_version)] = str(REF_data[i]["allele_count"])+"/"+str(REF_data[i]["total_count"])
                                    
            bar.next()

#----------------------------------------------------------
# OUTPUT SECTION
#----------------------------------------------------------
print("Writing to output ... ")
dataframe = pd.DataFrame.from_dict(data, 'index')
dataframe.to_csv(output_file, sep="\t", index = False)

dataframe.fillna(".", inplace=True)
dataframe.to_csv(nonan_output_file, sep="\t", index = False)

print("Finished.")
print("Lines that matches uORF IDs: "+str(uORF_id_counter))
print("uORF-relevant IDs for this chromosome: "+str(len(ids)))
print("Tested lines from refsnp json: "+str(counter))
print("uORF-relevant IDs without relevant study numbers: "+str(discarded))
print("uORF-relevant IDs with relevant study numbers: "+str(uORF_id_counter - discarded))

print("whole time: "+str(round(((time.time() - whole_start_time)/ 60),2))+" min")


            