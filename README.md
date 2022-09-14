# Supplemental scripts for uORFdb.
This repository contains supplemental scripts for the publication of uORFdb (https://www.bioinformatics.uni-muenster.de/tools/uorfdb).
Instructions on how to run the scripts can be found in the Supplementary Methods.

## Contents
- ### uORF detection
	- uORF_Finder
	- extract_uORFexon_regions.py

- ### Variant analysis
	- MODULES
	- ReformatMetadata.py
	- bam2fastq.sh
	- bwa.sh
	- clean.sh
	- sort.sh
	- merge.sh
	- markdups.sh
	- baserecalibrator.sh
	- applybqsr.sh
	- fingerprint.sh
	- mutect2.sh
	- getpileup.sh
	- contamination.sh
	- readorientation.sh
	- filtermutectcalls.sh
	- extract_uORFexon_POS.py
	- filter_and_annotate_VCFs.sh
	- VarPreprocessing.py
	- extract_sample_names.sh
	- CoveragePreprocessing.py
	- CoverageMapping.py
	- convert_JSON.py
	- AlleleFreqMapping.py
	- uORFVarAnalysis.py

- ### Preparations for the database insert
	- getGeneData.pl
	- makeTranscripts.pl
	- makePublFile.pl
	- makeTaxonomy.pl
	- makeVariations.pl

- ### The database
	- schema.sql
	- lib
	- importDB.pl	

- ### Database exports
	- makeRSS.pl
	- dumpRelations.sh
	- genomebrowser.pl	

## Cite
If you use these scripts in your research, please cite this repository.
