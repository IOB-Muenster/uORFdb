uORF FINDER v2

USAGE

	uORF_finder.py [-h] [--refseq REFSEQ] [--output OUTPUT] [--fasta FASTA] [-v] [-p] [-d] [-c] [--complement] [--duplicate] [--test21]

	optional arguments:
	  -h, --help           show this help message and exit
	  --refseq REFSEQ      path to refseq file, required
	  --output OUTPUT      path to output directory, required
	  --fasta FASTA        path to fasta file that contains seq of all chromosomes, required
	  -v, --verbose        verbose output, default = False
	  -p, --preprocessing  only run preprocessing, default = False
	  -d, --uorfdetection  only run uORF detection, default = False
	  -c, --chrom_output   enable one additional output file per chromosome, recommended!, default = False
	  --complement         only run chromosomes which were not processed yet (--chrom_output is set to True as well), default = False
	  --duplicate          includes duplicate transcript ids on different chromosomes, default = False
	  --test21             test run chr21 (does not work for all species, recommended for homo sapiens)
 
EXAMPLE COMMAND

	./uORF_finder.py \
	-v \
	-c \
	--refseq <refseq_path> \
	--fasta <fasta_path> \
	--output <output_dir> \
	--duplicate
	  