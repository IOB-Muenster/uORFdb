#!/usr/bin/env perl


# AUTHOR

# Felix Manske, felix.manske@uni-muenster.de


# COPYRIGHT

# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright
#  notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
#  notice, this list of conditions and the following disclaimer in the
#  documentation and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR AS IS AND ANY EXPRESS OR IMPLIED WARRANTIES,
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO,  PROCUREMENT  OF  SUBSTITUTE GOODS  OR  SERVICES;
# LOSS  OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
# OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.


#=================================================================================================================#
# Get gene metadata for uORFs, transcripts and exons
#-----------------------------------------------------------------------------------------------------------------#
#
# DESCRIPTION
#
# Uses the NCBI metadata files all.gene_info and gene2refseq to provide gene metadata
# for a file with transcript IDs. Matching ignores the version tag of IDs.
# Transcripts from pseudo genes (""type of gene"" column = pseudo) and suppressed transcripts
# will be automatically removed.
#
#
# USAGE
#
# getGeneData.pl -input INPUTF  -inputtype uorf | exon | transcript -gene GENECONFF | -prepare
#		[-assembly ASSEMBLY -taxid TAXID] [-verbose]
#
# The script annotates a file of uORFs, transcripts or exons (INPUTF) with gene
# metadata from a gene config file GENECONFF. You must specify the input file type as uorf, exon
# or transcript by specifying -inputtype. See EXAMPLES for the file format of each type.
# If GENECONFF does not exist, setting -prepare will create the file and continue with the
# downstream analysis.
# By specifying the optional parameter -taxid the analysis can be limited to a certain taxonomy ID
# (TAXID). Optionally setting -assembly will replace the assembly info from NCBI with the custom
# string ASSEMBLY. It is not recommended to replace the assembly info without specifying a taxonomy
# ID, unless all taxa in the INPUTF have the same assembly name. The choices for TAXID and ASSEMBLY
# will only affect the analysis. Preparation of the gene config file (-prepare) is not affected.
# The script will also create a blacklist file of suppressed transcripts or transcripts which belong
# to pseudo genes in the directory of the INPUTF. This file must be used for the downstream
# preparation of transcripts for the database.
# Use the flag -verbose to turn on detailed messages.
#
#
# EXAMPLES
#
# The script expects an INPUTF in uORF format to have the transcript ID and gene symbol in the
# 4th and 6th column, respectively (counting from 0). The transcript file must have the transcript
# ID in the very first column. No gene symbol is expected. For the exon file, the transcript ID
# must be in the second column (counting from 0). Again, no gene symbol is expected. All files must
# be tab-delimited.
#
# Prepare a full gene config file and conduct analysis only for human (9606) genes.
# Use the custom assembly string "hg38" in the annotated uORF output:
#
#	./getGeneData.pl -input INPUTF -inputtype uorf -prepare -assembly hg38 -taxid 9606
#
# Use a pre-computed gene config file and annotate all uORFs with genes. Use the
# assembly info from NCBI:
#
#	./getGeneData.pl -input INPUTF -inputtype uorf -gene GENECONFF
#
#
# OUTPUT
#
# The output will be written to the directory of the INPUTF as INPUTF.annotated. For each line in
# the INPUTF, the output will show the raw line from the INPUTF and the metadata. If no metadata was
# found, the value ERROR will be appended. If -prepare was specified, a gene config file will be created
# in the directory of the INPUTF (gene.config). This file is not limited to the taxid specified as -taxid
# and ignores the optional custom assembly string given via -assembly.
# The script will write transcript IDs that are suppressed and/or stem from genes of type pseudo to
# INPUTF.transcr.suppressed. It will write transcripts from INPUTF that could not be matched for
# other reasons to INPUTF.transcr.unmatched
# Temporary files will be created in the the directory of the INPUTF.
#
#
# CAVEATS
#
# The script uses a lot of RAM, since it attempts to load the NCBI metadata files into a hash. 60GB of RAM
# should be safe. This script will use the assembly info provided by NCBI which is the latest assembly given
# for a particular transcript. However, the same transcript may appear multiple times with different assemblies.
# This script will use the last seen assembly. For the analysis it is helpful to set -assembly to the
# assembly used for uORF calling.
#
#
# DEPENDENCIES
#
# *) -prepare only
# curl
# grep
# gunzip
# join
# sort
#=================================================================================================================#


use strict;
use warnings;
use Getopt::Long;
use File::Basename;


#------------------------------------------------------------------------------------------#
# Read in CLI arguments
#------------------------------------------------------------------------------------------#
my $inputF = "";
my $userAssembly = "";
my $userTaxId = "";
my $geneF = "";
my $inputType = "";

my $isPrepare = 0;
my $isVerbose = 0;

my $dlLink1 = "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/All_Data.gene_info.gz";
my $dlLink2 = "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz";

GetOptions ("input=s" => \$inputF,
			"assembly=s" => \$userAssembly,
			"taxid=i" => \$userTaxId,
			"gene=s" => \$geneF,
			"prepare" => \$isPrepare,
			"verbose" => \$isVerbose,
			"inputtype=s" => \$inputType,
			) or die("ERROR: Unexpected argument\n");

die "ERROR: Missing essential argument -input\n" if (not $inputF);
die "ERROR: Missing essential argument -gene\n" if (not $geneF and $isPrepare == 0);
die "ERROR: Invalid or empty gene config file ->$geneF<-\n" if ($geneF and not -s $geneF and $isPrepare == 0);

die "ERROR: Input file ->", $inputF, "<- does not exist" if (! -s $inputF);
die "ERROR: Please specify an input file type" if (not $inputType);

$inputType = lc($inputType);
if ($inputType ne "transcript" and $inputType ne "exon" and $inputType ne "uorf") {
	die "ERROR: Invalid input file type: ->$inputType<-. Must be ->exon<-, ->transcript<- or ->uorf<-\n";
}

my $outDir = dirname($inputF);
my $outF = $inputF . ".annotated";
my $supprF = $inputF . ".transcr.suppressed";
my $unmatchedF = $inputF . ".transcr.unmatched";


# Internal hashes
my %transcr2data = ();
my %suppressedIds = ();
my %unmatchedIds = ();


#------------------------------------------------------------------------------------------#
# Prepare the geneF is necessary
#------------------------------------------------------------------------------------------#
if ($isPrepare == 1) {
	print "INFO: Preparing gene config file\n";
	print "INFO: Processing gene_info\n";
		
	# Download the files, remove the header and sort by gene ID. Status code 0 is error for system
	system("curl -s -o - $dlLink1 | gunzip -c | grep -Ev \"^#\" | sort -k 2,2 > $outDir/all.gene_info") and die "ERROR: Could not download or prepare All_Data.gene_info.gz";
	print "INFO: Processing gene2refseq\n";
	system("curl -s -o - $dlLink2 | gunzip -c | grep -Ev \"^#\" | sort -k 2,2 > $outDir/gene2refseq") and die "ERROR: Could not download or prepare gene2refseq.gz";
	
	# Join the files by the gene ID
	print "INFO: Merging by gene ID\n";
	system("join -j 2 -t \"\t\" $outDir/all.gene_info $outDir/gene2refseq > $outDir/tmp.gene.config");
	
	# Only keep needed columns, reformat, deduplicate
	print "INFO: Reformatting\n";
	open (CONF, "<", "$outDir/tmp.gene.config") or die "ERROR: Could not open gene config file, $!";
	while(<CONF>) {
		chomp($_);
		
		print "DEBUG: Read ", $., " lines\n" if ($. % 500000 == 0 and $isVerbose == 1);
		
		my @splits = split("\t", $_, -1);
		my $ncbiID = $splits[18];
		$ncbiID = "" if ($ncbiID eq "-");
		
		next if (not $ncbiID);
		
		# Loose version tag for improved matching
		$ncbiID =~ s/\.\d+$//;
		
		my $taxId = $splits[1];
		my $assembly = $splits[27];
		my $chr = $splits[6];
		my $symb = $splits[2];
		my $type = $splits[9];
		my $offSymb = $splits[10];
		my $aliases = $splits[4];
		my $offName = $splits[11];
		my $otherNames = $splits[13];
		my $status = $splits[17];
		my $descr = $splits[8];
		my $geneId = $splits[0];
		
		# Set to "", if "-"
		$taxId = "" if ($taxId eq "-"); 
		$assembly = "" if ($assembly eq "-");
		$chr = "" if ($chr eq "-"); 
		$symb = "" if ($symb eq "-");
		$type = "" if ($type eq "-");
		$offSymb = "" if ($offSymb eq "-");
		$aliases = "" if ($aliases eq "-");
		$offName = "" if ($offName eq "-");
		$otherNames = "" if ($otherNames eq "-");
		$descr = "" if ($descr eq "-");
		$geneId = "" if ($geneId eq "-");
		
		next if (not $geneId);
		
		$chr = "chr".$chr if($chr and $taxId and $chr !~ m/^chr/ and $taxId == 9606);
		$offSymb = $symb if (not $offSymb);		
		
		if ($offSymb) {
			$offSymb =~ s/^ //;
			$offSymb =~ s/ $//;
		}
		else {
			die "ERROR: Missing symbol for gene ID $geneId";
		}				
				
		# I want to get all names and the description is fine, if
		# it does not already appear in the official names
		my $descrEscaped = quotemeta($descr);
		my $offNescaped = quotemeta($offName);

		if ($otherNames) {		
			if ($otherNames =~ m/^([^|]+\|)*$offNescaped(\|[^|]+)*$/i or not $offName) {
				$offName = $otherNames
			}
			else {
				$offName .= "|$otherNames";
			}
		}

		if ($descr) {
			if (not $offName) {
				$offName = $descr;
			}
			elsif ($offName !~ m/^([^|]+\|)*$descrEscaped(\|[^|]+)*$/i) {
				$offName .= "|$descr";
			}
		}		

		$offName =~ s/\|/#/g;
		$offName =~ s/^#+//;
		$offName =~ s/^ //;
		$offName =~ s/ $//;
				
		$aliases =~ s/\|/#/g;
		$aliases =~ s/^ //;
		$aliases =~ s/ $//;
				
		# Remove the symbol from the aliases, if necessary
		my $offSescaped = quotemeta($offSymb);
		$aliases =~ s/^$offSescaped#//;
		$aliases =~ s/#$offSescaped$//;
		$aliases =~ s/#$offSescaped#/#/;
		$aliases =~ s/^$offSescaped$//;
		
		if (exists $transcr2data{$ncbiID}){
			# Update empty values, warn on conflicting non-empty values and delete those.
			# I cannot change the assembly, since that changes frequently even for the
			# same transcripts. Assembly is now the last seen assembly.
			my @keys = ("taxid", "symbol", "aliases", "names", "geneid", "chr", "status", "type");
			my @values = ($taxId, $offSymb, $aliases, $offName, $geneId, $chr, $status, $type);
			
			# Only ambiguous transcript IDs are undef
			next if (not defined $transcr2data{$ncbiID});
			
			for (my $i = 0; $i <= $#keys; $i++) {
				my $key = $keys[$i];
				my $value = $values[$i];
				
				# The value in the hash is not the current value
				if ($transcr2data{$ncbiID}->{$key} ne $value) {
					# Skip, if current value is empty
					if (not $value) {
						next;
					}
					else {
						# Update value in hash, if that is empty
						if (not $transcr2data{$ncbiID}->{$key}) {
							$transcr2data{$ncbiID}->{$key} = $value
						}
						# Warn on conflicting non-empty value and set it to undef
						else {
							print "WARNING: Ambiguous $key for transcript $ncbiID " .
								"Current: $value; last seen: $transcr2data{$ncbiID}->{$key}\n";
							# Mark the transcript as ambiguous.	
							$transcr2data{$ncbiID} = undef;
							print "WARNING: Deleted entry\n";
							last;
						}
					}
				}
			}
		}
		else {
				$transcr2data{$ncbiID} = {"taxid" => $taxId, "assembly" => $assembly, "symbol" => $offSymb, "aliases" => $aliases,
					"names" => $offName, "geneid" => $geneId, "chr" => $chr, "status" => $status, "type" => $type};
		}
	}
	close(CONF);
	
	# Write the output
	print "INFO: Writing gene config file\n";
	
	open(FINALCONF, ">", "$outDir/gene.config") or die "ERROR: Could not open final gene config file, $!";
	
	print FINALCONF "#NCBI_ID\tTax_ID\tAssembly\tChr\tSymbol\tAliases\tNames\tGene_ID\tTranscrStatus\tGeneType\n";
	
	foreach my $ncbiID (keys(%transcr2data)) {
		# Remove ambiguous transcript IDs
		if (not defined $transcr2data{$ncbiID}) {
			delete $transcr2data{$ncbiID};
			next;
		}
		
		# Print the full gene config file with the data from NCBI. Useful for follow-up analyses.
		print FINALCONF $ncbiID, "\t", $transcr2data{$ncbiID}->{"taxid"}, "\t", $transcr2data{$ncbiID}->{"assembly"}, "\t",
			$transcr2data{$ncbiID}->{"chr"}, "\t", $transcr2data{$ncbiID}->{"symbol"}, "\t", $transcr2data{$ncbiID}->{"aliases"}, "\t",
			$transcr2data{$ncbiID}->{"names"}, "\t", $transcr2data{$ncbiID}->{"geneid"}, "\t", $transcr2data{$ncbiID}->{"status"}, "\t",
			$transcr2data{$ncbiID}->{"type"}, "\n";
			
		# If user decided to analyze a particular taxon and/or to replace the assembly string,
		# modify the entries in the hash and delete not needed entries.
		if ($userTaxId and $userTaxId != $transcr2data{$ncbiID}->{"taxid"}) {
			delete $transcr2data{$ncbiID};
			next;
		}
		# Add suppressed ncbiIDs to special hash
		if ($transcr2data{$ncbiID}->{"status"} =~ m/^SUPPRESSED$/i) {
			delete $transcr2data{$ncbiID};
			$suppressedIds{$ncbiID} = "";
			next;
		}
		# Add ncbiIDs from pseudogenes to special hash
		if ($transcr2data{$ncbiID}->{"type"} =~ m/^pseudo$/i) {
			delete $transcr2data{$ncbiID};
			$suppressedIds{$ncbiID} = "";
			next;
		}
		if ($userAssembly) {
			$transcr2data{$ncbiID}->{"assembly"} = $userAssembly;
		}
		
	}
	close(FINALCONF);
}
# The gene config file exists
else {
	print "INFO: Reading gene config file\n";
	
	open (CONF, "<", $geneF) or die "ERROR: Could not open gene config file, $!";
	while(<CONF>) {
		# Skip header
		next if ($_ =~ m/^#/);
		
		print "DEBUG: Read ", $., " lines\n" if ($. % 500000 == 0 and $isVerbose == 1);
		chomp($_);
		
		my ($ncbiID, $taxId, $assembly, $chr, $symbol, $aliases, $names, $geneId, $status, $type) = split("\t", $_, -1);
		
		# Filter by taxid (opt.) and replace assembly string (opt.)
		next if ($userTaxId and $taxId != $userTaxId);
		
		# Add suppressed ncbiIDs to special hash
		if ($status =~ m/^SUPPRESSED$/i) {
			$suppressedIds{$ncbiID} = "";
			next;
		}
		# Add ncbiIDs from pseudo genes to special hash
		if ($type =~ m/^pseudo$/i) {
			$suppressedIds{$ncbiID} = "";
			next;
		}
		
		$assembly = $userAssembly if ($userAssembly);
		
		
		if (exists $transcr2data{$ncbiID}) {
			die "ERROR: Ambiguous transcript ->$ncbiID<-\n";
		}
		else {
			$transcr2data{$ncbiID} = {"taxid" => $taxId, "assembly" => $assembly, "symbol" => $symbol, "aliases" => $aliases,
				"names" => $names, "geneid" => $geneId, "chr" => $chr, "status" => $status, "type" => $type};
		}
	}
	close(CONF);
}


#------------------------------------------------------------------------------------------#
# Read uORF file and add gene metadata.
#------------------------------------------------------------------------------------------#
print "INFO: Adding gene data\n";

open(OUT, ">", $outF) or die "ERROR: Could not open output file, $!";
open(IN, "<", $inputF) or die "ERROR: Could not open input file, $!";
while (<IN>) {
	# Skip header
	next if ($. < 2);
	
	chomp($_);
	
	my @splits = split("\t", $_, -1);
	my ($id, $symbol) = (undef, undef);
	
	# Get transcript ID and symbol depending on the file type
	if ($inputType eq "uorf") {
		$id = $splits[4];
		$symbol = $splits[6];
	}
	elsif ($inputType eq "transcript") {
		$id = $splits[0];
		$symbol = "NA";
	}
	elsif ($inputType eq "exon") {
		$id = $splits[2];
		$symbol = "NA";
	}
	
	if (not defined $id or not defined $symbol) {
		die "ERROR: Could not extract transcript ID or gene symbol from uORF file.\n";
	}
	
	# Delete version tag for easier lookup.
	$id =~ s/\.\d+$//;
	
	if (exists $transcr2data{$id}) {
		print OUT $_, "\t", $transcr2data{$id}->{"taxid"}, "\t", $transcr2data{$id}->{"assembly"}, "\t",
			$transcr2data{$id}->{"chr"}, "\t", $transcr2data{$id}->{"symbol"}, "\t",$transcr2data{$id}->{"aliases"}, "\t",
			$transcr2data{$id}->{"names"}, "\t", $transcr2data{$id}->{"geneid"}, "\n";
	}
	else {
		if(exists $suppressedIds{$id}) {
			print "WARNING: Transcript ->$id<- is suppressed by NCBI or gene is a pseudo gene. Removed all occurences.\n";
			next;
		}
		else {
			print "WARNING: Cannot find metadata for symbol ->$symbol<- and transcript ->$id<-\n";
			$unmatchedIds{$id} = "";
			print OUT $_, "\t", "ERROR", "\n";
		}
	}
}
close(IN);
close(OUT);


# Print all suppressed IDs for the selected taxon. Important for downstream processing
# of transcripts.
open(SUPPR, ">", $supprF) or die "ERROR: Could not open file for suppressed IDs, $!";
print SUPPR join ("\n", keys(%suppressedIds));
close(SUPPR);


# Print all unmatched IDs for the selected taxon. Important for downstream processing
# of transcripts.
open(UNM, ">", $unmatchedF) or die "ERROR: Could not open file for unmatched IDs, $!";
print UNM join ("\n", keys(%unmatchedIds));
close(UNM);

print "\nDONE\n"
