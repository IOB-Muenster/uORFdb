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


#=================================================================================================#
# Prepare transcript and exon files
#-------------------------------------------------------------------------------------------------#
#
# DESCRIPTION
#
# Prepare the transcript and exon files for uORFdb. Skip transcripts that were suppressed by NCBI.
#
#
# USAGE
#
# makeTranscripts.pl -transcr TRANSCRF -gene GENEF [-taxid TAXID] [-assembly ASSEMBLY]
# 		[-blacklist BLACKLISTF] [-verbose]
#
# Create a transcript and an exon file for uORFdb from TRANSCRF. This is the transcripts
# TSV file from the uORF finder. Optionally, transcripts that appear in BLACKLISTF can be excluded.
# Each record in the output files will be annotated with gene metadata from GENEF.
# By providing an ASSEMBLY, a custom assembly string is used as the assembly name, instead of
# NCBI's name. It is not recommended to replace the assembly without specifying a TAXID, unless
# all taxa in the TRANSCRF have the same assembly name.
#
#
# EXAMPLES
#
# Create and annotate transcript and exon files for human data. Use the custom assembly string hg38:
# 	
#	./makeTranscripts.pl -transcr transcripts.tsv -gene gene.config -taxid 9606
#		-assembly hg38 -blacklist blacklist.transcr
#
#
# OUTPUT
#
# Temporary transcript and exon files for the annotation with gene metadata (out.transcripts and
# out.exons) will be generated in the directory of the TRANSCRF. The annotated files for uORFdb
# out.transcript.annotated and out.exon.annotated will be created in the same directory.
# For both data types, *.suppressed and *.unmatched files will be created in the directory of the
# TRANSCRF. See getGeneData.pl for more details.
#
#
# CAVEATS
#
# To prepare a GENEF, run getGeneData.pl separately. For other caveats, see the section in
# getGeneData.pl.
#
#
# DEPENDENCIES
#
# *) getGeneData.pl (must be in the directory of this script)
#=================================================================================================#


use strict;
use warnings;
use Getopt::Long;
use File::Basename;


#-------------------------------------------------------------------------------------------------#
# Read args from CLI
#-------------------------------------------------------------------------------------------------#
my $transcrF = "";
my $blacklistF = "";
my $geneF = "";
my $userAssembly = "";
my $userTaxID = "";
my $isVerbose = 0;

GetOptions ("transcr=s" => \$transcrF,
			"blacklist=s" => \$blacklistF,
			"assembly=s" => \$userAssembly,
			"taxid=i" => \$userTaxID,
			"gene=s" => \$geneF,
			"verbose" => \$isVerbose,
			) or die("ERROR: Unexpected argument\n");

die "ERROR: Missing essential argument -transcr\n" if (not $transcrF);
die "ERROR: Missing essential argument -gene\n" if (not $geneF);

# Output files for transcript and exon
my $outDir = dirname($transcrF);
my $outTF = $outDir . "/" . "out.transcripts";
my $outEF = $outDir . "/" . "out.exons";

my $cwd = dirname($0);


#-------------------------------------------------------------------------------------------------#
# Read blacklist of NCBI IDs, if applicable
#-------------------------------------------------------------------------------------------------#
print "INFO: Reading blacklist\n";

my %blacklisted = ();
if (-s $blacklistF) {
	open(BLACKLIST, "<", $blacklistF) or die "ERROR: Could not open blacklist file, $!"; 
	while(<BLACKLIST>) {
		chomp($_);
		$blacklisted{$_} = "";
	}
	close(BLACKLIST);
}


#-------------------------------------------------------------------------------------------------#
# Process transcript file
#-------------------------------------------------------------------------------------------------#
print "INFO: Processing transcript file\n";

open(OUTT, ">", $outTF) or die "ERROR: Could not open output transcript file, $!";
print OUTT "ncbiID\tchrom\ttxStart\ttxLength\ttxLength_noint\tcdsStart\tcdsEnd\ttlsLen_noint\tstrand\tkozakContext\tkozakStrength\n";

open(OUTE, ">", $outEF) or die "ERROR: Could not open output exon file, $!";
print OUTE "exonStart\texonEnd\tncbiID\tstrand\n";

open(TRANSCR, "<", $transcrF) or die "ERROR: Could not open transcript file, $!"; 
while(<TRANSCR>) {
	# Skip header
	next if ($. < 2);
	
	chomp($_);
		
	my @splits = split("\t", $_, -1);
	
	my $ncbiID = $splits[1];
	
	# Apply blacklist of suppressed NCBI IDs.
	# Blacklist has no version tag.
	my $ID = $ncbiID =~ s/\.\d+$//r;
	if (exists $blacklisted{$ID}) {
		print "DEBUG: Skipped blacklisted transcript ->$ncbiID<-\n";
		next;
	}
	
	my $chr = $splits[2];
	my $strand = $splits[3];
	my $txStart = $splits[4];
	my $txEnd =  $splits[5];
	my $cdsStart = $splits[6];
	my $cdsEnd = $splits[7];
	my $exonStart = $splits[9];
	my $exonEnd = $splits[10];
	my $cdsKozak = $splits[26];
	$cdsKozak = "" if ($cdsKozak eq ".");
	my $cdsKozakStrength = $splits[27];
	$cdsKozakStrength = "" if ($cdsKozakStrength eq ".");
	
	my $selectionStatus = $splits[30];
	my $duplStatus = $splits[31];
	
	# Skip invalid chromosomes or duplicate IDs
	next if ($selectionStatus =~ m/invalid_chr/ or $duplStatus ne '.');
	
	# R* columns are (based on) relative coordinates. Introns already removed.
	my $RtxStart = $splits[17];
	my $RtxEnd = $splits[18];
	my $RtlsLen = $splits[25];
	
	my $txLen = $txEnd - $txStart;
	my $txLenNoInt = $RtxEnd - $RtxStart;
	
	
	#---------------------------------------------------------------------------------------------#
	# Write transcript file
	#---------------------------------------------------------------------------------------------#
	print OUTT $ncbiID, "\t", $chr, "\t", $txStart, "\t", $txLen, "\t", $txLenNoInt, "\t", $cdsStart, "\t", $cdsEnd,
		"\t", $RtlsLen, "\t", $strand, "\t", $cdsKozak , "\t", $cdsKozakStrength, "\n";
	
	
	#---------------------------------------------------------------------------------------------#
	# Write exon file
	#---------------------------------------------------------------------------------------------#
	# Remove trailing comma and print single pairs of start and stop
	$exonStart =~ s/,$//;
	$exonEnd =~ s/,$//;
	
	my @exStarts = split(",", $exonStart);
	my @exEnds = split(",", $exonEnd);
	
	die "ERROR: Number of exon starts and stops does not match\n" if (scalar(@exStarts) != scalar(@exEnds));
	
	for (my $i = 0; $i <= $#exStarts; $i++) {
		print OUTE $exStarts[$i], "\t", $exEnds[$i], "\t", $ncbiID, "\t", $strand, "\n";
	}
}
close(TRANSCR);
close(OUTT);
close(OUTE);


#---------------------------------------------------------------------------------------------#
# Annotate transcript and exon file with gene metadata
#---------------------------------------------------------------------------------------------#
print "INFO: Getting gene metadata\n";

# Return code 0 is error for system
my $commandTr = "$cwd/getGeneData.pl -input $outTF -gene $geneF -v -inputtype transcript";
my $commandEx = "$cwd/getGeneData.pl -input $outEF -gene $geneF -v -inputtype exon";

if ($userAssembly) {
	$commandTr .= " -assembly $userAssembly";
	$commandEx .= " -assembly $userAssembly";
}
if ($userTaxID) {
	$commandTr .= " -taxid $userTaxID";
	$commandEx .= " -taxid $userTaxID";
}

system($commandTr) and die "ERROR: Cannot annotate transcript file";
system($commandEx) and die "ERROR: Cannot annotate exon file";

print "\nDONE\n";	