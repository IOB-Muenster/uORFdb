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


#==========================================================================================#
# Create taxonomy input for uORFdb
#------------------------------------------------------------------------------------------#
#
# DESCRIPTION
#
# Read the taxonomy files from NCBI and create a file which can be inserted
# into uORFdb.
#
#
# USAGE
#
# makeTaxonomy.pl -i [TAXDIR] -o [OUTF]
#
# TAXDIR is the directory which contains the extracted taxdump data
# from NCBI (ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz). The output
# is written to OUTF. Temporary files will be created in the directory of
# the output file.
#
#
# OUTPUT
#
# A tab-separated file containing taxid, name and lineage. The full lineage from the NCBI
# taxonomy database is reported. If there are duplicate strings of taxon name + lineage,
# only one is kept, even if the taxonomy ID is different.
# The output file only contains data about taxa with the type "*scientific*" in names.dmp.
# Temporary files are created in the directory of the OUTF.
# 
#
# DEPENDENCIES
#
# GNU Awk 4.1.4 or later.
# NCBI Taxonomy Toolkit v0.8.0 or later (https://bioinf.shenwei.me/taxonkit/download/).
#==========================================================================================#


use strict;
use warnings;
use File::Basename;
use Data::Dumper;
use Getopt::Long;


my $inDir = "";
my $outF = "";
my $outDir = "";

GetOptions ("i=s" => \$inDir,
			"o=s" => \$outF,
			) or die("ERROR: Unexpected argument\n");

die "ERROR: Missing essential argument -i or -o\n" if (not $inDir or not $outF);
print "ERROR: ->" . $inDir . "<- is not a valid directory\n" if (not -d $inDir);

$outDir = dirname($outF);
print "ERROR: Attempting to write to invalid directory ->" . $outDir . "<-\n" if (not -d $outDir);

print "INFO: Reading taxonomy from ->" . $inDir . "<-\n";
print "INFO: Writing output to ->" . $outF . "<-\n";


#------------------------------------------------------------------------------#
# Extract all scientific names and their taxids
#------------------------------------------------------------------------------#
print "INFO: Extracting names for taxa\n";

my $namesF = "$inDir/names.dmp";
my $tmpOutF = $outDir . "/tmp." . basename($outF);

open(OUT, ">", $tmpOutF) or die "Could not open temporary output file $tmpOutF";
open(NAMES, "<", $namesF) or die "Could not open $namesF from NCBI taxdump.";
while(<NAMES>) {
	
	chomp($_);
	my ($id, $name, $name2, $type) = split(/\|/, $_);
	$id =~ s/\s+//g;
	$name =~ s/^\s+//;
	$name =~ s/\s+$//;
	$name =~ s/\s\s+/ /g;
	$type =~ s/\s+//g;
    
	if ($type =~ m/scientific/) {
		print OUT $id, "\t", $name, "\n"
	}
}
close(NAMES);
close(OUT);


#--------------------------------------------------------------------------------#
# Get the full lineage for each taxid
#--------------------------------------------------------------------------------#
print "INFO: Getting lineage\n";
my $tmpOutF2 = $outDir . "/tmp2." . basename($outF);

my $command = "taxonkit --data-dir $inDir lineage " . $tmpOutF . " > $tmpOutF2";
	
# Status code 0 is error for system()
system($command) and die "ERROR: Could not get lineage with taxonkit";
system("rm $tmpOutF") and die "ERROR: Could not remove temporary output file";


#---------------------------------------------------------------------------------#
# Only keep unique combinations of taxon + lineage, even if taxid is different
#---------------------------------------------------------------------------------#
print "INFO: Making lineage unique\n";

my $tmpOutF3 = $outDir . "/tmp3." . basename($outF);
$command = "awk -F\"\t\" 'BEGIN{OFS = \"\t\"} !seen[\$2, \$3]++' " . $tmpOutF2 ." > $tmpOutF3";

# Status code 0 is error for system()
system($command) and die "ERROR: Could not get unique lineages";
system("rm $tmpOutF2") and die "ERROR: Could not remove temporary output file";


#--------------------------------------------------------------------------------#
# Polish the lineage: Loose the last rank which is already in 2nd column
#--------------------------------------------------------------------------------#
print "INFO: Polishing lineage\n";
open(OUT, ">", $outF) or die "Could not open output file $outF";
open(TMP, "<", $tmpOutF3) or die "Could not open temporary file $tmpOutF3";
while(<TMP>) {
	
	chomp($_);
	
	my @lineSplits = split("\t", $_, -1);
	my @lineageSplits = ("");
	@lineageSplits = split(";", $lineSplits[$#lineSplits]) if ($lineSplits[$#lineSplits]);
	
	# If the taxon in the second column is repeated in the lineage,
	# strip it from the lineage.
	if ($lineSplits[1] eq $lineageSplits[$#lineageSplits]) {
		splice(@lineageSplits, $#lineageSplits, 1)
	}
	
	print OUT join("\t", @lineSplits[0..1]), "\t", join(";", @lineageSplits), "\n"
}
close(TMP);
close(OUT);

# Status code 0 is error for system()
system("rm $tmpOutF3") and die "ERROR: Could not remove temporary output file";

print "\nDONE\n";
