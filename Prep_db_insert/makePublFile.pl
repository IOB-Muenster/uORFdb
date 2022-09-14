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


#======================================================================================#
# Translate publication file format of original uORFdb to new file format
#--------------------------------------------------------------------------------------#
# DESCRIPTION
#
# Translate the tab-delimited publications file from the original uORFdb into a
# publication file for the updated version. Get metadata for genes and publications,
# check constraints and optionally skip database duplicates.
#
#
# USAGE
#
# makePublFile.pl -p PUBLF [-w WHITELIST] [-b BLACKLIST] [--no-check-db] [--debug]
#
# Process original uORFdb file PUBLF, perform minimal checks on lines in WHITELIST
# and skipp lines in BLACKLIST. See EXAMPLES for more details. By default,
# publications and genes in the new uORFdb database are skipped. If 
# --no-check-db is set, the program does not filter publications and genes which 
# already occur in the database (useful for a fresh insert).
#
#
# OUTPUT
#
# All output files are created in the directory of the PUBLF.
#
# Since constraints are checked and metadata is updated by querying the NCBI APIs,
# the script will output two files that can be used to update the original uORFdb:
#	pubs_old and authors.pm.
#
# The following files provide input data for the updated version of uORFdb:
#	pubs and pubs.genes.
#
# The following error files are created: pubs.err for errors in publications,
# entrez.err for errors in ENTREZ gene IDs and genbank.err for errors in
# GenBank accessions.
#
#
# EXAMPLES
#
# The script checks the database file for several logical constraints: Some fields
# in uORFdb provide information about the presence or absence of a topic in the
# individual publication (boolean fields). The scrip checks that the fields only
# have the values "yes" (or "1"), "no" or ("0"), "na" (or "-") or no content
# (case independent): boolean filter.
# The script checks a set of essential dependencies between the fields of the uORFdb1
# database file. An entry passes this dependency filter, if all fields of a
# dependency are provided or if all are missing. The fields and dependencies are:
# If there is a taxonomy ID, there must be a taxon and vice versa. If a gene symbol
# is provided, a taxon, a taxonomy ID and either a gene ID or a GenBank ID must be
# present. The same applies, if a GenBank ID or a gene ID are found.
# The following pairs must have the same numbers of entries, if they are not empty:
#		1) Taxon and taxonomy ID
#		2) gene name, gene symbol, gene ID and GenBank ID
#		3) gene symbol, gene ID and GenBank ID
#		4) gene ID and GenBank ID.
# 
# Each entry has to have a single PubMed ID that is a number. Duplicate PubMed IDs in
# the database file are reported as errors. The gene ID has to be a number or a list
# of numbers separated by semicolon. Any entry failing these filters is written to the
# separate output file "pubs.err".
#
# A whitelist file may be supplied to check only a very basic set of constraints for a
# publication (boolean filter and verifying that the PubMed ID is a single number).
# Lines in the whitelist starting with "#" are skipped and can be used for comments.
# The whitelist is useful, in case the algorithm is too strict for some gene names
# (e.g. names containing ","). Whitelisted lines are written to the output file after
# booleans have been fixed (!) and the PubMed ID has been checked. The whitelist file
# must contain one SHA256 sum per line. The SHA256 sums for each line include the "\n"
# and can be obtained from the error file. Likewise, a blacklist file can be defined.
# It has the same structure as the whitelist file. Blacklisted lines will be directly
# skipped after reading the line. Once entries in the error file have been edited, it
# is helpful to blacklist the old versions of the entries. Thus, if the old versions
# of the entries should ever appear again in the old version of uORFdb, they would be
# directly rejected by this script.
#
#
# CAVEATS
#
# If a whitelisted line fails the boolean adaptation or the check of the PubMed
# ID, it will still end up in the error file. This is intended, since all of these
# processes are essential for the db.
#
#
# DEPENDENCIES:
#
# DBI
# Digest::SHA
# JSON::Tiny
# LWP::UserAgent
# uorfdb::genes
# uorfdb::db
# uorfdb::publications
#--------------------------------------------------------------------------------------#


use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use Digest::SHA qw(sha256_hex);


# Own mods
use lib '../Database/lib/';
use uorfdb::genes qw(checkEntrezID filterGenes lookupGenbank lookupEntrez lookupAssembly getAssemblyByELink);
use uorfdb::publications qw(fixBooleans checkEssential checkCount checkPubmed findDuplPubmed getMetaData getAuthors filterPubs);


# Read the checksum list for white-/blacklisting
sub readList {
	my $inF = $_[0];
	my %checksums = %{$_[1]};
	
	open(LIST, "<", $inF) or die "Could not open checksums file, $inF";
	while(<LIST>) {
		# Skip comments
		next if ($_ =~ m/^#/);
		
		chomp($_);
		$checksums{$_} = undef;
	}
	close(LIST);
	
	return \%checksums;
}


my $publF = "";
my $whitelF = "";
my $blacklF = "";
my $isDebug = undef;
my $noCheckDB = undef;
GetOptions ("p=s" => \$publF,
			"w=s" => \$whitelF,
			"b=s" => \$blacklF,
			"debug" => \$isDebug,
			"no-check-db" => \$noCheckDB,
			) or die("ERROR: Unexpected argument\n");

die "ERROR: Missing essential argument -p\n" if (not $publF);

my $cwd = dirname($publF);
my $outF = $cwd . "/pubs";
my $errF = $cwd . "/pubs.err";
my $authorsF = $cwd . "/authors.pm";

my $errC = 0;
my $whitelistC = 0;
my $blacklistC = 0;
my $filterC = 0;

my %whitelist = ();
my %blacklist = ();
my $pubmedsR = {};
my $assembliesR = {};

#-----------------------------------------------------------------------------------------------------#
# Check constraints and filter offending entries
#-----------------------------------------------------------------------------------------------------#
# Read a list of whitelisted lines. For these lines, only the booleans will be fixed
# and the PubMed id will be checked => no downstream validations.
if ($whitelF) {
	%whitelist = %{readList($whitelF, \%whitelist)};
}

# Read a list of blacklisted lines. These lines will be removed.
if ($blacklF) {
	%blacklist = %{readList($blacklF, \%blacklist)};
}

open(ERROR, ">", $errF) or die "Could not open output file $errF";
open(PUB, "<", $publF) or die "Could not open publication file $publF";
while(<PUB>) {
	
	my $checksum = sha256_hex($_);
	
	# Skip blacklisted lines
	if (exists $blacklist{$checksum}) {
		$blacklistC++;
		next;
	}
	
	chomp($_);
	my $errMsg = "";
	my $isWhitelErr = 0;
	
	# Otherwise split looses empty elements at the end of the line	
	my @splits =split("\t", $_, -1);
	my $pubmedid = $splits[40];
	
	for (my $i = 0; $i <= $#splits; $i++) {
		($splits[$i], $errMsg) = @{publications::fixBooleans($splits[$i], $i, $errMsg)};
		
		# Skip further validation for whitelisted files
		if (not exists $whitelist{$checksum}) {
			($splits[$i], $errMsg) = @{publications::checkEssential($splits[$i], $i, \@splits, $errMsg)};
		}
	}
	
	$errMsg = publications::checkPubmed($pubmedid, $errMsg);
	
	# Skip further validation for whitelisted files
	if (not exists $whitelist{$checksum}) {
		($pubmedsR, $errMsg) = publications::findDuplPubmed($pubmedid, $pubmedsR, $errMsg);
		$errMsg = publications::checkCount(\@splits, $errMsg);
		$errMsg = genes::checkEntrezID($splits[5], $errMsg);
	}
	else {
		$whitelistC++;
		if (not exists $pubmedsR->{$pubmedid}) {
			$pubmedsR->{$pubmedid} = {'lines' => [[]]};
		}
	}
	
	if ($errMsg) {		
		$errC++;
		
		chomp($errMsg);
		
		if ($errMsg =~ m/boolean/i) {
			if (exists $whitelist{$checksum}) {
				print "WARNING: Cannot be whitelisted due to boolean error.\n"
			}
			
			$errMsg .= "\nWARNING: Cannot be whitelisted due to boolean error.";
			$isWhitelErr = 1;
		}
		# Error about (duplicate) PubMed Id cannot be whitelisted
		if ($errMsg =~ m/pubmed/i) {
			if (exists $whitelist{$checksum}) {
				print "WARNING: Cannot be whitelisted due to critical PubMed ID error.\n"
			}
			
			$errMsg .= "\nWARNING: Cannot be whitelisted due to critical PubMed ID error.";
			$isWhitelErr = 1;
		}
		
		
		print ERROR $_, "\n", $errMsg;
		
		if ($isWhitelErr ne 1) {
			print ERROR "\nTo whitelist use: $checksum"
		}
		
		print ERROR "\n----------------------------------------------------------------------------------------------------------------\n";
	}
	else {		
		# Add corrected line to %pubmeds
		if (not exists $pubmedsR->{$pubmedid}) {
			die "ERROR: Internal error. Unexpected PubMed ID ->$pubmedid<-"
		}
		else {
			if (not exists $pubmedsR->{$pubmedid}->{'lines'}) {
				$pubmedsR->{$pubmedid} = {'lines' => [\@splits]};
			}
			else {
				push (@{$pubmedsR->{$pubmedid}->{'lines'}}, \@splits)
			}
		}
	}	
}
close(PUB);
close(ERROR);


print "INFO: $errC records need manual intervention at $errF\n" if ($errC > 0);
print "INFO: Only minimal check performed on " . $whitelistC . " records, due to whitelist\n" if (%whitelist);
print "INFO: " . $blacklistC . " records removed, due to blacklist.\n";

if (keys(%{$pubmedsR})) {
	#-----------------------------------------------------------------------------------------------------#
	# Get publication date and authors for PASS publications
	#-----------------------------------------------------------------------------------------------------#
	print "INFO: Retrieving publication metadata.\n";
	
	# Hashes for gene information
	my %entrez = ();
	my %genbanks = ();
	
	# Update year and get authors
	$pubmedsR = publications::getMetaData($pubmedsR);
	my $authorsR = publications::getAuthors($pubmedsR);
	
	# check, if pubs exist in db
	if (defined $noCheckDB) {
		print "INFO: Did not check publications for database duplicates.\n"
	}
	else {
		($pubmedsR, $filterC) = publications::filterPubs($pubmedsR, $isDebug);
		print "INFO: " . $filterC . " records filtered as database duplicates\n";
	}
	
	# This file is needed for the new uORFdb
	open(OUT, ">", $outF) or die "Could not open output file $outF";
	
	# For compatibility with the old uORFdb, a second publication file is written.
	open(DBOLD, ">", $outF . "_old") or die "Could not open output file $outF" . "_old";
	
	# Authors is a file for old uORFdb that is in Perl module format
	open(AUTHORS, ">", $authorsF) or die "Could not open author output file $authorsF";
	print AUTHORS 'package authors;' . "\n" . 'our %values = (' . "\n";
	
	foreach my $pubmedid (keys(%{$pubmedsR})) {
		my @lines = @{$pubmedsR->{$pubmedid}->{'lines'}};
		my $wasPrinted = 0;
		my $out = "";
		my $outDBone = "";
	
		for (my $i = 0; $i <= $#lines; $i++) {
			
			my $lineR = $lines[$i];
			
			# If @$lineR is empty, the PubMed ID only exists in the error output or is a
			# database duplicate --> skip here
			next if (not @{$lineR});
			$out .= join("\t", @{$lineR});
			
			my $lastLineIdx = $#{$lineR};
			$outDBone .= join("\t", @{$lineR}[0..$lastLineIdx-16]);
		
			# Just publication year, not date
			$outDBone .= "\t" . (split("/", @{$lineR}[$lastLineIdx-15]))[0];
			
			$outDBone .= "\t" . join("\t", @{$lineR}[$lastLineIdx-14..$lastLineIdx-13]);
			
			
			# Get Entrez and or GenBank ID for publications which are not in the
			# error output
			if ($lineR->[5] !~ m/^\s*$/) {
				foreach my $id (split(";", $lineR->[5])) {
					$id =~ s/\s+//g;
					next if ($id !~ m/^\d+$/);
					$entrez{$id} = []
				}
			}
			if ($lineR->[6] !~ m/^\s*$/) {
				foreach my $id (split(";", $lineR->[6])) {
					$id =~ s/\s+//g;
					$genbanks{$id} = []
				}
			}
			
			# Needed for old uORFdb.
			# Print authors only for records that had no errors.
			if ($wasPrinted == 0) {
				print AUTHORS "\t", $pubmedid, " => \"", $authorsR->{$pubmedid}, "\",\n";
				$wasPrinted = 1
			}
			
			$out .= "\n";
			$outDBone .= "\n";
		}
		print OUT $out;
		
		# Output with a format compatible with old uORFdb
		print DBOLD $outDBone;
	}
		
	print AUTHORS ');' . "\n" . '1;';
	
	close(AUTHORS);
	close(OUT);
	close(DBOLD);

	# Sort publication output. Return code 0 is error for system.
	system('sort -t ' ."'\t' -k43,43nr -k41,41nr $outF > $outF.sorted") and die "ERROR: Could not sort publications";
	system("mv $outF.sorted $outF") and die "ERROR: Could not write sorted publications";
	system('sort -t ' ."'\t' -k43,43nr -k41,41nr $outF" . "_old > $outF"."_old.sorted") and die "ERROR: Could not sort publications";
	system("mv $outF" ."_old.sorted $outF" . "_old") and die "ERROR: Could not write sorted publications";

	print "DEBUG: Stored " . keys(%{$pubmedsR}) . " PubMed IDs \n" if ($isDebug);
	print "DEBUG: Stored " . (keys(%entrez) + keys(%genbanks)) . " gene accessions \n" if ($isDebug);
	
	#-----------------------------------------------------------------------------------------------------#
	# Get gene information for PASS publications
	#-----------------------------------------------------------------------------------------------------#
	print "INFO: Retrieving gene metadata.\n";
	
	my $entrezR = \%entrez;
	my $genbanksR = \%genbanks;
	
	# check, if genes exist in db
	if (defined $noCheckDB) {
		print "INFO: Did not check genes for database duplicates.\n"
	}
	else {
		# I don't want to get gene infos on genes which are already in the db.
		# -> Remove those.
		($entrezR, $genbanksR, $filterC) = genes::filterGenes($entrezR, $genbanksR);
		print "INFO: " . $filterC . " genes filtered as database duplicates\n";
	}
	
	# Get infos on missing genes
	($genbanksR, $assembliesR) = genes::lookupGenbank($genbanksR, $assembliesR);
	($entrezR, $assembliesR) = genes::lookupEntrez($entrezR, $assembliesR);
	$assembliesR = genes::lookupAssembly($assembliesR);
	
	
	$outF = $cwd . "/pubs.genes";
	$errF = $cwd . "/genbank.err";
	open(OUT, ">", $outF);
	open(ERROR, ">", $errF);
	foreach my $id (keys(%{$genbanksR})) {
		if (@{$genbanksR->{$id}}) {		
			my $assembly = $genbanksR->{$id}->[1] || "";
			
			if ($assembly !~ m/^\s*$/) {
				if (exists $assembliesR->{$assembly}) {
					$assembly = $assembliesR->{$assembly} || "";
					$genbanksR->{$id}->[1] = $assembly;
				}
				else {
					# Unresolved assembly accesion --> ERROR.
					die "ERROR: Unexpected assembly accession ->$assembly<-" if ($assembly =~ m/^GC[F_A]_[0-9]+\.[0-9]+$/)
				}
			}
			
			print OUT join("\t", @{$genbanksR->{$id}}), "\t\t", $id, "\n";
		}
		else {
			print ERROR $id, "\n";
		}
	}
	close(ERROR);
	
	$errF = $cwd . "/entrez.err";
	open(ERROR, ">", $errF);
	foreach my $id (keys(%{$entrezR})) {
		if (@{$entrezR->{$id}}) {
			my $assembly = $entrezR->{$id}->[1] || "";
			
			if ($assembly !~ m/^\s*$/) {
				if (exists $assembliesR->{$assembly}) {
					$assembly = $assembliesR->{$assembly};
					$entrezR->{$id}->[1] = $assembly;
				}
				else {
					# Unresolved assembly accesion --> ERROR.
					die "ERROR: Unexpected assembly accession ->$assembly<-" if ($assembly =~ m/^GC[F_A]_[0-9]+\.[0-9]+$/)
				}
			}
			print OUT join("\t", @{$entrezR->{$id}}), "\t", $id, "\t\n";
		}
		else {
			print ERROR $id, "\n";
		}
	}
	close(ERROR);
	close(OUT);
}
else {
	print "INFO: No records remaining\n"
}

print "DONE"
