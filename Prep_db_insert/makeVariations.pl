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
# Prepare a variation file for uORFdb
#======================================================================================#
# DESCRIPTION
#
# Read a file with annotated uORF variants and format it in a way that is compatible
# with the database.
#
#
# USAGE
#
# makeVariations.pl -input INFILE -blacklist BLACKLISTF
#
#	OR
#
# makeVariations.pl -i INFILE -b BLACKLISTF
#
# Read a TSV file with annotated uORF variants. Output a file in a database-friendly
# format in the directory of the input variant file (as out.INFILE). Please be aware
# of the limitations in the allele frequency calculation (see CAVEATS). Variants can
# be skipped based on transcript IDs in the BLACKLISTF. The version tag of transcripts
# is ignored.
#
#
# OUTPUT
#
# A tab-separated file (out.INFILE) containing variants split by allele
# with the following columns:
#
#	*) Variant position (0-based)
#	*) Allele (either alt or ref)
#	*) Column indicating, if the allele is a reference allele ("t") or not ("f").
#	*) Column indicating the length of the variation
#	*) Location of the variant, e.g. CDS (empty, if ref allele)
#	*) New CDS distance of the mutated uORF (empty, if ref allele)
#	*) New length of the mutated uORF (empty, if ref allele)
#	*) dbSNP ID matching position + allele (allele considered before splitting).
#		Empty, if ref allele.
#	*) ClinVar ID matching position + allele (allele considered before splitting).
#		Empty, if ref allele.
#	*) Frequency of allele in the cancers and optionally in the reference databases
#		(names: see following columns).
#	*) Frequency names of the same order and length as the previous column.
#	*) List of alt sequences caused by mutation in the uORF (empty, if ref allele)
#	*) Names of alt sequences of the same length and order as the previous column.
#	*) Effect of the alt sequences. Same length and order as the previous column.
#	*) List of patients affected by the variation.
#	*) List of cancers for the patients (same length and order as the previous column).
#	*) Chromosome of the reference uORF
#	*) NCBI ID of the transcript that the reference uORF resides on.
#	*) Start position of the reference uORF (0-based)
#	*) Stop position of the reference uORF (0-based)
#	*) Full nucleotide sequence of the reference uORF
#
#
# CAVEATS
#
# The input variants file must contain a header line and must only have one alt allele
# per line. Also, only one dbSNP ID and one ClinVar ID are allowed.
# The allele frequencies are inferred from the genotypes of affected patients and a
# total allele count which is hardcoded. See %counts hash.
# We use the sex of the patients to calculate the allele frequencies for the X and Y
# chromosome (with regard to the individual and the total allele counts). Again, the 
# total allele counts are hardcoded (see %countsX and %countsY hashes). It is assumed
# that the female genotype is XX and the male genotype is XY.
#--------------------------------------------------------------------------------------#


use strict;
use warnings;
use File::Basename;
use Getopt::Long;


#----------------------------------------------------------------#
# Functions
#----------------------------------------------------------------#
# Translate variant effects to more db friendly terms (non-redundant)
sub translateEffects {
	my %dict = (
		"aTIS->aTISpos"						=>	"changed position",
		"aTIS->uSTARTloss"					=>	"loss",
		"uAUG->uSTARTloss"					=>	"loss",
		"uSTOP->uSTOPdownstream"				=>	"downstream uSTOP",
		"NA->uSTOPinSeq"					=>	"upstream uSTOP",
		"uSTOP->uSTOPinSeq"					=>	"upstream uSTOP",
		"uSTOP->loss"						=>	"loss",
		
	);
	
	my $effect = $_[0];
	
	if (exists $dict{$effect}) {
		return $dict{$effect};
	}
	else {
		return $effect;
	}
}


#----------------------------------------------------------------#
# Hardcoded total allel counts
#----------------------------------------------------------------#
# Counts for autosomes
my %counts = (
			"TCGA-BRCA"	=> 117 * 2,
			"TCGA-COAD"	=> 113 * 2,
			"TCGA-LAML"	=> 48 * 2,
			"TCGA-LUAD"	=> 142 * 2,
			"TCGA-PRAD"	=> 122 * 2,
			"TCGA-SKCM"	=> 135 * 2,
			"Total"		=> 677 * 2
);

# Counts for chrX
my %countsX = (
			"TCGA-BRCA"	=> 233,
			"TCGA-COAD"	=> 170,
			"TCGA-LAML"	=> 69,
			"TCGA-LUAD"	=> 217,
			"TCGA-PRAD"	=> 122,
			"TCGA-SKCM"	=> 182,
			"Total"		=> 993
);

# Counts for chrY
my %countsY = (
			"TCGA-BRCA"	=> 1,
			"TCGA-COAD"	=> 56,
			"TCGA-LAML"	=> 27,
			"TCGA-LUAD"	=> 67,
			"TCGA-PRAD"	=> 122,
			"TCGA-SKCM"	=> 88,
			"Total"		=> 361,
);


#----------------------------------------------------------------#
# Read args from CLI
#----------------------------------------------------------------#
my $varF = "";
my $blacklistF = "";

GetOptions ("input=s" 		=> \$varF,
		"blacklist=s"	=> \$blacklistF
	) or die("ERROR: Unexpected argument\n");
			
die "ERROR: Missing essential argument -input\n" if (not $varF);

# Output file
my $outVarF = dirname($varF) . "/out." . basename($varF);


#----------------------------------------------------------------#
# Read blacklist (optional)
#----------------------------------------------------------------#
my %blacklist = ();
my $blacklistC = 0;

if ($blacklistF) {
	open(BLACKLIST, "<", $blacklistF) or die "ERROR: Cannot open blacklist file ->$blacklistF<-";
	while(<BLACKLIST>) {
		chomp($_);
		$_ =~ s/\.\d+$//;
		next if ($_ =~ m/^#/);
		$blacklist{$_} = undef;
	}
	close(BLACKLIST);
}

#----------------------------------------------------------------#
# Read variation file
#----------------------------------------------------------------#
open(OUTVAR, ">", $outVarF) or die "ERROR: Cannot open output variation file ->$outVarF<-";
open(VAR, "<", $varF) or die "ERROR: Cannot open variation file ->$varF<-";
VARIANT: while(<VAR>) {
	# Skip header
	next VARIANT if ($. < 2);
	chomp($_);
	
	my @splits = split("\t", $_);
	my ($pos, $ref, $alt, $cdsDist, $uorfLen, $startEffect, $stopEffect, $kozakEffect, $seqEffect, $exacs, $gnomads, 
		$topmedV2s, $topmedV3s, $tGTs, $cancers, $patients, 
		$dbsnp, $clinvar, $id, $location, $seq, $altseq, $altstart, $altstop, $altkozak, $altstrength, $sexes) = 
		@splits[3, 4, 5, 94, 92, 75, 76, 77, 78, 14, 15,
			16, 17, 19, 57, 58,
			9, 10, 73, 74, 79, 80, 82, 84, 86, 88, 98];
		
	die "ERROR: Only one alt allele allowed per line, but found ->$alt<-" if ($alt =~ m/[;,]/);
	
	die "ERROR: No uORF reference sequence" if ($seq eq ".");
	die "ERROR: No patient" if ($patients eq ".");
	die "ERROR: No cancer type" if ($cancers eq ".");
	
	# Variant position in input is 1-based. In the database, it is 0-based
	$pos--;
	
	# Use the longer allele for variant length
	my $length = length($ref);
	$length = length($alt) if (length($alt) > $length);


	#----------------------------------------------------------------#
	# Extract uORF metadata from variantion ID
	#----------------------------------------------------------------#
	if ($id =~ m/^.*(chr.+_.+_[0-9]+_[0-9]+)_.+$/) {
		$id = $1
	}
	else {
		die "ERROR: Unexpected ID ->$id<-";
	}	
	
	die "ERROR; No variant ID" if ($id eq "." or not $id);
	$id =~ s/NM_/NM-/;
	
	my ($chr, $ncbiid, $uStart, $uStop) = split("_", $id);
	$ncbiid =~ s/NM-/NM_/;
	
	# Filter variants whose transcripts appear on the blacklist.
	if (%blacklist) {
		my $tmpID = $ncbiid =~ s/\.\d+$//r;
		if (exists $blacklist{$tmpID}) {
			$blacklistC++;
			next;
		}
	}
	
	
	# Diamond is a placeholder for Indels that introduce a frame shift. They haven't been analyzed, yet.
	foreach my $value ($stopEffect, $kozakEffect, $altkozak, $altstrength) {
		if ($value eq "<>") {
			print "INFO: Skipping unannotated indel variant\n\t$_\n";
			next VARIANT 
		}
	}
	# In these columns, a diamond is illegal.
	foreach my $value ($cdsDist, $uorfLen, $startEffect, $seqEffect) {
		if ($value eq "<>") {
			die "ERROR: Unexpected diamond placeholder in column ->$value<- and line ->$.<-";
		}
	}
	

	
	#----------------------------------------------------------------#
	# Create lists of frequencies (cancer + refdbs) for ref and alt.
	# Create a name list.
	#----------------------------------------------------------------#	
	my %freqs = ($alt => {'afs' => "", 'names' => ""},
		$ref => {'afs' => "", 'names' => ""}); 
	
	$cancers =~ s/ //g;
	$tGTs =~ s/ //g;
	$sexes =~ s/ //g;
	
	my @cancers = split(";", $cancers);
	my @gts = split(";", $tGTs);
	my @sexes = split(";", $sexes);
	die "ERROR: Genotype list and cancer list must have the same length" if (@cancers != @gts);
	die "ERROR: Genotype list and sexes list must have the same length" if (@sexes != @gts);
	
	my %afs = ();
	
	# Use the genotype and the total allele counts from the cancer types
	# to infer the allele frequencies for each cancer type.
	# Use the information about the sex of the patient for the allele
	# frequencies on the X and Y chromosome.
	if (@cancers) {
		for (my $i = 0; $i <= $#cancers; $i++) {
			my $cancer = $cancers[$i];
			my $gt = $gts[$i];
			$gt =~ s/\/|\|//g;
			my $sex = $sexes[$i];
			my $af = 0;
			
			if ($sex ne "male" and $sex ne "female") {
				die "ERROR: Unknown sex ->$sex<-"
			}
			# May be a bug or due to pseudoautosomal regions
			if ($chr eq "chrY" and $sex eq "female") {
				print "WARNING: Found female patient with variant on chrY";
				$af = 0;
			}
			# Male patient on chrY or any patient on any other chromosome
			else {
				# Would not hold for females and chrY, but has been excluded before
				if ($gt eq "01") {
					$af = 1
				}
				# Would not hold for females and chrY, but has been excluded before
				elsif ($gt eq "10") {
					$af = 1
				}
				elsif ($gt eq "00") {
					$af = 0
				}
				elsif ($gt eq "11") {
					# Autosome
					if ($chr ne "chrX" and $chr ne "chrY") {
						$af = 2
					}
					elsif ($chr eq "chrX") {
						if ($sex eq "female") {
							$af = 2
						}
						# Male with 1/1 genotype: Bug or pseudoautosomal region
						else {
							print "WARNING: Genotype 1/1 for male patient on chrX\n";
							$af = 1
						}
					}
					# Would not hold for females, but females with chrY variants
					# have been excluded before.
					else {
						print "WARNING: Genotype 1/1 for male patient on chrY\n";
						$af = 1
					}
				}
				else {
					die "ERROR: Unrecognized genotype ->$gt<-";
				}
			}
			
			if (exists $afs{"Total"}) {
				# Autosome
				if ($chr ne "chrX" and $chr ne "chrY") {
					$afs{"Total"} += sprintf('%.6f', $af / $counts{"Total"})
				}
				elsif ($chr eq "chrX") {
					$afs{"Total"} += sprintf('%.6f', $af / $countsX{"Total"})
				}
				else {
					$afs{"Total"} += sprintf('%.6f', $af / $countsY{"Total"})
				}
			}
			else {
				# Autosome
				if ($chr ne "chrX" and $chr ne "chrY") {
					$afs{"Total"} = sprintf('%.6f', $af / $counts{"Total"})
				}
				elsif ($chr eq "chrX") {
					$afs{"Total"} = sprintf('%.6f', $af / $countsX{"Total"})
				}
				else {
					$afs{"Total"} = sprintf('%.6f', $af / $countsY{"Total"})
				}
			}
			
			
			if (exists $counts{$cancer}) {
				# Autosome
				if ($chr ne "chrX" and $chr ne "chrY") {
					$af = sprintf('%.6f', $af / $counts{$cancer});
				}
				elsif ($chr eq "chrX") {
					$af = sprintf('%.6f', $af / $countsX{$cancer});
				}
				else {
					$af = sprintf('%.6f', $af / $countsY{$cancer});
				}
			}
			else {
				die "ERROR: Unrecognized cancer type ->$cancer<-";
			}
			
			if (exists $afs{$cancer}) {
				$afs{$cancer} += $af
			}
			else {
				$afs{$cancer} = $af
			}
		}
		
		$freqs{$alt}->{'afs'} = join("; ", values(%afs)) . "; ";
		$freqs{$alt}->{'names'} = join("; ", keys(%afs)) . "; "
	}
	
	# Sometimes, I only have an older version of TopMed
	if ($topmedV3s eq ".") {
		$topmedV3s = $topmedV2s;
	}
	# The upstream script should only give me one of TopMed V2 and V3.
	# It should be the frequency that is given in dbSNP VCF.
	else {
		if ($topmedV2s ne ".") {
			die "ERROR: Only one of TopMed V2 and V3 is allowed, but found both.\n\t$_";
		}
	}
	
	my %refdbs = ("gnomAD" => $gnomads, "ExAC" => $exacs, "TopMed" => $topmedV3s);
	foreach my $refdb (keys(%refdbs)) {
		if ($refdbs{$refdb} ne ".") {
			my ($refc, $altc, $totalc) = split(",", $refdbs{$refdb});
			my $reffreq = sprintf('%.6f', $refc / $totalc);
			my $altfreq = sprintf('%.6f', $altc / $totalc);
			
			$freqs{$ref}->{'afs'} .= $reffreq . "; ";
			$freqs{$alt}->{'afs'} .= $altfreq . "; ";
			
			$freqs{$ref}->{'names'} .= $refdb . "; ";
			$freqs{$alt}->{'names'} .= $refdb . "; ";
		}
	}
	foreach my $allele ($ref, $alt) {
		$freqs{$allele}->{'names'} =~ s/; $//;
		$freqs{$allele}->{'afs'} =~ s/; $//;
	}

	
	#----------------------------------------------------------------#
	# Create lists of alt sequences for start codon, stop codon,
	# kozak context and fullseq. Create a names and effects list.
	#----------------------------------------------------------------#
	my %altSeqs = (	"start codon" 		=> [$startEffect, $altstart],
			"stop codon"		=> [$stopEffect, $altstop],
			"Kozak context"		=> [$kozakEffect, $altkozak],
			"sequence"			=> [$seqEffect, $altseq]);

	my $altNames = "";
	my $altSeqs = "";
	my $altEffects = "";
	
	foreach my $feature (keys(%altSeqs)) {
		my $effect = $altSeqs{$feature}->[0];
		$effect =~ s/_/ /g;
		$effect = translateEffects($effect);
		my $seq = $altSeqs{$feature}->[1];
		
		if ($effect eq ".") {
			if ($seq ne ".") {
				if ($altSeqs{"start codon"}->[0] !~ m/loss/) {
					die "ERROR: No variant effect, but sequence has changed, line ->$.<-"
				}
			}
			else {
				# Kozak should be loss, if start is lost
				if ($feature eq "Kozak context") {
					if ($altSeqs{"start codon"}->[0] =~ m/loss/) {
						$effect = "loss";
						$seq = "";
					}
					else {
						next;
					}
				}
				else {
					next;
				}
			}
		}
		else {
			if ($seq eq "." and $effect !~ m/pos/) {
				die "ERROR: Variant effect ->$effect<-, but sequence has not changed, line ->$.<-"
			}
		}
		$seq = "" if ($seq eq ".");
		$effect = "" if ($effect eq ".");
		
		# I don't like that sometimes, it says loss in the sequence
		if ($effect =~ m/loss/) {
			$seq = "";
		}
		if ($seq eq "loss") {
			$seq = "";
			$effect = "loss";
		}
		
		$altNames .= $feature . "; ";
		$altSeqs .= $seq . "; ";
		$altEffects .= $effect . "; "
	}
	$altNames =~ s/; $//;
	$altSeqs =~ s/; $//;
	$altEffects =~ s/; $//;
	
	
	#----------------------------------------------------------------#
	# Fix (empty) values
	#----------------------------------------------------------------#
	$cdsDist = "" if ($cdsDist eq "." or $cdsDist =~ m/loss/);
	$uorfLen = "" if ($uorfLen eq "." or $uorfLen =~ m/loss/);
	$dbsnp = "" if ($dbsnp eq ".");
	$clinvar = "" if ($clinvar eq ".");
	
	$patients =~ s/;$//;
	$patients =~ s/;/; /g;
	$cancers =~ s/;$//;
	$cancers =~ s/;/; /g;
	$dbsnp =~ s/;$//;
	$clinvar =~ s/;$//;
	$location =~ s/,$//;
	
	die "ERROR: Only one dbSNP ID allowed, but got ->$dbsnp<-" if ($dbsnp =~ m/;/);
	die "ERROR: Only one ClinVar ID allowed, but got ->$clinvar<-" if ($clinvar =~ m/;/);

	
	print OUTVAR $pos, "\t", $ref, "\t", "t", "\t", $length, "\t", "\t\t\t\t\t", $freqs{$ref}->{'afs'}, "\t", $freqs{$ref}->{'names'}, "\t\t\t\t", $patients, "\t", $cancers, "\t", $chr, "\t", $ncbiid, "\t", $uStart, "\t", $uStop, "\t", $seq, "\n";
	print OUTVAR $pos, "\t", $alt, "\t", "f", "\t", $length, "\t", $location, "\t", $cdsDist, "\t", $uorfLen, "\t", $dbsnp, "\t", $clinvar, "\t", $freqs{$alt}->{'afs'}, "\t", $freqs{$alt}->{'names'}, "\t", $altSeqs, "\t", $altNames, "\t", $altEffects, "\t", $patients, "\t", $cancers, "\t", $chr, "\t", $ncbiid, "\t", $uStart, "\t", $uStop, "\t", $seq, "\n";
}
close(VAR);
close(OUTVAR);

print "->$blacklistC<- entries skipped due to blacklist.\n" if ($blacklistC);
print "\nDONE\n";
