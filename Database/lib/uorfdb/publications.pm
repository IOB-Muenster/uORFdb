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


package publications;
use Exporter;
use DBI;
use Data::Dumper;
use Encode qw(decode encode);
use JSON::Tiny qw(decode_json);
use LWP::UserAgent;

# Own mods
use lib '../';
use uorfdb::db;

@ISA = ('Exporter');
@EXPORT      = ();
@EXPORT_OK   = qw(getFieldName fixBooleans checkEssential checkCount checkPubmed findDuplPubmed getMetaData getAuthors filterPubs);


# Returns the field name for an index
sub getFieldName {
	my $idx = $_[0];
	
	my @names = (
		"Version", #0
        "Taxon", #1
        "Taxonomy ID", #2
        "Gene Name", #3
        "Official Symbol", #4
        "Entrez Gene ID", #5
        "GenBank ID", #6
        "uORFs per Transcript", #7
        "Alternative Promoters", #8
        "Alternative Splicing", #9
        "Tissue Specific uORFs", #10
        "non-AUG uORFs", #11
        "Number", #12
        "Length", #13
        "Distance from 5'-CAP", #14
        "Distance From uORF-Stop to CDS", #15
        "CDS Overlap", #16
        "RNA Secondary Structure", #17
        "CDS Repression", #18
        "CDS Induction", #19
        "Start Site Selection", #20
        "Nonsense-Mediated Decay", #21
        "mRNA Destabilisation", #22
        "Ribosome Load", #23
        "Ribosome Pausing / Stalling", #24
        "Ribosome Shunting", #25
        "Kozak Consensus Sequence", #26
        "Translational Status", #27
        "Termination (context)", #28
        "uORF RNA/Peptide Sequence", #29
        "Regulatory Sequence Motif", #30
        "Co-factor / Ribosome Interaction", #31
        "Disease-related uORFs", #32
        "Acquired Mutations /SNPs", #33
        "Mouse Models", #34
        "Ribosome Profiling", #35
        "Bioinformatics / Array / Screens", #36
        "Proteomics", #37
        "Methods", #38
        "Reviews", #39
        "Pub Med ID", #40
        "Comment", #41
        "Year", #42
        "PDF", #43
        "Free Text" #44
	);
	
	if ($idx =~ m/\|/) {
		my @idxs = split(/\|/, $idx);
		my @res = ();
		
		foreach my $idx (@idxs) {
			if (exists $names[$idx]) {
				push(@res, $names[$idx]);
			}
			else {
				die "ERROR: No field with index $idx\n";
			}
		}
		return "->" . join(" or ", @res) . "<-";
	}
	else {
		if (exists $names[$idx]) {
			return "->" . $names[$idx] . "<-";
		}
		else{
			die "ERROR: No field with index $idx\n";
		}
	}
}


# Adapt boolean expressions to uORFdb2 format
sub fixBooleans {
	my $fieldCont = $_[0];
	my $fieldIdx = $_[1];
	my $errMsg = $_[2];
	
	# Indices of fields with boolean values
	my %booleans = map {$_ => undef} (8..39, 43);
	
	# The field is a boolean field
	if (exists $booleans{$fieldIdx}) {	
		if ($fieldCont =~ m/^yes$/i) {
			$fieldCont = 1;
		}
		elsif ($fieldCont =~ m/^no$/i) {
			$fieldCont = 0;
		}
		elsif ($fieldCont =~ m/^na$/i) {
			$fieldCont = "-";
		}
		elsif ($fieldCont =~ m/^$/i) {
			$fieldCont = "";
		}
		elsif ($fieldCont !~ m/^[01-]$/i) {
			$errMsg .= "ERROR: Column " . getFieldName($fieldIdx) .
				" marked as boolean, but found expression ->$fieldCont<-\n";
		}
	}
	
	return [$fieldCont, $errMsg];
}


# Check that all essential field combinations are there.
# It is OK, if all fields from a combination are missing.
sub checkEssential {
	my $fieldCont = $_[0];
	my $fieldIdx = $_[1];
	my @fields = @{$_[2]};
	my $errMsg = $_[3];
	
	# If a have a taxon [1], I want to have the taxon id [2] (and vice versa)
	# If I have a gene symbol [4], I want to have the gene id [5] or GenBank id [6] (and vice versa)
	# and the taxon [1] and its taxid [2]. BUT: It's OK to just have a taxon and its id
	# without any gene information.
	# If I have a GenBank ID [6] or a gene id [5], I need to know the taxonomy (taxid [2] + taxon [1])
	# and the symbol [4].
	my %missing = ( 1 => ['2'], 2 => ['1'], 4 => ['1', '2', '5|6'], 5 => ['1', '2', '4'], 6 => ['1', '2', '4']);
	
	if (exists $missing{$fieldIdx}) {
		foreach my $depField (@{$missing{$fieldIdx}}) {
			
			# If my key field is missing, I don't care about my dependency
			next if ($fieldCont =~ m/^$/);
			
			# 1|2 indicates that 1 or 2 or both must be present
			if ($depField =~ m/\|/) {
				my @depFields = split(/\|/, $depField);
				my $isErr = 1;
				foreach my $dep (@depFields) {
					my $depFieldCont = $fields[$dep];
					if ($depFieldCont !~ m/^$/) {
						$isErr = 0;
					}
				}
				$errMsg .= "ERROR: Dependency " . getFieldName($depField) . " of essential field " .
					getFieldName($fieldIdx) . " is empty\n" if($isErr == 1);
				
			}
			else {
				my $depFieldCont = $fields[$depField];
				if ($depFieldCont =~ m/^$/) {
					$errMsg .= "ERROR: Dependency " . getFieldName($depField) . " of essential field " .
						getFieldName($fieldIdx) . " is empty\n";
				}
			}
		}
	}
	
	return [$fieldCont, $errMsg];
}


# Check number of elements in dependent fields
sub checkCount {
	my @fields =  @{$_[0]};
	my $errMsg = $_[1];
	
	# Indeces of fields which should have an equal count of
	# elements
	my %equals = (1 => [2], 3 => [4, 5, 6], 4 => [5, 6], 5 => [6]);
	
	# I don't care, if one field is empty and the other is filled.
	# This is the task of checkEssential.
	# Non-empty fields should have the same number of elements
	# (separated by ";").
	foreach my $equal (keys(%equals)) {
		my $field = $fields[$equal];
		next if ($field =~ m/^$/);
		
		foreach my $dep (@{$equals{$equal}}) {
			my $depField = $fields[$dep];
			next if ($depField =~ m/^$/);
			
			if ((split("; ", $field)) != (split("; ", $depField))) {
				$errMsg .= "ERROR: Fields in columns " . getFieldName($dep) . " and " .getFieldName($equal) .
					" must have equal numbers of elements\n"
			}
		}
	}
	return $errMsg;
}


# One PubMed ID per line and valid ID (digits only)
sub checkPubmed {
	my $pubmedId = $_[0];
	my $errMsg = $_[1];
	
	if ($pubmedId =~ m/^\s*$/) {
		$errMsg .= "ERROR: PubMed ID must not be empty.\n"
	}
	else {
		$errMsg .= "ERROR: Invalid PubMed ID ->$pubmedId<-. Hint: Only one ID per line is allowed.\n" if($pubmedId !~ m/^\d+$/);
	}
	
	return $errMsg;
}


# Find entries with the same PubMed ID and raise error
sub findDuplPubmed {
	my $pubmedId = $_[0];
	my $pubmedsR = $_[1];
	my $errMsg = $_[2];
	
	return ($pubmedsR, $errMsg) if ($pubmedId =~ m/^\s*$/);
	
	if (exists $pubmedsR->{$pubmedId}) {
		$errMsg .= "ERROR: PubMed ID ->$pubmedId<- not unique\n";
	}
	else {
		$pubmedsR->{$pubmedId} = {'lines' => [[]]};
	}
	
	return($pubmedsR, $errMsg)
}


# Retrieve metadata for the publications from the NCBI Literature
# Citation Exporter API
sub getMetaData {
	my $pubsR = $_[0];
	my @queries = keys(%{$pubsR});
	
	# Only send 200 ids at once
	while(@queries) {
		my $query = join(",", splice(@queries, 0, 200));
		
		# Send to the API
		my $url = "https://api.ncbi.nlm.nih.gov/lit/ctxp/v1/pubmed/?format=ris";
		my $ua = LWP::UserAgent->new(timeout => 20);
		my $response = $ua->get($url."&id=$query");
		
		if (not $response->is_success) {
			print Dumper($response);
			die "ERROR: Cannot connect to website\n"
		}
		else {						
			if (exists $response->{'error'}) {
				die "API ERROR: " . $response->{'error'};
			}
			
			my $citation = $response->content;
			$citation = decode('UTF-8', $citation);
			
			my @citations = split("\n", $citation, -1);
			
			my $pubType = "";
			my $authors = "";
			my $title = "";
			my $year = "";
			my $abstract = "";
			my $startPage = "";
			my $endPage = "";
			my $volume = "";
			my $issue = "";
			my $id = "";
			my $doi = "";
			my $journalFull = "";
			my $journalAbbr = "";
			my $publisher = "";
			my $country = "";
			
			foreach my $line (@citations) {
				$line =~ s/\r$//;
				next if (not $line);
				
				my @words = split(" ", $line, -1);
				my @tags = splice(@words, 0, 2);
				my $tag = $tags[0] || undef;
				my $value = join(" ", @words);
				$value = encode('UTF-8', $value);
				
				die "ERROR: Line without RIS tag detected: ->" . $line . "<-\n" if (not defined $tag);
				
				# Skip keywords
				next if ($tag eq "KW");
				
				# Extract publication metadata
				if ($tag eq "TY") {
					$pubType = $value;
				}
				if ($tag eq "AU") {
					$authors .= "; " . $value;
				}
				elsif ($tag eq "T1") {
					$title = $value;
				}
				elsif ($tag eq "Y1") {
					$year = $value;
				}
				elsif ($tag eq "AB") {
					$abstract = $value;
				}
				elsif ($tag eq "SP") {
					$startPage = $value;
				}
				elsif ($tag eq "EP") {
					$endPage = $value;
				}
				elsif ($tag eq "VL") {
					$volume = $value;
				}
				elsif ($tag eq "IS") {
					$issue = $value;
				}
				elsif ($tag eq "AN") {
					$id = $value;
				}
				elsif ($tag eq "DO") {
					$doi = $value;
				}
				elsif ($tag eq "JF") {
					$journalFull = $value;
				}
				elsif ($tag eq "J2") {
					$journalAbbr = $value;
				}
				elsif ($tag eq "PB") {
					$publisher = $value;
				}
				elsif ($tag eq "CY") {
					$country = $value;
				}
				# This tag marks the end of the record.
				elsif ($line =~ m/^ER\b/) {
					
					die "ERROR: Publication type is mandatory" if (not $pubType);
					$authors =~ s/^; //;
					
					foreach my $lineR (@{$pubsR->{$id}->{'lines'}}) {
						next if (not @{$lineR});
						$lineR->[42] = $year;
						push(@{$lineR}, $pubType, $authors, $title, $abstract, $startPage, $endPage,
							$volume, $issue, $doi, $journalFull, $journalAbbr, $publisher, $country);
					}
					
					# Reset values for next publication	
					foreach my $var ($id, $pubType, $authors, $title, $year, $abstract, $startPage, $endPage,
						$volume, $issue, $doi, $journalFull, $journalAbbr, $publisher, $country) {
						$var = "";
					}
				}
			}
			
		# Sleep 1 sec due to API limits
		sleep(1);
		}
	}
	return ($pubsR);
}

# Get abbreviated author names for old uORFdb
sub getAuthors {
	my $pubsR = $_[0];
	my @queries = keys(%{$pubsR});
	
	my %authors = ();
	
	# Only send 500 at once (API limit)
	while(@queries) {
		my $query = join("&", splice(@queries, 0, 500));
		my %queries = ('id' => $query);
		
		# Send to the API
		my $url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json";
		my $ua = LWP::UserAgent->new(timeout => 20);
		my $response = $ua->post($url, \%queries);
		
		if (not $response->is_success) {
			print Dumper($response);
			die "ERROR: Cannot connect to website\n"
		}
		else {
			if ($response->content) {
				$response = decode_json ($response->content);
			}
			else {
				die "API ERROR: No or invalid response";
			}
			
			if (exists $response->{'error'}) {
				die "API ERROR: " . $response->{'error'};
			}
			
			#Loop over the returned PubMed ids
			my @returnedIds = @{$response->{'result'}->{'uids'}};
			
			# I cannot say that I always have 500 query elements. I just have not more than 500 elements.
			my $queryC = scalar((split("&", $query)));
			print "ERROR: API returned unexpected number of ids \n" if (@returnedIds != $queryC);
			
			foreach my $id (@returnedIds) {
				# May happen, if a queried PubMed ID starts with a 0, but is returned without 0.
				die "ERROR: Unexpected id ->$id<- returned by API" if (not exists $pubsR->{$id});
				
				# The PubMed ID may be invalid
				if (exists $response->{'result'}->{$id}->{'error'}) {
					die "ERROR: API error $response->{'result'}->{$id}->{'error'} for PubMed id $id\n";
				}
				
				my $authorsStr = "";
				
				# Get all authors and get UTF-8 right
				foreach my $author (@{$response->{'result'}->{$id}->{'authors'}}) {
					my $name = encode("UTF-8", $author->{'name'});
					$authorsStr .= $name . ", ";
				}
				
				$authorsStr =~ s/, $//;
				
				# *Team., -> *Team; as on PubMed website
				$authorsStr =~ s/\.\,/;/g;
				
				$authors{$id} = $authorsStr;
			}
		}
	}
	return \%authors;
}


# Retain only those publication entries, which are abscent from uORFdb
sub filterPubs {
	my $pubmedsR = $_[0];
	my $isDebug = $_[1];
	
	my $filterC = 0;
	
	# uORFdb field name in publication table => column index in line
	# Don't query year.
	my %fieldsOfInterest = (
		 "gene_name_paper"			=> 3,
		 "number"					=> 12,
		 "length"					=> 13,
		 "alt_promoters"			=> 8,
		 "alt_splicing"				=> 9,
		 "tissue_spec_uorfs"		=> 10,
		 "non_aug_uorfs"			=> 11,
		 "dist_5_cap"				=> 14,
		 "dist_uorf_stop_cds"		=> 15,
		 "cds_overlap"				=> 16,
		 "rna_sec_struct"			=> 17,
		 "cds_repression"			=> 18,
		 "cds_induction"			=> 19,
		 "start_site_sel"			=> 20,
		 "nonsense_mediated_decay"	=> 21,
		 "mrna_destab"				=> 22,
		 "ribosome_load"			=> 23,
		 "ribosome_pause_stall"		=> 24,
		 "ribosome_shunting"		=> 25,
		 "kozak_cons_seq"			=> 26,
		 "transl_status"			=> 27,
		 "term_context"				=> 28,
		 "uorf_rna_pept_seq"		=> 29,
		 "regul_seq_motif"			=> 30,
		 "cof_ribosome_interact"	=> 31,
		 "disease_rel_uorfs"		=> 32,
		 "acq_mutations_snps"		=> 33,
		 "mouse_models"				=> 34,
		 "ribosome_prof"			=> 35,
		 "bioinf_array_screens"		=> 36,
		 "proteomics"				=> 37,
		 "methods"					=> 38,
		 "reviews"					=> 39,
		 "id_publication"			=> 40,
	);
	
	my @ids = keys(%{$pubmedsR});
	
	# Connect to db
	my $dbh = uorfdb::db::connect();
		
	print "INFO: Filtering existing publications in uORFdb\n";
	ID: foreach my $id (@ids) {
			LINE: foreach my $lineR (@{$pubmedsR->{$id}{'lines'}}) {
				next LINE if (not @{$lineR});
				
				my $fields = "";
				my @values = ();
				foreach my $field (keys(%fieldsOfInterest)) {
					my $value = $lineR->[$fieldsOfInterest{$field}];
					if ($value =~ m/^\s*$/) {
						$fields .= "$field is NULL AND ";
					}
					else {
						$fields .= "$field = ? AND ";
						push (@values, $value);
					}
				}
								
				$fields =~ s/ AND $//;

				my $sth = $dbh->prepare('SELECT id_publication FROM content WHERE ' . $fields);
				$sth->execute(@values);
				my @selects = @{$sth->fetchall_arrayref};
				
				# Just empty the lines which are already in uORFdb; not the whole PubMed ID!
				foreach my $res (@selects) {
					if (exists $pubmedsR->{$res->[0]}) {
						print "DEBUG: Publication with PubMed ID ->" . $res->[0] ."<- already in database\n" if (defined $isDebug);
						$filterC++;
						$lineR = [];
					}
				}
			}
	}
	
	$dbh->disconnect or die $DBI::errstr;
	
	return $pubmedsR, $filterC;	
}

1;
