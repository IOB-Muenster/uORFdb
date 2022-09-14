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


package genes;
use Exporter;
use DBI;
use Data::Dumper;
use JSON::Tiny qw(decode_json);
use LWP::UserAgent;

# Own mods
use lib '../';
use uorfdb::db;

@ISA = ('Exporter');
@EXPORT      = ();
@EXPORT_OK   = qw(checkEntrezID filterGenes lookupGenbank lookupEntrez lookupAssembly getAssemblyByELink);


sub checkEntrezID {
	my $id = $_[0];
	my $errMsg = $_[1];
	
	if ($id =~ m/^\s*$/) {
		return $errMsg;
	}
	else {
		foreach my $id (split(";", $id)) {
			$id =~ s/ //g;
			# "#" is a placeholder for unknown
			if ($id !~ m/^\d+$/ and $id ne "#") {
				$errMsg .= "ERROR: Entrez gene ID ->$id<- appears to be invalid."
			}
		}
	}
	return $errMsg;
}


# Retain only Entrez gene or GenBank IDs which are abscent from uORFdb
sub filterGenes {
	my $entrezR = $_[0];
	my $genbanksR = $_[1];
	
	my @entrez = keys(%{$entrezR});
	my @genbanks = keys(%{$genbanksR});
	
	my $filterC = 0;
	
	# Connect to db
	my $dbh = uorfdb::db::connect();	
	my $sth = "";
	
	if (@entrez) {
		print "INFO: Getting missing Entrez gene IDs from the db\n";
		my $sth = $dbh->prepare('SELECT geneid FROM gene WHERE geneid IN (?' . (',?' x (scalar(@entrez) - 1)) . ')');
		$sth->execute(@entrez);
		
		my @selects = @{$sth->fetchall_arrayref};
		foreach my $res (@selects) {
			if (exists $entrezR->{$res->[0]}) {
				delete $entrezR->{$res->[0]};
				$filterC++;
			}
			# If the id does not exist, it was already deleted.
			# This means, that the ID lives twice inside the db.
			# However, gene IDs must be unique:
			# https://doi.org/10.1093/nar/gkq1237
			else {
				die "ERROR: ID ->" . $res->[0] . "<- is not unique.";
			}
		}
	}
	
	if (@genbanks) {
		print "INFO: Getting missing Nucleotide accessions from the db\n";
		$sth = $dbh->prepare('SELECT geneid FROM gene WHERE genbank_id IN (?' . (',?' x (scalar(@genbanks) - 1)) . ')');
		$sth->execute(@genbanks);
		
		my @selects = @{$sth->fetchall_arrayref};
		foreach my $res (@selects) {
			if (exists $genbanksR->{$res->[0]}) {
				delete $genbanksR->{$res->[0]};
				$filterC++;
			}
			# If the id does not exist, it was already deleted.
			# This means, that the ID lives twice inside the db.
			# However, the genbank id aka the NCBI nucleotide
			# accession must also be unique in the db:
			# https://www.ncbi.nlm.nih.gov/genbank/samplerecord/#AccessionB
			else {
				die "ERROR: ID ->" . $res->[0] . "<- is not unique."
			}
		}
	}
	
	$dbh->disconnect or die $DBI::errstr;
	
	return ($entrezR, $genbanksR, $filterC);	
}


# Get gene metadata by NCBI accession ( = GenBank ID)
sub lookupGenbank {
	my $genbanksR = $_[0];
	my $assembliesR = $_[1];	
	
	my @queries = keys(%{$genbanksR});
	my $url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&retmode=json";
	
	# Only send 500 at once (API limit)
	while(@queries) {
		my $query = join("&", splice(@queries, 0, 500));
		my %queries = ('id' => $query);
		
		# Send to the API
		my $ua = LWP::UserAgent->new(timeout => 20);
		my $response = $ua->post($url, \%queries);
		
		if (not $response->is_success) {
			print Dumper($response);
			die "ERROR: Cannot connect to website"
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
			
			#Loop over the returned ids
			for (my $i = 0; $i <= $#{$response->{'result'}->{'uids'}}; $i ++) {
				
				my $id = $response->{'result'}->{'uids'}->[$i];
				my $res = $response->{'result'}->{$id};
				
				if (exists $res->{'error'}) {
					print "API ERROR: $res->{'error'} for GenBank id $id\n";
					$genbanksR->{$id} = [];
					next;
				}
				my $taxid = $res->{'taxid'} || "";
				my $assembly = $res->{'assemblyacc'} || "";
				my $names = $res->{'title'} || "";
				
				# Subs provide further information.
				# Subtype is a list of information types (=~ hash key).
				# Subname is a list of data for the types (=~ hash value).
				my $subtype = $res->{'subtype'} || "";
				my $subname = $res->{'subname'} || "";
				
				my $chrom = "";
				
				if ($subtype and $subname) {
					my @types = split(/\|/, $subtype);
					my @names = split(/\|/, $subname);

					foreach (my $i = 0; $i <= $#types; $i++) {
						if ($types[$i] =~ m/^chromosome$/i) {
							$chrom = $names[$i];
						}
					}
				}
				
				
				# Genbank ids in the publication file may or may not have a version tag
				my $genbankId = $res->{'accessionversion'} || "";
				
				if ($assembly !~ m/^\s*$/) {
					$assembly = getAssemblyByELink($assembly);
					
					if ($assembly =~ m/^GC[F,A]_[0-9]+\.[0-9]+$/) {
						$assembliesR->{$assembly} = undef
					}
				}
				else {
					$assembly = getAssemblyByELink($genbankId) if ($genbankId !~ m/^\s*$/);
					
					if ($assembly =~ m/^GC[F,A]_[0-9]+\.[0-9]+$/) {
						$assembliesR->{$assembly} = undef
					}
				}
				
				if (exists $genbanksR->{$genbankId}) {
					push(@{$genbanksR->{$genbankId}}, ($taxid, $assembly, $chrom, "", $names, ""));
				}
				else {
					$genbankId = $res->{'caption'} || "";
					if (exists $genbanksR->{$genbankId}) {
						push(@{$genbanksR->{$genbankId}}, ($taxid, $assembly, $chrom, "", $names, ""));
					}
					else {
						die "ERROR: Could not get GenBank id from uid ->$id<-";
					}
					
				}
			}
		}
		# Sleep 700 msec due to API limits
		select(undef, undef, undef, 0.7);
	}
	
	return ($genbanksR, $assembliesR);
}


# Get gene metadata by Entrez gene ID
sub lookupEntrez {
	my $entrezR = $_[0];
	my $assembliesR = $_[1];	
	
	
	my @queries = keys(%{$entrezR});
	my $url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&retmode=json";
	
	# Only send 500 at once (API limit)
	while(@queries) {
		my $query = join("&", splice(@queries, 0, 500));		
		my %queries = ('id' => $query);
		
		# Send to the API
		my $ua = LWP::UserAgent->new(timeout => 20);
		my $response = $ua->post($url, \%queries);
		
		if (not $response->is_success) {
			print Dumper($response);
			die "ERROR: Cannot connect to website"
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
			
			#Loop over the ids
			for (my $i = 0; $i <= $#{$response->{'result'}->{'uids'}}; $i ++) {
				
				my $id = $response->{'result'}->{'uids'}->[$i];
				my $res = $response->{'result'}->{$id};
				
				if (exists $res->{'error'}) {
					print "API ERROR: $res->{'error'} for Entrez gene id $id\n";
					$entrezR->{$id} = [];
					next;
				}
				
				my $taxid = $res->{'organism'}->{'taxid'} || "" ;
				
				# This appears to always give me the latest assembly accession
				my $assembly = $res->{'locationhist'}->[0]->{'assemblyaccver'} || "";
				
				
				# Sometimes, the gene database does not provide an assembly accession
				# in the locationhist. However, it gives a nuccore id as
				# {'genomicinfo'}->[0]->{'chraccver'}. I have to use
				# a custom elink workaround to lookup the nuccore id.
				# This directly returns the assembly name
				if (not $assembly) {
					# It's OK to have multiple locations, as long as the chraccver
					# is the same.
					my %locations = ();
					
					foreach my $loc (@{$res->{'genomicinfo'}}) {
						$locations{$loc->{'chraccver'}} = undef;
					}
					die "ERROR: Ambiguous locations for one gene ->$id<-" if (keys(%locations) > 1);
					
					$assembly = $res->{'genomicinfo'}->[0]->{'chraccver'} || "";
					
					if ($assembly !~ m/^\s*$/) {
						$assembly = getAssemblyByELink($assembly);			
					}
					if ($assembly =~ m/^GC[F,A]_[0-9]+\.[0-9]+$/) {
						$assembliesR->{$assembly} = undef
					}
				}
				else {
					$assembliesR->{$assembly} = undef if ($assembly !~ m/^\s*$/);
				}
				
				my $chrom = $res->{'chromosome'} || "";
				$chrom = "chr".$chrom if($chrom !~ m/^chr/ and $taxid == 9606);
				
				my $offSymb = $res->{'nomenclaturesymbol'} || "";
				$offSymb = $res->{'name'} || "" if (not $offSymb);
				
				$offSymb =~ s/^ //;
				$offSymb =~ s/ $//;
				die "ERROR: Missing symbol for geneid ->$id<-" if (not $offSymb);
				
				
				# I want to get all names and the description is fine, if
				# it does not already appear in the official names
				my $offName = $res->{'nomenclaturename'} || "";
				my $otherNames = $res->{'otherdesignations'}|| "" ;
				my $descript = $res->{'description'} || "";
				my $descrEscaped = quotemeta($descript);
				my $offNescaped = quotemeta($offName);
				
				if ($otherNames) {		
					if ($otherNames =~ m/^([^|]+\|)*$offNescaped(\|[^|]+)*$/i or not $offName) {
						$offName = $otherNames
					}
					else {
						$offName .= "|$otherNames";
					}
				}
		
				if ($descript) {
					if (not $offName) {
						$offName = $descript;
					}
					elsif ($offName !~ m/^([^|]+\|)*$descrEscaped(\|[^|]+)*$/i) {
						$offName .= "|$descript";
					}
				}
				
				$offName =~ s/\|/#/g;
				$offName =~ s/^#+//;
				$offName =~ s/^ //;
				$offName =~ s/ $//;
				
				my $aliases = $res->{'otheraliases'} || "";
				$aliases =~ s/, /#/g;
				$aliases =~ s/^ //;
				$aliases =~ s/ $//;
				
				# Remove the symbol from the aliases, if necessary
				my $offSescaped = quotemeta($offSymb);
				$aliases =~ s/^$offSescaped#//;
				$aliases =~ s/#$offSescaped$//;
				$aliases =~ s/#$offSescaped#/#/;
				$aliases =~ s/^$offSescaped$//;
				
				push (@{$entrezR->{$id}},($taxid, $assembly, $chrom, $offSymb, $offName, $aliases));
			}
		}
		# Sleep 700 msec due to API limits
		select(undef, undef, undef, 0.7);
	}
	
	return ($entrezR, $assembliesR);
}


# Get assembly name from from an assembly accession
sub lookupAssembly {
	my $assembliesR = $_[0];

	my @queries = map {$_ .= "[Assembly Accession]"} keys(%{$assembliesR});
	
	# The assembly accession is unique, so no need to use the
	# taxid: https://www.ncbi.nlm.nih.gov/assembly/model/
	my $url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&retmode=json";
	
	my @ids = ();
	
	# Only send 10 at once (API limit)
	while(@queries) {
		my $query = join(" OR ", splice(@queries, 0, 10));		
		
		my %queries = ('term' => $query);
		
		# Send to the API
		my $ua = LWP::UserAgent->new(timeout => 20);
		my $response = $ua->post($url, \%queries);
		
		if (not $response->is_success) {
			print Dumper($response);
			die "ERROR: Cannot connect to website"
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
			if (exists $response->{'esearchresult'}->{'error'}) {
				die "API error with assemblies ->$query<-";
			}
			
			# Getting the ids for all assembly accessions at once has the
			# negative effect, that I have no idea which id belongs to
			# which assembly accession. This, however, will be resolved
			# downstream.
			push(@ids, @{$response->{'esearchresult'}->{'idlist'}});
		}
		
		# Sleep 700 msec due to API limits
		select(undef, undef, undef, 0.7);
	}
	
	# Sleep 700 msec due to API limits
	select(undef, undef, undef, 0.7);
	
	# Get the assembly names and assembly accessions and accessions
	# from the ids from the previous step.
	$url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&retmode=json";
	@queries = @ids;
	
	# Only send 500 at once (API limit)
	while(@queries) {
		my $query = join("&", splice(@queries, 0, 500));		
		
		my %queries = ('id' => $query);
		
		# Send to the API
		my $ua = LWP::UserAgent->new(timeout => 20);
		my $response = $ua->post($url, \%queries);
		
		if (not $response->is_success) {
			print Dumper($response);
			die "ERROR: Cannot connect to website"
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
			
			#Loop over the returned ids
			for (my $i = 0; $i <= $#{$response->{'result'}->{'uids'}}; $i ++) {
				
				my $id = $response->{'result'}->{'uids'}->[$i];
				my $res = $response->{'result'}->{$id};
				
				if (exists $res->{'error'}) {
					die "API ERROR: $response->{'result'}->{$id}->{'error'} for assembly id ->$id<-";
					next;
				}
				
				my $assemblyN = $res->{'assemblyname'} || "";
				my $assemblyAcc = $res->{'assemblyaccession'} || "";
				
				if (not $assemblyN or not $assemblyAcc) {
					die "API ERROR: Missing data for for assembly id ->$id<-"
				}
				
				if (exists $assembliesR->{$assemblyAcc}) {
					$assembliesR->{$assemblyAcc} = $assemblyN;
				}
				else {
					die "API ERROR: Unexpected assembly accession ->$assemblyAcc<- returned";
				}
			}
		}
		# Sleep 700 msec due to API limits
		select(undef, undef, undef, 0.7);
	}
	return $assembliesR;
}


# Get assembly via elink
sub getAssemblyByELink {
	my $assembly = $_[0];
	
	#-----------------------------------------------------#
	# Get the assembly id via elink
	#-----------------------------------------------------#
	my $url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nuccore&db=assembly&retmode=json&id=" . $assembly;
	
	# Send to the API
	my $ua = LWP::UserAgent->new(timeout => 20);
	my $response = $ua->get($url);
	
	if (not $response->is_success) {
		print Dumper($response);
		die "ERROR: Cannot connect to website or invalid ID"
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
		
		my @res = ();
		if (exists $response->{'linksets'}->[0]->{'linksetdbs'}) {
			@res = @{$response->{'linksets'}->[0]->{'linksetdbs'}->[0]->{'links'}};
		}
		
		if (@res) {
			$assembly = $res[0];
		}
		else {
			return "";
		}
	}
	# Sleep 700 msec due to API limits
	select(undef, undef, undef, 0.7);
	
	
	#-----------------------------------------------------#
	# Get the assembly name from the assembly id.
	#-----------------------------------------------------#
	$url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&retmode=json&id=" . $assembly;
	
	# Send to the API
	$ua = LWP::UserAgent->new(timeout => 20);
	$response = $ua->get($url);
	
	if (not $response->is_success) {
		print Dumper($response);
		die "ERROR: Cannot connect to website"
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
		
		die "ERROR: Ambiguous assembly id ->$assembly<-" if (@{$response->{'result'}->{'uids'}} > 1);
		
		my $id = $response->{'result'}->{'uids'}->[0];
		my $res = $response->{'result'}->{$id};
			
		if (exists $res->{'error'}) {
			die "API ERROR: $res->{'error'} for assembly id ->$id<-";
			next;
		}
		
		# There will be an entry for the 'latestaccession', if the current assembly
		# is not the latest.
		if ($res->{'latestaccession'}) {
			die "ERROR: The current assembly id ->$id<- is not the latest for accession ->$assembly<-";
		}
		
		my $assemblyN = $res->{'assemblyname'} || "";
		my $assemblyAcc = $res->{'assemblyaccession'} || "";
		
		if (not $assemblyN or not $assemblyAcc) {
			die "API ERROR: Missing data for for assembly id ->$id<-"
		}
		return $assemblyN
	}
}


1;
