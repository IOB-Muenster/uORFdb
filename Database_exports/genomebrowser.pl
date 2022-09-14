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


use strict;

use Sys::Hostname;
use DBI;
use CGI;

use lib '../Database/lib/';
use uorfdb::db;

# Get gene, transcript and uORF data for one or multiple id_gene(s)
# Note: Exons may not overlap
sub getGTUdata {
	my %ids = %{$_[0]};
	my $dbh = $_[1];

	my %res            = ();
	my $globalTax      = "";
	my $globalAssembly = "";

	my $sids = join(",", keys(%ids));

	# Get results from db as arrayRef of rows
	# LEFT OUTER JOIN does not loose information for rows which don't appear in all tables.
	# ASSUMPTION: If an element is present in a table, I can find all the necessary information for the custom tracks
	my $resultsR =
	  $dbh->selectall_arrayref("select g.id, t.id, g.id_taxonomy, g.assembly, g.symbol, "
		  . "g.chrom, t.strand, t.ncbiid, t.pos, t.length, t.cds_start, t.cds_end, e.startpos, "
		  . "e.endpos, u.id, CONCAT(t.ncbiid, '_', u.startcodon, '.', CASE WHEN t.strand = '+' THEN "
		  . "(select count(*)+1 from uorf u2 where u2.id_transcript = u.id_transcript and u2.startcodon "
		  . "= u.startcodon and u2.startpos < u.startpos) WHEN t.strand = '-' THEN "
		  . "(select count(*)+1 from uorf u2 where u2.id_transcript = u.id_transcript and u2.startcodon = "
		  . "u.startcodon and u2.endpos > u.endpos) END)"
		  . "as name,  u.startpos, u.endpos, u.rframe FROM gene g LEFT OUTER JOIN transcript t "
		  . "ON g.id = t.id_gene LEFT OUTER JOIN exon e ON e.id_transcript = t.id LEFT OUTER JOIN uorf u "
		  . "ON u.id_transcript = t.id where g.id in ($sids)");

	foreach my $rowR (@{$resultsR}) {
		my (
			$idGene,  $idTranscr, $taxonomy,     $assembly,   $symbol,    $chrom,
			$strand,  $ncbiid,    $transcrStart, $transcrLen, $cdsStart,  $cdsEnd,
			$exStart, $exEnd,     $idUORF,       $uorfName,   $uorfStart, $uorfEnd,
			$uframe
		) = @{$rowR};

		$globalTax = $taxonomy                           if (not $globalTax);
		die "Genes must only have one source organism" if ($taxonomy ne $globalTax);

		$globalAssembly = $assembly                      if (not $globalAssembly);
		die "Genes must only have one source assembly" if ($assembly ne $globalAssembly);

		# A gene without a transcript
		if (not defined $idTranscr) {
			die "Cannot create track for gene $idGene. No transcripts.";
		}

		# A gene with a transcript, but without exons and uorfs
		elsif (defined $idTranscr and not defined $exStart and not defined $uorfStart) {
			my $transcrEnd = $transcrStart + $transcrLen;

			if (exists $res{$idGene}) {
				if (exists $res{$idGene}{$idTranscr}) {
					next if (exists $res{$idGene}{$idTranscr}{'data'});

					$res{$idGene}{$idTranscr}{'data'} = {
						'name'    => $ncbiid,
						'chr'     => $chrom,
						'gsymbol' => $symbol,
						'strand'  => $strand,
						'start'   => $transcrStart,
						'end'     => $transcrEnd,
						'cds'     => [$cdsStart, $cdsEnd]
					};
				}
				else {
					$res{$idGene}{$idTranscr} = {
						'data' => {
							'name'    => $ncbiid,
							'chr'     => $chrom,
							'gsymbol' => $symbol,
							'strand'  => $strand,
							'start'   => $transcrStart,
							'end'     => $transcrEnd,
							'cds'     => [$cdsStart, $cdsEnd]
						}
					};
				}
			}
			else {
				$res{$idGene} = {
					$idTranscr => {
						'data' => {
							'name'    => $ncbiid,
							'chr'     => $chrom,
							'gsymbol' => $symbol,
							'strand'  => $strand,
							'start'   => $transcrStart,
							'end'     => $transcrEnd,
							'cds'     => [$cdsStart, $cdsEnd]
						}
					}
				};
			}
		}

		# A gene with a transcript and uorfs, but without exons
		elsif (defined $idTranscr and not defined $exStart and defined $uorfStart) {
			my $transcrEnd = $transcrStart + $transcrLen;

			if (exists $res{$idGene}) {
				if (exists $res{$idGene}{$idTranscr}) {
					next if (exists $res{$idGene}{$idTranscr}{'data'});

					$res{$idGene}{$idTranscr}{'data'} = {
						'name'    => $ncbiid,
						'chr'     => $chrom,
						'gsymbol' => $symbol,
						'strand'  => $strand,
						'start'   => $transcrStart,
						'end'     => $transcrEnd,
						'cds'     => [$cdsStart, $cdsEnd]
					};

					if (exists $res{$idGene}{$idTranscr}{'uorf'}) {
						next if (exists $res{$idGene}{$idTranscr}{'uorf'}{$idUORF});
						$res{$idGene}{$idTranscr}{'uorf'}{$idUORF} = [$uorfStart, $uorfEnd, $uframe, $uorfName];
					}
					else {
						$res{$idGene}{$idTranscr}{'uorf'} = {$idUORF => [$uorfStart, $uorfEnd, $uframe, $uorfName]};
					}
				}
				else {
					$res{$idGene}{$idTranscr} = {
						'data' => {
							'name'    => $ncbiid,
							'chr'     => $chrom,
							'gsymbol' => $symbol,
							'strand'  => $strand,
							'start'   => $transcrStart,
							'end'     => $transcrEnd,
							'cds'     => [$cdsStart, $cdsEnd]
						},
						'uorf' => {$idUORF => [$uorfStart, $uorfEnd, $uframe, $uorfName]}
					};
				}
			}
			else {
				$res{$idGene} = {
					$idTranscr => {
						'data' => {
							'name'    => $ncbiid,
							'chr'     => $chrom,
							'gsymbol' => $symbol,
							'strand'  => $strand,
							'start'   => $transcrStart,
							'end'     => $transcrEnd,
							'cds'     => [$cdsStart, $cdsEnd]
						},
						'uorf' => {$idUORF => [$uorfStart, $uorfEnd, $uframe, $uorfName]}
					}
				};
			}
		}

		# A gene with a transcript with exons, but without uorfs
		elsif (defined $idTranscr and defined $exStart and not defined $uorfStart) {
			my $transcrEnd = $transcrStart + $transcrLen;

			# Exon positions must be relative to transcription start
			my $exLen = $exEnd - $exStart;
			$exStart = $exStart - $transcrStart;

			if (exists $res{$idGene}) {
				if (exists $res{$idGene}{$idTranscr}) {
					if (exists $res{$idGene}{$idTranscr}{'data'}) {
						if (not exists $res{$idGene}{$idTranscr}{'data'}{'exStarts'}{$exStart}) {
							$res{$idGene}{$idTranscr}{'data'}{'exons'}{$exStart} = $exLen;
						}
					}
					else {
						$res{$idGene}{$idTranscr}{'data'} = {
							'name'    => $ncbiid,
							'chr'     => $chrom,
							'gsymbol' => $symbol,
							'strand'  => $strand,
							'start'   => $transcrStart,
							'end'     => $transcrEnd,
							'cds'     => [$cdsStart, $cdsEnd],
							'exons'   => {$exStart => $exLen}
						};
					}
				}
				else {
					$res{$idGene}{$idTranscr} = {
						'data' => {
							'name'    => $ncbiid,
							'chr'     => $chrom,
							'gsymbol' => $symbol,
							'strand'  => $strand,
							'start'   => $transcrStart,
							'end'     => $transcrEnd,
							'cds'     => [$cdsStart, $cdsEnd],
							'exons'   => {$exStart => $exLen}
						}
					};
				}
			}
			else {
				$res{$idGene} = {
					$idTranscr => {
						'data' => {
							'name'    => $ncbiid,
							'chr'     => $chrom,
							'gsymbol' => $symbol,
							'strand'  => $strand,
							'start'   => $transcrStart,
							'end'     => $transcrEnd,
							'cds'     => [$cdsStart, $cdsEnd],
							'exons'   => {$exStart => $exLen}
						}
					}
				};
			}
		}

		# A gene with a transcript, uorfs and exons
		else {
			my $transcrEnd = $transcrStart + $transcrLen;

			# Exon positions must be relative to transcription start
			my $exLen = $exEnd - $exStart;
			$exStart = $exStart - $transcrStart;

			if (exists $res{$idGene}) {
				if (exists $res{$idGene}{$idTranscr}) {
					if (exists $res{$idGene}{$idTranscr}{'data'}) {
						if (not exists $res{$idGene}{$idTranscr}{'data'}{'exStarts'}{$exStart}) {
							$res{$idGene}{$idTranscr}{'data'}{'exons'}{$exStart} = $exLen;
						}
					}
					else {
						$res{$idGene}{$idTranscr}{'data'} = {
							'name'    => $ncbiid,
							'chr'     => $chrom,
							'gsymbol' => $symbol,
							'strand'  => $strand,
							'start'   => $transcrStart,
							'end'     => $transcrEnd,
							'cds'     => [$cdsStart, $cdsEnd],
							'exons'   => {$exStart => $exLen}
						};
					}
					if (exists $res{$idGene}{$idTranscr}{'uorf'}) {
						next if (exists $res{$idGene}{$idTranscr}{'uorf'}{$idUORF});

						$res{$idGene}{$idTranscr}{'uorf'}{$idUORF} = [$uorfStart, $uorfEnd, $uframe, $uorfName];
					}
					else {
						$res{$idGene}{$idTranscr}{'uorf'} = {$idUORF => [$uorfStart, $uorfEnd, $uframe, $uorfName]};
					}
				}
				else {
					$res{$idGene}{$idTranscr} = {
						'data' => {
							'name'    => $ncbiid,
							'chr'     => $chrom,
							'gsymbol' => $symbol,
							'strand'  => $strand,
							'start'   => $transcrStart,
							'end'     => $transcrEnd,
							'cds'     => [$cdsStart, $cdsEnd],
							'exons'   => {$exStart => $exLen}
						},
						'uorf' => {$idUORF => [$uorfStart, $uorfEnd, $uframe, $uorfName]}
					};
				}
			}
			else {
				$res{$idGene} = {
					$idTranscr => {
						'data' => {
							'name'    => $ncbiid,
							'chr'     => $chrom,
							'gsymbol' => $symbol,
							'strand'  => $strand,
							'start'   => $transcrStart,
							'end'     => $transcrEnd,
							'cds'     => [$cdsStart, $cdsEnd],
							'exons'   => {$exStart => $exLen}
						},
						'uorf' => {$idUORF => [$uorfStart, $uorfEnd, $uframe, $uorfName]}
					}
				};
			}
		}

		# Check which query ids could not be found
		if (exists $ids{$idGene}) {
			delete $ids{$idGene};
		}
	}

	die "ERROR: Id(s) " . join(", ", keys(%ids)) . " are not present in the gene table of the database" if (keys(%ids));

	return \%res, $globalTax, $globalAssembly;
}

# Write HTML for custom track details page
sub writeHtml {

	# The URL of this script with query params
	my $url = $_[0];

	my $dloadName = $_[1];
	my $assembly = $_[2];

	my $html = '<P><H2>Export sequences</H2><P>The following links will take you to the UCSC Table Browser<br>'
		. 'where you can either save the sequences to your PC or send<br>'
		. 'them to Galaxy. All form fields will already be filled out. You<br>'
		. 'just need to select the correct track from the "track" menu, the<br>'
		. 'correct table from the "table" menu and click on get output.<P>'
		. '<a href="BASEURL.FOOBAR&hgta_outFileName=uORFdb_trackdwnl.fa">Download FASTA sequences</a><P>'
		. '<a href="BASEURL.FOOBAR&boolshad.sendToGalaxy=1">Export sequences to Galaxy</a><P><P>';

	my $baseUrl = 
		'http://genome.ucsc.edu/cgi-bin/hgTables?ignoreCookie=1&hgta_group=user&hgta_outputType=sequence' .
		'&db=' .$assembly .
		'&hgt.customText=' . $url;

	# Put table browser URLs into the HTML
	$html =~ s/BASEURL\.FOOBAR/$baseUrl/g;

	return $html;
}

# Create tracks for transcripts and their uORFs
sub drawTracks {
	my %ids          = %{$_[0]};
	my $colorTranscr = $_[1];
	my $colorUORF    = $_[2];

	# What should be highlighted? Transcript, uORF?
	# Which id and in which color?
	my $view           = $_[3];
	my $highlightId    = $_[4];
	my $highlightColor = $_[5];

	my $dbUrl = $_[6];

	# Draw uORFs with introns.
	my $uorfIntrons = $_[7];
	
	my $assembly = $_[8];

	# Add this relative offset to display region to make the output nicer
	my $offsetView = 0.0016;

	# Display this region in the Genome Browser
	my $minView   = -1;
	my $maxView   = -1;
	my $chromView = "";

	my $params    = "browser hide all\n";
	my $trackHead = "";
	my $trackData = "";

	my $wasHighlighted = 0;

	foreach my $idGene (keys(%ids)) {
		foreach my $idTranscr (keys(%{$ids{$idGene}})) {
			$params .= $trackHead . $trackData;
			$trackHead = "";
			$trackData = "";

			my $data = $ids{$idGene}->{$idTranscr}->{'data'};

			my $gsymbol  = $data->{'gsymbol'};
			my $transcrN = $data->{'name'};
			my $chrom    = $data->{'chr'};

			# Genome Browser wants chromosome names that start with "chr"
			$chrom = "chr" . $chrom if ($chrom !~ m/^chr/);

			my $start    = $data->{'start'};
			my $end      = $data->{'end'};
			my $strand   = $data->{'strand'};
			my $cdsStart = $data->{'cds'}->[0];
			my $cdsEnd   = $data->{'cds'}->[1];

			# Sort exon starts and lenghts by exon start position
			my @exStarts = ();
			my @exLens   = ();

			if (exists $ids{$idGene}->{$idTranscr}->{'data'}->{'exons'}) {
				@exStarts = sort {$a <=> $b} keys(%{$ids{$idGene}->{$idTranscr}->{'data'}->{'exons'}});
				@exLens   = map  {$ids{$idGene}->{$idTranscr}->{'data'}->{'exons'}->{$_}} @exStarts;
			}

			my $exCount = @exStarts;
			$exCount = "" if ($exCount == 0);

			# Get global region to display in the Genome Browser
			die "ERROR: Custom track only accepts one chromosome" if ($chromView and $chrom ne $chromView);
			$chromView = $chrom;

			$minView = $start if ($minView == -1 or $minView > $start);
			$maxView = $end   if ($maxView < $end);

			# Print track headers once per gene with custom html details page
			my $html = writeHtml($dbUrl, $transcrN, $assembly);

			$trackHead .=
qq(track name=$transcrN	description="Transcript $transcrN of gene $gsymbol and its uORF(s)"	visibility=1	type=bedDetail);
			$trackData .=
			  "\titemRgb=On\n#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd"
			  . "\titemRgb\tblockCount\tblockSizes\tblockStarts\tid\thtml\n";

			# Use special color to highlight special transcripts and use "full" instead of "dense" for that track
			my $itemRGBtranscr = $colorTranscr;
			if ($highlightId and $view) {
				if ($view eq "transcript" and $highlightId == $idTranscr) {
					$trackHead =~ s/1$/2/;
					$itemRGBtranscr = $highlightColor;
					$wasHighlighted = 1;
				}
			}

			$trackData .=
			  "$chrom\t$start\t$end\t$transcrN\t0\t$strand\t$cdsStart\t$cdsEnd\t$itemRGBtranscr\t"
			  . $exCount . "\t"
			  . join(",", @exLens) . "\t"
			  . join(",", @exStarts) . "\t"
			  . "$transcrN\t$html\n";

			if (exists $ids{$idGene}->{$idTranscr}->{'uorf'}) {
				foreach my $idUORF (sort {$a <=> $b} keys(%{$ids{$idGene}->{$idTranscr}->{'uorf'}})) {
					my $itemRGBuorf = $colorUORF;
					if ($highlightId and $view) {
						if ($view eq "uorf" and $highlightId == $idUORF) {
							$trackHead =~ s/1$/2/;
							$itemRGBuorf    = $highlightColor;
							$wasHighlighted = 1;
						}
					}
					my $uorf = $ids{$idGene}->{$idTranscr}->{'uorf'}->{$idUORF};

					my $uorfStart = $uorf->[0];
					my $uorfEnd   = $uorf->[1];
					my $uorfFrame = $uorf->[2];
					my $uorfName  = $uorf->[3];

					# Get exon overlaps in uORFs, if desired.
					my @exUorfStarts = ();
					my @exUorfLens   = ();

					if ($uorfIntrons) {

						# Check, if exon start occur in uORF
						for (my $i = 0 ; $i <= $#exStarts ; $i++) {

							# Coordinates relative to transcr. start
							my $exStart = $exStarts[$i];
							my $exLen   = $exLens[$i];

							# Absolute genome coordinates
							my $absExStart = $exStart + $start;
							my $absExEnd   = $exLen + $absExStart;

							# Coordinates relative to uORF start
							$exStart = $absExStart - $uorfStart;

							# Exon starts before uORF
							if ($absExStart < $uorfStart) {

								# and is terminated before uORF
								if ($absExEnd < $uorfStart) {
									next;
								}

								# and overlaps with uORF
								else {
									push(@exUorfStarts, 0);
									$exStart = 0;
								}
							}

							# Exon starts in uORF
							elsif ($absExStart < $uorfEnd) {

								# I cannot let a track begin with an intron.
								# Thus, the beginning of the track must be the first exon.
								if (not @exUorfStarts) {
									$exStart   = 0;
									$uorfStart = $absExStart;
								}
								push(@exUorfStarts, $exStart);
							}

							# Exon is located downstream, no overlap.
							# I cannot let a track end with an intron.
							# Thus, the end of the track must be the end of
							# the last overlapping exon.
							else {
								$uorfEnd = $uorfStart + $exUorfLens[$#exUorfLens] + $exUorfStarts[$#exUorfStarts];
								last;
							}

							# The uORF may terminate before the exon.
							# So I need to check, if the full exon lenght
							# can be applied.
							my $relUorfEnd = $uorfEnd - $uorfStart;
							$exLen = $absExEnd - ($uorfStart + $exStart);
							if ($relUorfEnd - $exStart < $exLen) {
								push(@exUorfLens, ($relUorfEnd - $exStart));
								last;
							}
							else {
								push(@exUorfLens, $exLen);
							}
						}
					}

					# Mark the uORF in the CDS region with a thicker bar
					my $thickStart = 0;
					my $thickEnd   = 0;

					# The uORF is completely upstream of CDS
					if ($cdsStart > $uorfEnd) {

						# No thick bar
						$thickStart = $uorfStart;
						$thickEnd   = $uorfStart;
					}

					# The uORF overlaps the CDS in parts
					elsif ($cdsStart == $uorfEnd) {
						$thickStart = $uorfEnd;
						$thickEnd   = $uorfEnd;
					}
					else {
						# The uORF starts upstream of the CDS
						if ($cdsStart > $uorfStart) {
							$thickStart = $cdsStart;

							# The uORF ends within the CDS
							if ($cdsEnd >= $uorfEnd) {
								$thickEnd = $uorfEnd;
							}

							# The uORF continues after the CDS
							else {
								$thickEnd = $cdsEnd;
							}
						}

						# The uORF starts with(in) the CDS or after the CDS
						else {
							if ($uorfStart <= $cdsEnd) {
								$thickStart = $uorfStart;

								# The uORF ends with(in) the CDS
								if ($uorfEnd <= $cdsEnd) {
									$thickEnd = $uorfEnd;
								}

								# The uORF ends after the CDS
								else {
									$thickEnd = $cdsEnd;
								}
							}

							# The uORF starts after the CDS
							else {
								# No thick bar
								$thickStart = $uorfStart;
								$thickEnd   = $uorfStart;
							}
						}
					}
					if ($uorfIntrons) {

						# No overlapping exons, no track
						if (@exUorfLens) {
							$trackData .= "$chrom\t$uorfStart\t$uorfEnd\t$uorfName\t0\t$strand\t$thickStart\t$thickEnd\t$itemRGBuorf\t"
							  . scalar(@exUorfLens) . "\t" . join(",", @exUorfLens) . "\t" . join(",", @exUorfStarts)
							  . "\t$uorfName\t$html" . "\n";
						}
					}
					else {
						$trackData .= "$chrom\t$uorfStart\t$uorfEnd\t$uorfName\t0\t$strand\t$thickStart\t$thickEnd\t$itemRGBuorf\t1\t"
						  . ($uorfEnd - $uorfStart) . "\t0\t\t$uorfName\t$html\n";
					}
				}
			}
		}
	}
	die "ERROR: The id_$view ->$highlightId<- which should have been highlighted was not found"
	  if ($highlightId and $wasHighlighted == 0);

	$params .= $trackHead . $trackData;

	# Add one to browser view start pos. Genome Browser interprets this number as 1-based, but it is 0-based.
	# The rest of the custom tracks is interpreted as 0-based.
	$minView += 1;

	# Make display nicer by fitting the margins
	$offsetView = ($maxView - $minView) * $offsetView;
	$minView -= $offsetView;
	$maxView += $offsetView;

	$params = "browser position $chromView:" . int($minView) . "-" . int($maxView) . "\n" . $params;

	return $params;
}

#========================================#
# Main
#========================================#
sub makeBrowser {

	my $dbh;

	my $action      = $_[0];
	my $cgi         = $_[1];
	my $dbUrl       = $_[2];
	my $uorfIntrons = $_[3] || "";

	$dbh = uorfdb::db::connect;

	die "Cannot connect to database" if (!$dbh);

	my $view           = "";
	my $highlightId    = "";
	my %ids            = ();
	my $params         = "";
	my $track          = "";
	my $link           = "http://genome.ucsc.edu/cgi-bin/hgTracks?";
	my $colorTranscr   = '0,0,0';
	my $colorUORF      = '102,178,255';
	my $highlightColor = '255,0,0';

	# Screen width of client can be used to fit the track window of the Genome Browser to the screen
	my $screenWidth = undef;

	if ($cgi->param("width")) {
		$screenWidth = $cgi->param("width");
	}

	foreach my $idType ("id_gene", "id_transcript", "id_uorf", "id_mutation") {

		# View can be detected by type of parameter that I get
		if ($cgi->param($idType)) {
			$view = $idType;
			$view =~ s/^id_//;

			# Only one id is allowed
			%ids = map {$_ => undef} (split(",", $cgi->param($idType)));
			if (keys(%ids) > 1) {
				die "Only 1 id is allowed";
			}

			# Get the id_gene, if necessary
			if ($view ne "gene") {
				$highlightId = (keys(%ids))[0];
				if ($view eq "transcript") {
					my $resultsR = $dbh->selectall_arrayref("SELECT id_gene from transcript where id = $highlightId");

					if ($resultsR->[0]->[0]) {
						%ids = ($resultsR->[0]->[0] => undef);
					}
					else {
						die "Transcript error";
					}
				}
				elsif ($view eq "uorf") {
					my $resultsR = $dbh->selectall_arrayref(
"SELECT t.id_gene from transcript t INNER JOIN uorf u ON u.id_transcript = t.id where u.id = $highlightId"
					);

					if ($resultsR->[0]->[0]) {
						%ids = ($resultsR->[0]->[0] => undef);
					}
					else {
						die "uORF error";
					}
				}
				elsif ($view eq "mutation") {
					my $resultsR = $dbh->selectall_arrayref(
						"SELECT t.id_gene from transcript t INNER JOIN uorf u ON u.id_transcript = t.id "
						  . "INNER JOIN mutation m ON m.id_uorf = u.id where m.id = $highlightId");

					if ($resultsR->[0]->[0]) {
						%ids = ($resultsR->[0]->[0] => undef);
					}
					else {
						die "Mutation error";
					}
				}
			}

			# Get gene, transcript and uORF metadata
			my ($idsR, $taxon, $assembly) = getGTUdata(\%ids, $dbh);

			# Draw the tracks
			if ($action eq "data") {
				$params = drawTracks($idsR, $colorTranscr, $colorUORF, $view, $highlightId, $highlightColor, $dbUrl, $uorfIntrons, $assembly);
				return $params;
			}
			else {
				# View a specific assembly
				if ($assembly) {
					$link .= "db=$assembly&ignoreCookie=1";
				}

				# View the latest assembly for a taxon
				else {
					die "Taxon <> 9606" if ($taxon != 9606);
					$taxon = "human";
					$link .= "org=$taxon&ignoreCookie=1";
				}

				# If we know the client's screen width use it to adapt the Genome Browser's track view
				$link .= "&pix=$screenWidth" if (defined $screenWidth);
				$link .= "&hgt.customText=" . $dbUrl . "$idType=" . $cgi->param($idType) . "%26action=data";
				$link .= "%26width=$screenWidth" if (defined $screenWidth);
				return $link;
			}
		}
	}
}

sub main {
	eval {
		my $cgi    = new CGI;
		my $action = $cgi->param("action");

		# The base URL of uORFdb Genome Browser services
		my $dbUrl = "FOOBAR/genomebrowser.pl%3F";

		# Called by webbrowser with action = link => provides link for GenomeBrowser
		if ($action eq "link") {
			print makeBrowser("link", $cgi, $dbUrl . "");
		}

		# Called by GenomeBrowser with action = data => prints track data
		elsif ($action eq "data") {

			# The URL of this script with query params
			$dbUrl = $cgi->url(-query => 1);
			$dbUrl =~ s/;/%26/g;
			print makeBrowser("data", $cgi, $dbUrl, "uorfIntrons");
		}
		else {
			die "Unknown action ->$action<-";
		}
	};
	if ($@) {
		print "-" . $@;
	}
}

main();
