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


#===========================================================================#
#			Create an RSS 2.0 feed from uORFdb
#===========================================================================#
#
# DESCRIPTION:
#
# 	Create an RSS 2.0 feed with the publications from the database.
#
#
# USAGE:
#
#	makeRSS.pl --out FEED.rss [--image IMAGELINK] [--help]
#
#	Create an RSS 2.0 feed FEED.rss containing the publications from the
#	database. Use IMAGELINK as the image for the feed.
#	IMAGELINK, must be a web-accessible link to an image file.
#	
#
# CAVEATS:
#
#	Due to problems with CDATA and XML::RSS, the CDATA tag is lost when
#	parsing and writing an RSS feed. Thus, the feed is not updated, but
#	always created from scratch.
#	XML::RSS needs to be modified to allow dc tags in RSS 2.0: Add 
#	"$self->_out_dc_elements($item);" to the function _out_item_2_0_tags
#	in XML/RSS/Private/Output/V2_0.pm.
#
#============================================================================#


use strict;
use warnings;

use lib '../Database/lib/';
use uorfdb::db;
use XML::RSS::Private::Output::V2_0;

use DBI;
use Getopt::Long;
use POSIX qw(strftime locale_h);
use XML::RSS;


# Translate field names from db to human-readable strings
my %fieldNames = (
			"number" 				=>	"number",
			"length" 				=>	"length",
			"alt_promoters"				=>	"alternative promotors",
			"alt_splicing" 				=>	"alternative splicing",
			"tissue_spec_uorfs" 			=>	"tissue-specific uORFs",
			"non_aug_uorfs" 			=> 	"non-AUG uORFs",
			"dist_5_cap"				=> 	"distance from 5' cap",
			"dist_uorf_stop_cds" 			=> 	"distance from uORF-STOP to CDS",
			"cds_overlap" 				=> 	"CDS overlap",
			"rna_sec_struct"			=> 	"RNA secondary structure",
			"cds_repression"			=> 	"CDS repression",
			"cds_induction"				=> 	"CDS induction",
			"start_site_sel"			=> 	"start site selection",
			"nonsense_mediated_decay" 		=> 	"nonsense-mediated decay",
			"mrna_destab"				=>	"mRNA destabilization",
			"ribosome_load" 			=> 	"ribosome load",
			"ribosome_pause_stall"			=> 	"ribosome pausing/stalling",
			"ribosome_shunting"			=> 	"ribosome shunting",
			"kozak_cons_seq"			=>  "Kozak consensus sequence",
			"transl_status"				=> 	"translational status",
			"term_context"				=> 	"termination context",
			"uorf_rna_pept_seq"			=>	"uORF/RNA peptide sequence",
			"regul_seq_motif"			=>	"regulatory sequence motif",
			"cof_ribosome_interact"			=>	"co-factor/ribosome interaction",
			"disease_rel_uorfs"			=>	"disease-related uORFs",
			"acq_mutations_snps"			=>	"acquired mutations/SNPs",
			"mouse_models"				=>	"mouse models",
			"ribosome_prof"				=> 	"ribosome profiling",
			"bioinf_array_screens"			=> 	"Bioinformatics/arrays/screens",
			"proteomics"				=> 	"proteomics",
			"methods"				=>	"methods",
			"reviews"				=>	"review"
);


# Convert a timestamp (time since epoch) to a date
sub tsToDate {
	my $ts = $_[0];
	
	# Get current time in RFC 822 format
	my $old_locale = setlocale(LC_TIME, "C");
	my $date = strftime("%a, %d %b %Y %H:%M:%S %z", localtime($ts));
	setlocale(LC_TIME, $old_locale);
	
	return $date;	
}


# Get publications from database
sub getPubs {
	my %fieldNames = %{$_[0]};
	
	my %pubmeds = ();
	my %times = ();
	my $dbh = uorfdb::db::connect();
	
	my $sth = "";
	$sth = $dbh->prepare("select vc.*, c.ts, p.abstract, p.doi, p.publisher from v_contents vc inner join change c on ".
		"c.id = vc.id_change inner join publication p on p.id = vc.id_publication;");
	$sth->execute();
	
	# Get PubMed ID and metadata from db.
	# Loop over all matching rows and merge the
	# tags for duplicate (due to join) PubMed IDs.
	my $rowR = $sth->fetchrow_hashref() || undef;	
	while (defined $rowR) {
		my $id = $rowR->{'id_publication'};
		my $ts = $rowR->{'ts'};
		my $title = $rowR->{'title'} || "";
		my $abstract = $rowR->{'abstract'} || "";
		my $doi = $rowR->{'doi'} || "";
		my $publisher =$rowR->{'publisher'} || "";
		my $author = $rowR->{'authors'};
		my @authors = split(";", $author);
		
		$pubmeds{$id} = {} if (not exists $pubmeds{$id});
		
		$pubmeds{$id}->{"author"} = \@authors;
		$pubmeds{$id}->{"title"} = $title;
		$pubmeds{$id}->{"abstract"} = $abstract;
		$pubmeds{$id}->{"doi"} = $doi;
		$pubmeds{$id}->{"publisher"} = $publisher;
		$pubmeds{$id}->{"timestamp"} = $ts;
		
		# Time since epoch to human readable date
		$times{$ts} = tsToDate($ts);
		
		if (not exists $pubmeds{$id}->{"tags"}) {
			$pubmeds{$id}->{"tags"} = {};
		}

		foreach my $field (keys(%fieldNames)) {
				next if (not defined $rowR->{$field});
				
				my $tag = $fieldNames{$field};
				if (not exists $pubmeds{$id}->{"tags"}->{$tag}) {
					$pubmeds{$id}->{"tags"}->{$tag} = 1;
				}
				else {
					$pubmeds{$id}->{"tags"}->{$tag}++;
				}
		}
		$rowR = $sth->fetchrow_hashref() || undef;
	}	
	$dbh->disconnect;
	return \%pubmeds, \%times;
}


#---------------------------------------------------------------------------#
# Read arguments from CLI.
#---------------------------------------------------------------------------#
my $rssF = "";
my $isHelp = 0;
my $image = "";

my $helpMsg = <<EOF;
#===========================================================================#
#			Create an RSS 2.0 feed from uORFdb
#===========================================================================#

 DESCRIPTION:

 	Create an RSS 2.0 feed with the publications from the database.


 USAGE:

	makeRSS.pl --out FEED.rss [--image IMAGELINK] [--help]

	Create an RSS 2.0 feed FEED.rss containing the publications from the
	database. Use IMAGELINK as the image for the feed.
	IMAGELINK, must be a web-accessible link to an image file.
	

 CAVEATS:

	Due to problems with CDATA and XML::RSS, the CDATA tag is lost when
	parsing and writing an RSS feed. Thus, the feed is not updated, but
	always created from scratch.
	XML::RSS needs to be modified to allow dc tags in RSS 2.0: Add 
	"\$self->_out_dc_elements(\$item);" to the function _out_item_2_0_tags
	in XML/RSS/Private/Output/V2_0.pm.

#============================================================================#
EOF

# Get CLI arguments
GetOptions ("out:s"		=> \$rssF,
			"image:s"	=> \$image,
           	'help|?'	=> \$isHelp) or die "ERROR: Unknown option\n$helpMsg";

if ($isHelp) {
	die $helpMsg;
}

my $link = "https://www.bioinformatics.uni-muenster.de/tools/uorfdb/index.pl?";
my $rss = XML::RSS->new(version => '2.0');
my $pubmedsR = "";
my $timesR = "";


#---------------------------------------------------------------------------#
# Start a new feed
#---------------------------------------------------------------------------#
print "INFO: Creating feed ->" . $rssF . "<-\n";
	
$rss->add_module(prefix=>'dc', uri=>'http://purl.org/dc/elements/1.1/');
$rss->channel(
   	title       		=> "uORFdb publications",
   	link        		=> $link,
   	description 		=> "All publications in the uORFdb database.",
   	language 		=> "en",
   	managingEditor		=> 'bioinfo@uni-muenster.de (Wojciech Makalowski)',
   	webMaster		=> 'bioinfo@uni-muenster.de (Norbert Grundmann)',
   	generator		=> "Built using makeRSS.pl by Felix Manske using XML::RSS (https://metacpan.org/pod/XML::RSS)",
   	docs			=> 'https://www.rssboard.org/rss-specification',
   	ttl			=> 1440,
   		
 );

#---------------------------------------------------------------------------#
# Get publications from the database
#---------------------------------------------------------------------------#
if ($image) {	
	$rss->image(
		title  => "uORFdb publications",
		url    => $image,
   		link   => $link
 	);
}

# Get publication from uORFdb.
($pubmedsR, $timesR) = getPubs(\%fieldNames);

if (not keys(%{$pubmedsR})) {
	die "ERROR: No publications in database."
}
else {
	print "INFO: ->" . scalar(keys(%{$pubmedsR})) . "<- publications in the feed.\n";
}

my %pubmeds = %{$pubmedsR};
my %times = %{$timesR};


#---------------------------------------------------------------------------#
# Add publications to the feed
#---------------------------------------------------------------------------#
print "INFO: Adding to feed.\n";

# The biggest timestamp in the database is the publication date of the
# feed. Should avoid reloading by clients, if the feed content has not
# changed since the last time.
my $feedDate = (sort { $b <=> $a } keys(%times))[0];
$feedDate = $times{$feedDate};
$rss->channel(lastBuildDate	=> $feedDate);

foreach my $pubmedID (sort { $a <=> $b } keys(%pubmeds)) {
	# This should introduce a line break and bold print. It is ignored in Zotero.
	# I put the blank before br to get a least a whitespace in Zotero.
	my $descr .= '<![CDATA[<b>Tags from uORFdb:</b> <br>]]>';
	my $tag = "";
	
	# +/- columns from the database
	if (keys(%{$pubmeds{$pubmedID}->{"tags"}})) {
		$tag .= join("; ", sort keys(%{$pubmeds{$pubmedID}->{"tags"}}));
	}
	else {
		$tag = "NONE"
	}
	
	$descr .= $tag;
	$descr .= ' <![CDATA[<br><br><b>Abstract:</b> <br>]]>';
	$descr .= $pubmeds{$pubmedID}->{"abstract"} if ($pubmeds{$pubmedID}->{"abstract"});
	
	my $pubmedLink = 'https://pubmed.ncbi.nlm.nih.gov/' . $pubmedID;
	
	$rss->add_item(
				   # mode is a command for the XML::RSS module
		mode 		=> 'insert',
		title 		=> $pubmeds{$pubmedID}->{"title"},
		link  		=> $pubmedLink,
		guid			=> $pubmedLink,
		description 	=> $descr,
		pubDate 		=> $times{$pubmeds{$pubmedID}->{"timestamp"}},
		categorie 	=> 'Publication',
		dc => {
			title 		=> $pubmeds{$pubmedID}->{"title"},
			creator 	=> $pubmeds{$pubmedID}->{"author"},
			source 		=> $pubmeds{$pubmedID}->{"publisher"},
			identifier 	=> "doi:" . $pubmeds{$pubmedID}->{"doi"},
			subject		=> $tag
		},
	);
}
$rss->save($rssF);

print "\nDONE\n";
