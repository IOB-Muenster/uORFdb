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
use warnings;
use DBI;
use Data::Dumper;
use Time::Piece;
use Scalar::Util qw(reftype);
use Getopt::Long;
use Digest::SHA qw(sha256_hex);

# Own mods
use lib 'lib/';
use uorfdb::Import qw(readFile add2table getKey insertHash insertRowwise deleteFromTable getWhoWhen);
use uorfdb::db;
use uorfdb::utils;
use uorfdb::config;

#-----------------------------------------------------#
# Read CLI
#-----------------------------------------------------#

my $usage = <<'EOF';
#======================================================================================================#
# Intialize uORFdb from input files.
#------------------------------------------------------------------------------------------------------#

# DESCRIPTION

	Fill one or multiple relations with one or more input files.
	Foreign keys will be assigned on the fly. Multiple ways to
	manipulate the input data on the fly are supported.


# USAGE

	All input files are tab-delimited.
	
	./importDB.pl
		--main FILE			Output file of the uORF-finder + gene metadata
		--gene FILE			File just containing gene metadata
		--transcript FILE	Transcript file with gene metadata
		--exon FILE			Exon file with gene metadata
		--variant FILE		Output of makeVariations.pl
		--taxonomy FILE		File with taxid, name and lineage
		--publication FILE	Output file from makePublFile.pl
		[--help]


# EXAMPLES

	The treatment of the input data for the single relations is defined
	by the data structure in the tables hash. Each relation has its own subhash
	defining the input file, column names in the database and the indices of the
	respective columns in the input file. You can decide to select the complement
	of the column indices by using the mode 'ignore'. This will load all but the
	specified columns from the file. Otherwise, choose 'include'. Optionally, a
	unique constraint on a set of columns can be enforced prior to insertion by
	specifying the respective column indices in the input file. Exact duplicates
	are dropped regardless of the unique constraint. However, the unique contraint
	will make sure that the input file contains at most exact duplicates.
	Rows can be optionally skipped, if one the fields defined in the NOTEMPTYIDXs
	is empty.
	By default terminal whitespaces are stripped from the fields. In general this
	is a good idea, but not if they encode empty elements. In that case, set 
	"keep_whitespaces" to keep them.
	
		'file' => [
			FILE: FILEPATH: string,
			COLUMNS: COLUMN_NAMES_IN_DB: array,
			COLUMNIDXs: COLUMN_IDXS_IN_FILE: array,
			EXCLUDE/INCLUDE COLUMNIDXS: 'ignore' OR 'include',
			[ UNIQUEIDX: UNIQUE_IDXS: array ]
			[ NOTEMPTYIDX: FIELDS_CANNOT_BE_EMPTY: array ]
			KEEPTERMINALWHITESPACES: keep_whitespaces OR ''
		]
		
	The hash key assignes an internal hash name for the data processing.
	
		'hash' => INTERNAL_HASH: string

	The pre function modifies columns before foreign key lookup, but after reading
	the file. This is useful, if, for example, the input contains a field/ fields
	with multiple values that need to be split before insert or foreign key lookup.
	Split is hardcoded to split by '; '. It only works, if the specified fields (tcol)
	have the same number of entries. The first entry from the first field will be paired
	with the first entry from all other fields that should be split.
	If exec is split_idx, the target column is split and a new column 'idx' is created. When a
	row is split, the items are assigned with an index relative to that row.
	If the INTERNAL_HASH only contains one column the unique function can be selected
	which will unify the column values. If the INTERNAL_HASH only contains one column, empty
	fields can be removed using the "rm_empty" function.
	
		'pre' => [
			{'tcol' => TARGET_COLUMN(S): array,
			'exec' => 'split' | 'split_idx'}
			], | {'exec' => 'unique'} | {'exec' => 'rm_empty'}

	The "foreign" key holds an array of foreign keys that should be looked up. The general
	logic is:
	
		SELECT tcol FROM ttab WHERE qcol = ? and qcol = ?
		
	The hash structure for a foreign key FOREIGN1 looks like this:
	
		FOREIGN1 => {
			'qcol' => QUERY_COLUMN(S): array,
			'tcol' => TARGET_COLUMN: string,
			'ttab' => TARGET_TABLE: string,
			[ 'no_delete'=> undef ],
			'mode' => 'first-strict' | 'first-relaxed' | 'all-strict' | 'all-relaxed'
					(default: 'first-strict'),
			[ 'mode_or' => QUERY_COLUMN_NAME(S): array ],
			[ 'keep_undef' => 1 ]
		}

	Query columns (qcol) need to exist in the internal hash, in order to be queried.
	Thus, generally, they need to be specified in the input file (see above). If
	the qcols are not needed for the database insert, they should be deleted
	(this is the default behaviour of the foreign key lookup). If the columns
	need to be kept for the database insert, 'no_delete'=> undef must be set. The 'mode'
	allows you to specify which results should be returned as the foreign key from
	tcol: first-* or all-*. In first-* mode, only the first foreign key for a query is
	returned. Besides, the query must only yield one foreign key! Otherwise, the program
	will crash. In all-* mode, an array of foreign keys is returned. In *-strict mode,
	the qcols must all have defined values. When 'mode_or' is specified together with
	*-strict, at least one column in the 'mode_or' array must be defined. The mode_or
	function is internally implemented as  * where qcol1 = ? or qcol2 = ? , so multiple
	values can match the condition. All query columns are queried using "or", not just the
	columns in the 'mode_or' array. In *-relaxed mode, the qcols may have undefined/empty
	values. However, if all qcols are undefined, no foreign key lookup is performed and
	the entry is removed from the internal hash (by default). Also, by default, entries
	that did not match a foreign key are removed. However, by specifying 'keep_undef'
	empty foreign keys (internally: value undef) are kept.
	From the structure of the lookup, it is not possible to have two or more foreign keys
	in a hash that have query columns with the same name. To enable this, you may append
	'#' followed by a counter to make them unique (during input reading and foreign key
	lookup). For the lookup, the program will automatically strip the '#' and anything after
	it from the query column. Thus, column names in the database cannot contain a hash.
	
	The mod function allows to modify columns after foreign key lookup, but before the insert.
	COLNAME can be an arbitrary name that is not used in the calculation; only for status
	messages. NEW_COLNAME in the other examples is used as the new column name in the internal
	hash and should thus not be arbitrary.
	The strings in the target column may be converted to upper case
	
		{ COLNAME: str => {
			'tcol' => TARGET_COLUMN(S): array, 
			'exec' =>['upper']}
		}
			
	or a column may be completely deleted from the INTERNAL_HASH.
		
		{ COLNAME: str => {
			'tcol' => TARGET_COLUMN(S): array,
			'exec' =>['delete']}
		}
		
	A new column may be intialized by adding, subtracting or concatenating
	existing columns. If N columns are concatenated, the delimitor (e.g. ',')
	must be repeated N-1 times in exec.
	
		{ NEW_COLNAME: str => {
			'tcol' => TARGET_COLUMNS: array, 
			'exec' => ['+'] | ['-'] | ['.",".']}
		}
		
	A static value may also be used to initialize a new column. The second value in tcol
	simply defines the number of fields in the new column. The first value in tcol provides
	the static value to be used.
		
		{ NEW_COLNAME: str => {
			'tcol' => [STATIC_VALUE, TARGET_COLUMN], 
			'exec' =>['init_static']}
		}

	A new column may also be created by simply copying values from exactly one existing column
		
		{ NEW_COLNAME: str => {
			'tcol' => TARGET_COLUMN: array, 
			'exec' =>[]}
		}

	You may also copy exactly one existing column to a new hash, by passing the target hash
	reference to exec.
		
		{ NEW_COLNAME: str => {
			'tcol' => TARGET_COLUMN: array, 
			'exec' => [ TARGET_HASH_REF ] }l
		}
	
	By setting exec to 'chunk_seqs', you can split the value(s) of one or more target columns
	into chunks. The fields of each individual row in the target columns must all produce the
	same number of parts. The chunksize must be specified for every target column by appending
	it to "chunk_seqs" using "_": E.g.: If you want to chunk two columns into parts of 10 and 20,
	respectively: chunk_seqs_10_20. The modification is performed directly on the target column(s).
	A new column called 'part' holds the index of the chunks.
		
		{ parts: str => {
				'tcol' => TARGET_COLUMN: array, 
				'exec' =>[chunk_seqs_size]}
		}

	By default values from the INTERNAL_HASH are now inserted into the database as 
	completely new entries. The insert will fail, if any constraints of the database
	are violated, e.g. entries already exist. If entries already exist, you may want
	to specify one of two behaviours: Skip ('ignore') or update the complete row with
	the new values. Existing entries must be uniquely identifyably by the target columns
	(conflict_target). Internally, ignore is implemented as DO NOTHING conflict_action
	in the PostgreSQL INSERT. Update is implemented as DO UPDATE SET conflict_action in
	INSERT.
	
		
		'insert' => {
			columns => TARGET COLUMN(S): array, 
			behaviour => 'ignore' OR 'update'
		}
	
	The post function starts after the insert. It can be used to insert an array of data
	into a new column TARGET_COLUMN in another internal TARGET_HASH.
	
		'post' => [{ 
			'insert' => {
				tcol => TARGET_COLUMN: array,
				data => DATA: array,
				thash => TARGET_HASH: hash reference }
		}]


# LIMITATIONS
	
	Column names in the database may not contain '#', as names containing that character are treated
	in a special way by this script (see foreign key lookup in EXAMPLES).
	

# DEPENDENCIES

	CGI
	DBI
	Perl >= v5.26.1
	URI::Escape;
	uorfdb::config
	uorfdb::db
	uorfdb::Import.pm
	uorfdb::utils;

EOF
;

my $help = 0;
my $debug = 0;

# Input files
my $mainF = undef;
my $geneF = undef;
my $transF = undef;
my $exonF = undef;
my $taxonomyF = undef;
my $publF = undef;
my $varF = undef;

my $argC = @ARGV;

GetOptions ("main:s" => \$mainF,
	    "gene:s" => \$geneF,
            "transcript:s"   => \$transF,
            "exon:s"   => \$exonF,
            "taxonomy:s"   => \$taxonomyF,
            "publication:s"   => \$publF,
            "variant:s" => \$varF,
           'help|?' => \$help) or die;
           	
if ($help > 0 or $argC == 0) {
	print $usage;
	exit 0;
}       


# Hashes for each table
my %genes = ();
my %genes_only = ();
my %reference_transcr = ();
my %transcripts = ();
my %exons = ();
my %reference_uorfkozak = ();
my %reference_uorftype = ();
my %uorfs = ();
my %uorfs_seq = ();
my %taxonomies = ();
my %authors = ();
my %publications = ();
my %contents = ();
my %pubauths = ();
my %pubgenes = ();
my %pubtaxa = ();
my %pubcontents = ();
my %variants = ();
my %uorfvars = ();
my %reference_names = ();
my %reference_effects = ();
my %altseqs = ();
my %reference_freqnames = ();
my %frequencies = ();


#-----------------------------------------------------#
# Configure input
#-----------------------------------------------------#
my %tables = (
	'gene' => {
		'file' => [$transF, ['chrom', 'id_taxonomy', 'assembly', 'symbol', 'aliases', 'names', 'geneid'], [1, 11, 12, 14, 15, 16, 17],
			"include", [17]],
		'hash' => \%genes,		
	},
	'gene_only' => {
		'file' => [$geneF, ['id_taxonomy', 'assembly', 'chrom', 'symbol', 'names', 'aliases', 'geneid', 'genbank_id'],
					[0, 1, 2, 3, 4, 5, 6, 7], "include", [6, 7]],
		'hash' => \%genes_only,
		'insert' => {
                        'columns' => ['geneid', 'genbank_id'],
                        'behaviour' => 'ignore'
                },
	},
	'reference_transcr' => {
		'file' => [$transF, ['name'],
			[10], "include", [], [10]],
		'hash' => \%reference_transcr,
		'insert' => {
			'columns' => ['name'], 
			'behaviour' => 'ignore'
		},
	},
	'transcript' => {
		'file' => [$transF, ['ncbiid', 'pos', 'length', 'length_noint', 'cds_start', 'cds_end', 'tlslength', 'strand', 'kozakcontext', 'name', 'geneid'],
			[0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 17], "include", [0,17]],
		'hash' => \%transcripts, 
		'foreign' => [{'id_gene'=> { 'qcol' => ['geneid'], 'tcol' => 'id', 'ttab' => 'gene'}},
					{'id_kozakstrength'=> { 'qcol' => ['name'], 'tcol' => 'id', 'ttab' => 'reference', 'keep_undef' => 1}}]
	},
	'exon' => {
		'file' => [$exonF, ['startpos', 'endpos', 'ncbiid', 'geneid'], [0, 1, 2, 10], "include"],
		'hash' => \%exons,
		'foreign' => [{'id_gene'=> { 'qcol' => ['geneid'], 'tcol' => 'id', 'ttab' => 'gene'}},
		{'id_transcript'=> { 'qcol' => ['id_gene', 'ncbiid'], 'tcol' => 'id', 'ttab' => 'transcript'}}]
	},
	'reference_uorfkozak' => {
		'file' => [$mainF, ['name'],
			[20], "include", [], [20]],
		'hash' => \%reference_uorfkozak,
		'insert' => {
			'columns' => ['name'], 
			'behaviour' => 'ignore'
		},
	},
	'reference_uorftype' => {
		'file' => [$mainF, ['name'],
			[28], "include", [], [28]],
		'hash' => \%reference_uorftype,
		'insert' => {
			'columns' => ['name'], 
			'behaviour' => 'ignore'
		},
	},
	'uorf' => {
		'file' => [$mainF, ['startpos', 'endpos', 'startcodon', 'stopcodon', 'cdsdist', 'cdsdist_noint',
			'length', 'length_noint', 'fivepdist', 'kozakcontext', 'name', 'ucscparams', 'rframe', 'geneid',
			'ncbiid', 'name#2'], [7, 8, 9, 10, 15, 16, 13, 14, 11, 19, 20, 3, 18, 37, 4, 28], "include"],
		'hash' => \%uorfs,
		'foreign' => [{'id_gene'=> { 'qcol' => ['geneid'], 'tcol' => 'id', 'ttab' => 'gene', 'no_delete'=> undef}},
			{'genomeV'=> { 'qcol' => ['geneid'], 'tcol' => 'assembly', 'ttab' => 'gene'}},
			{'id_transcript'=> { 'qcol' => ['id_gene', 'ncbiid'], 'tcol' => 'id', 'ttab' => 'transcript'}},
			{'id_kozakstrength'=> { 'qcol' => ['name'], 'tcol' => 'id', 'ttab' => 'reference', 'keep_undef' => 1}},
			{'id_type'=> { 'qcol' => ['name#2'], 'tcol' => 'id', 'ttab' => 'reference', 'keep_undef' => 1}}
			],
		'mod' => [
			{'kozakcontext' => {'tcol' => ['kozakcontext'], 'exec' =>['upper']}},
			{'ucscparams' => {'tcol' => ['genomeV', 'ucscparams', 'startpos', 'endpos'], 'exec' => ['.",".','.",".', '.",".']}},
			{'genomeV' => {'tcol' => ['genomeV'], 'exec' => ['delete']}},
			]
	},
	'seq' => {
		'file' => [$mainF, ['startpos', 'endpos', 'startcodon', 'cdsdist', 'cdsdist_noint',
			'length', 'length_noint', 'fivepdist', 'kozakcontext', 'rframe', 'geneid',
			'ncbiid', 'nucseq', 'aminoseq'],
			[7, 8, 9, 15, 16, 13, 14, 11, 19, 18, 37, 4, 23, 24], "include"],
		'hash' => \%uorfs_seq,
		'foreign' => [{'id_gene'=> { 'qcol' => ['geneid'], 'tcol' => 'id', 'ttab' => 'gene'}},
			{'id_transcript'=> { 'qcol' => ['id_gene', 'ncbiid'], 'tcol' => 'id', 'ttab' => 'transcript'}},
			{'id_uorf' => {'qcol' => ['startpos', 'endpos', 'startcodon', 'cdsdist', 'cdsdist_noint',
			'length', 'length_noint', 'fivepdist', 'kozakcontext', 'rframe', 'id_transcript'], 'tcol' => 'id', 'ttab' => 'uorf'}}
			],
		'mod' => [{ 'parts' => {'tcol' => ['nucseq', 'aminoseq'], 'exec' => ['chunk_seqs_1023_341']}}]
	}, 
	'taxonomy' => {
		'file' => [$taxonomyF, ['id', 'name', 'lineage'], [0, 1, 2],
			"include", [0]],
		'hash' => \%taxonomies
	},
	'author' => {
		'file' => [$publF, ['name'], [46],
			"include"], #, [46]
		'pre' => [{'tcol' => ['name'], 'exec' => 'split'},
				{'exec' => 'unique'}], # unify author name
		'insert' => { # only insert authors that are not yet present in the db.
			'columns' => ['name'], 
			'behaviour' => 'ignore'
		},
		'hash' => \%authors
	},
	'publication' => {
		'file' => [$publF, ['id', 'pubtype', 'title', 'pubdate', 'abstract', 'startpage', 'endpage', 'volume', 'issue', 'doi',
			'journal', 'journal_abbr', 'publisher', 'country'],
			[40, 45, 47, 42, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57], "include", [40]],
		'insert' => {
			'columns' => ['id'], 
			'behaviour' => 'ignore'
		},
		'hash' => \%publications
	},
	'content' => {
		'file' => [$publF, ['gene_name_paper', 'number', 'length', 'alt_promoters', 'alt_splicing', 'tissue_spec_uorfs',
		'non_aug_uorfs', 'dist_5_cap', 'dist_uorf_stop_cds', 'cds_overlap', 'rna_sec_struct', 'cds_repression', 'cds_induction', 'start_site_sel',
		'nonsense_mediated_decay', 'mrna_destab', 'ribosome_load', 'ribosome_pause_stall', 'ribosome_shunting', 'kozak_cons_seq', 'transl_status',
		'term_context', 'uorf_rna_pept_seq', 'regul_seq_motif', 'cof_ribosome_interact', 'disease_rel_uorfs', 'acq_mutations_snps', 'mouse_models',
		'ribosome_prof', 'bioinf_array_screens', 'proteomics', 'methods', 'reviews', 'id_publication'], 
			[3, 12, 13, 8, 9, 10, 11, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40], "include"],
		'hash' => \%contents,
	},
	'pubauthor' => {
		'file' => [$publF, ['id_publication', 'name'], 
			[40, 46], "include"],
		'pre' => [{'tcol' => ['name'], 'exec' => 'split_idx'}],
		'foreign' => [{'id_author'=> { 'qcol' => ['name'], 'tcol' => 'id', 'ttab' => 'author', 'mode' => 'first-strict'}}],
		'hash' => \%pubauths
	},
	'contgene' => {
		'file' => [$publF, ['gene_name_paper', 'number', 'length', 'alt_promoters', 'alt_splicing', 'tissue_spec_uorfs',
			'non_aug_uorfs', 'dist_5_cap', 'dist_uorf_stop_cds', 'cds_overlap', 'rna_sec_struct', 'cds_repression', 'cds_induction', 'start_site_sel',
			'nonsense_mediated_decay', 'mrna_destab', 'ribosome_load', 'ribosome_pause_stall', 'ribosome_shunting', 'kozak_cons_seq', 'transl_status',
			'term_context', 'uorf_rna_pept_seq', 'regul_seq_motif', 'cof_ribosome_interact', 'disease_rel_uorfs', 'acq_mutations_snps', 'mouse_models',
			'ribosome_prof', 'bioinf_array_screens', 'proteomics', 'methods', 'reviews','genbank_id', 'geneid', 'id_publication'],
			[3, 12, 13, 8, 9, 10, 11, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 6, 5, 40], "include"],
		'pre' => [{'tcol' => ['geneid', 'genbank_id'], 'exec' => 'split'}],
		# Geneid is unique across taxa, says NCBI. Also, GenBank accessions are unique.
		'foreign' => [
			{'id_gene'=> { 'qcol' => ['geneid', 'genbank_id'], 'tcol' => 'id', 'ttab' => 'gene', 'mode' => 'all-strict','mode_or' => ['geneid', 'genbank_id']}},
			{'id_content'=> { 'qcol' => ['gene_name_paper', 'number', 'length', 'alt_promoters', 'alt_splicing', 'tissue_spec_uorfs',
			'non_aug_uorfs', 'dist_5_cap', 'dist_uorf_stop_cds', 'cds_overlap', 'rna_sec_struct', 'cds_repression', 'cds_induction', 'start_site_sel',
			'nonsense_mediated_decay', 'mrna_destab', 'ribosome_load', 'ribosome_pause_stall', 'ribosome_shunting', 'kozak_cons_seq', 'transl_status',
			'term_context', 'uorf_rna_pept_seq', 'regul_seq_motif', 'cof_ribosome_interact', 'disease_rel_uorfs', 'acq_mutations_snps', 'mouse_models',
			'ribosome_prof', 'bioinf_array_screens', 'proteomics', 'methods', 'reviews', 'id_publication'], 'tcol' => 'id', 'ttab' => 'content', 'mode' => 'first-relaxed'}}],
		'hash' => \%pubgenes
	},
	'conttax' => {
		'file' => [$publF, ['id_taxonomy', 'gene_name_paper', 'number', 'length', 'alt_promoters', 'alt_splicing', 'tissue_spec_uorfs',
			'non_aug_uorfs', 'dist_5_cap', 'dist_uorf_stop_cds', 'cds_overlap', 'rna_sec_struct', 'cds_repression', 'cds_induction', 'start_site_sel',
			'nonsense_mediated_decay', 'mrna_destab', 'ribosome_load', 'ribosome_pause_stall', 'ribosome_shunting', 'kozak_cons_seq', 'transl_status',
			'term_context', 'uorf_rna_pept_seq', 'regul_seq_motif', 'cof_ribosome_interact', 'disease_rel_uorfs', 'acq_mutations_snps', 'mouse_models',
			'ribosome_prof', 'bioinf_array_screens', 'proteomics', 'methods', 'reviews', 'id_publication'],
			[2, 3, 12, 13, 8, 9, 10, 11, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40], "include"],
		'pre' => [{'tcol' => ['id_taxonomy'], 'exec' => 'split'}],
		'foreign' => [
			{'id_content'=> { 'qcol' => ['gene_name_paper', 'number', 'length', 'alt_promoters', 'alt_splicing', 'tissue_spec_uorfs',
			'non_aug_uorfs', 'dist_5_cap', 'dist_uorf_stop_cds', 'cds_overlap', 'rna_sec_struct', 'cds_repression', 'cds_induction', 'start_site_sel',
			'nonsense_mediated_decay', 'mrna_destab', 'ribosome_load', 'ribosome_pause_stall', 'ribosome_shunting', 'kozak_cons_seq', 'transl_status',
			'term_context', 'uorf_rna_pept_seq', 'regul_seq_motif', 'cof_ribosome_interact', 'disease_rel_uorfs', 'acq_mutations_snps', 'mouse_models',
			'ribosome_prof', 'bioinf_array_screens', 'proteomics', 'methods', 'reviews', 'id_publication'], 'tcol' => 'id', 'ttab' => 'content', 'mode' => 'all-relaxed'}}],
		'hash' => \%pubtaxa
	},
	'variant' => {
		'file' => [$varF, ['chrom', 'pos', 'allele', 'isref', 'length', 'dbsnpid', 'clinvarid'],
			[16, 0, 1, 2, 3, 7, 8], "include"],
		'hash' => \%variants
	},
	'reference_location' => {
		'file' => [$varF, ['name'],
			[4], "include", [], [4]],
		'hash' => \%reference_names,
		'insert' => {
			'columns' => ['name'], 
			'behaviour' => 'ignore'
		},
	},
	'uorfvar' => { 
		'file' => [$varF, ['chrom', 'pos', 'allele', 'isref', 'length', 'name', 'dbsnpid', 'clinvarid', 'alt_cds_dist', 'alt_uorf_length', 'ncbiid', 'chrom#2', 'startpos', 'endpos', 'nucseq'],
			[16, 0, 1, 2, 3, 4, 7, 8, 5, 6, 17, 16, 18, 19, 20], "include"],
		'hash' => \%uorfvars,
		'foreign' => [{'id_transcript'=> { 'qcol' => ['ncbiid'], 'tcol' => 'id', 'ttab' => 'transcript'}},
			{'id_uorf'=> { 'qcol' => ['id_transcript', 'chrom#2', 'startpos', 'endpos', 'nucseq'], 'tcol' => 'id', 'ttab' => 'v_uorfs_fullseq'}},
			{'id_variant'=> { 'qcol' => ['chrom', 'pos', 'allele', 'isref', 'length', 'dbsnpid', 'clinvarid'], 'tcol' => 'id', 'ttab' => 'variant', 'mode' => 'first-relaxed'}},
			{'id_location'=> { 'qcol' => ['name'], 'tcol' => 'id', 'ttab' => 'reference', 'keep_undef' => 1}}]
	},
	'reference_type' => {
		'file' => [$varF, ['name'],
			[12], "include", [], [12]],
		'hash' => \%reference_names,
		'pre' => [{'tcol' => ['name'], 'exec' => 'split'},
				{'exec' => 'unique'}],
		'insert' => {
			'columns' => ['name'], 
			'behaviour' => 'ignore'
		},
	},
	'reference_effect' => {
		'file' => [$varF, ['name'],
			[13], "include", [], [13], "keep_whitespaces"],
		'hash' => \%reference_effects,
		'pre' => [{'tcol' => ['name'], 'exec' => 'split'},
				{'exec' => 'rm_empty'},
				{'exec' => 'unique'}],
		'insert' => {
			'columns' => ['name'],
			'behaviour' => 'ignore'
		},
	},
	'uorfvaraltseq' => {
		'file' => [$varF, ['chrom', 'ncbiid', 'chrom#2', 'startpos', 'endpos', 'nucseq', 'pos', 'allele', 'isref', 'length', 'dbsnpid', 'clinvarid', 'name', 'name#2'],
			[16, 17, 16, 18, 19, 20, 0, 1, 2, 3, 7, 8, 12, 13], "include", [],[], "keep_whitespaces"],
		'hash' => \%altseqs,
		'pre' => [{'tcol' => ['name', 'name#2'], 'exec' => 'split'}],
		'foreign' => [{'id_type'=> { 'qcol' => ['name'], 'tcol' => 'id', 'ttab' => 'reference'}},
			{'id_effect'=> { 'qcol' => ['name#2'], 'tcol' => 'id', 'ttab' => 'reference', 'keep_undef' => 1}},
			{'id_transcript'=> { 'qcol' => ['ncbiid'], 'tcol' => 'id', 'ttab' => 'transcript'}},
			{'id_uorf'=> { 'qcol' => ['id_transcript', 'chrom#2', 'startpos', 'endpos', 'nucseq'], 'tcol' => 'id', 'ttab' => 'v_uorfs_fullseq'}},
			{'id_variant'=> { 'qcol' => ['chrom', 'pos', 'allele', 'isref', 'length', 'dbsnpid', 'clinvarid'], 'tcol' => 'id', 'ttab' => 'variant', 'mode' => 'first-relaxed'}},
			{'id_uorfvar'=> { 'qcol' => ['id_variant', 'id_uorf'], 'tcol' => 'id', 'ttab' => 'uorfvar'}},
			],
	},
	'altseq' => {
		'file' => [$varF, ['chrom', 'ncbiid', 'chrom#2', 'startpos', 'endpos', 'nucseq#2', 'pos', 'allele', 'isref', 'length', 'dbsnpid', 'clinvarid', 'nucseq', 'name', 'name#2'],
			[16, 17, 16, 18, 19, 20, 0, 1, 2, 3, 7, 8, 11, 12, 13], "include", [], [], "keep_whitespaces"],
		'hash' => \%altseqs,
		'pre' => [{'tcol' => ['nucseq', 'name', 'name#2'], 'exec' => 'split'}],
		'foreign' => [{'id_type'=> { 'qcol' => ['name'], 'tcol' => 'id', 'ttab' => 'reference'}},
			{'id_effect'=> { 'qcol' => ['name#2'], 'tcol' => 'id', 'ttab' => 'reference', 'keep_undef' => 1}},
			{'id_transcript'=> { 'qcol' => ['ncbiid'], 'tcol' => 'id', 'ttab' => 'transcript'}},
			{'id_uorf'=> { 'qcol' => ['id_transcript', 'chrom#2', 'startpos', 'endpos', 'nucseq#2'], 'tcol' => 'id', 'ttab' => 'v_uorfs_fullseq'}},
			{'id_variant'=> { 'qcol' => ['chrom', 'pos', 'allele', 'isref', 'length', 'dbsnpid', 'clinvarid'], 'tcol' => 'id', 'ttab' => 'variant', 'mode' => 'first-relaxed'}},
			{'id_uorfvar'=> { 'qcol' => ['id_variant', 'id_uorf'], 'tcol' => 'id', 'ttab' => 'uorfvar'}},
			{'id_uorfvaraltseq'=> { 'qcol' => ['id_uorfvar', 'id_type', 'id_effect'], 'tcol' => 'id', 'ttab' => 'uorfvaraltseq', 'mode' => 'first-relaxed'}},
			],	
		'mod' => [{ 'parts' => {'tcol' => ['nucseq'], 'exec' => ['chunk_seqs_1023']}}]
	},
	'reference_freqname' => {
		'file' => [$varF, ['name'],
			[10], "include", [], [10]],
		'hash' => \%reference_freqnames,
		'pre' => [{'tcol' => ['name'], 'exec' => 'split'},
				{'exec' => 'unique'}],
		'insert' => {
			'columns' => ['name'], 
			'behaviour' => 'ignore'
		},
	},
	'frequency' => {
		'file' => [$varF, ['chrom', 'pos', 'allele', 'isref', 'length', 'dbsnpid', 'clinvarid', 'name', 'freq'],
			[16, 0, 1, 2, 3, 7, 8, 10, 9], "include"],
		'hash' => \%frequencies,
		'pre' => [{'tcol' => ['name', 'freq'], 'exec' => 'split'}],
		'foreign' => [{'id_reference'=> { 'qcol' => ['name'], 'tcol' => 'id', 'ttab' => 'reference'}},
		{'id_variant'=> { 'qcol' => ['chrom', 'pos', 'allele', 'isref', 'length', 'dbsnpid', 'clinvarid'], 'tcol' => 'id', 'ttab' => 'variant', 'mode' => 'first-relaxed'}}],
		'insert' => {
			'columns' => ['id_variant', 'id_reference'], 
			'behaviour' => 'ignore'
		},
	},
	'reference_cancer' => {
		'file' => [$varF, ['name'],
			[15], "include", [], [15]],
		'hash' => \%reference_freqnames,
		'pre' => [{'tcol' => ['name'], 'exec' => 'split'},
				{'exec' => 'unique'}],
		'insert' => {
			'columns' => ['name'], 
			'behaviour' => 'ignore'
		},
	},
	'patient' => {
		'file' => [$varF, ['accession', 'name'],
			[14, 15], "include"],
		'hash' => \%frequencies,
		'pre' => [{'tcol' => ['accession', 'name'], 'exec' => 'split'}],
		'foreign' => [{'id_cancer'=> { 'qcol' => ['name'], 'tcol' => 'id', 'ttab' => 'reference'}}],
		'insert' => {
			'columns' => ['accession', 'id_cancer'], 
			'behaviour' => 'ignore'
		},
	},
	'patvariant' => {
		'file' => [$varF, ['chrom', 'pos', 'allele', 'isref', 'length', 'dbsnpid', 'clinvarid', 'accession', 'name'],
			[16, 0, 1, 2, 3, 7, 8, 14, 15], "include"],
		'hash' => \%frequencies,
		'pre' => [{'tcol' => ['accession', 'name'], 'exec' => 'split'}],
		'foreign' => [{'id_cancer'=> { 'qcol' => ['name'], 'tcol' => 'id', 'ttab' => 'reference'}},
			{'id_patient'=> { 'qcol' => ['accession', 'id_cancer'], 'tcol' => 'id', 'ttab' => 'patient'}},
			{'id_variant'=> { 'qcol' => ['chrom', 'pos', 'allele', 'isref', 'length', 'dbsnpid', 'clinvarid'], 'tcol' => 'id', 'ttab' => 'variant', 'mode' => 'first-relaxed'}}]
	}
);
#=========================================================#

# Connect to db
my $dbh = uorfdb::db::connect() || undef;
die "ERROR: Undefined or empty database handle" if (not defined $dbh or not $dbh);

print "INFO: Successfully opened ->" . $uorfdb::config::DB . "<-\n";

my $changeh = $dbh->prepare("INSERT INTO change (username, ip, ts) VALUES (?,?,?) RETURNING id;");

# Order of tables must reflect interdependencies by foreign keys
TABLE: foreach my $table ('taxonomy', 'author', 'gene', 'gene_only', 'reference_transcr', 'transcript', 'exon', 'reference_uorfkozak',
	'reference_uorftype', 'uorf', 'seq', 'publication', 'content', 'contgene', 'conttax', 'pubauthor', 'variant', 'reference_location', 'uorfvar',
	'reference_type', 'reference_effect', 'uorfvaraltseq', 'altseq', 'reference_freqname', 'frequency', 'reference_cancer', 'patient', 'patvariant') {
	
	# Switch to indicate, if rowwise insert must be used
	my $insertRow = 0;
	
	# Names of columns with array references
	my %refCols = ();
	
	#-----------------------------------------------------#
	# Get data from input files for all tables
	#-----------------------------------------------------#	
	
	my $file = $tables{$table}{'file'};
	
	my $hashR = $tables{$table}{'hash'};
	
	my $fileN = $file->[0];
	
	if (defined $fileN){
		
		print "INFO: Processing data for table ->$table<-\n";
		
		die "ERROR: File $fileN does not exist or is empty." if (! -s $fileN);
			
		my $colNamesR = $file->[1];
		my $colR = $file->[2];
		my $command = $file->[3];
		my $uniqueIdxR = $file->[4];
		my $notEmptyIdxR = $file->[5];
		my $keepWhitespaces = $file->[6];
					
		# Read from file
		my $dataR = Import::readFile($fileN, $colNamesR, $colR, $command, $uniqueIdxR, $notEmptyIdxR, $keepWhitespaces);
		
		# Push data to table
		$hashR = Import::add2table($hashR, $dataR);
		
		
	}
	next TABLE if (keys(%{$hashR}) == 0);
	
	print "INFO: Processing data for table ->$table<-\n" if (not defined $fileN);
	
	
	#-----------------------------------------------------#
	# Do preprocessing on the data
	#-----------------------------------------------------#
	if (exists $tables{$table}{'pre'}) {
		my @pres = @{$tables{$table}{'pre'}};

		foreach my $preR (@pres) {
			my $exec = $preR->{'exec'};
			
			if ($exec =~ m/^split/) {
				my @tcols = @{$preR->{'tcol'}};
				my %splits = map {$_ => undef} @tcols;
				my @splitList = @tcols;
				my %tmp = ();
				
				# Index field for split columns?
				my $isIdx = 0;
				$isIdx = 1 if ($exec =~ m/idx$/);
				
				# Go through all rows
				foreach (my $i = 0; $i <= $#{$hashR->{$splitList[0]}}; $i++) {
					
					my $elemC = 0;
					my $lastElemC = undef;
					
					# Check, if column can be split (not null) and how many elements I will have after splitting
					foreach my $split (keys(%splits)) {
						# If the value is null, don't split that row for this field -> repeat it like the others
						if (not defined $hashR->{$split}->[$i]) {
							delete $splits{$split};
							next;
						}
						# Amount of elements after splitting
						$elemC = (split("; ", $hashR->{$split}->[$i], -1));
						$lastElemC = $elemC if (not defined $lastElemC);
						die "ERROR: Split columns must have the same number of elements" if (not $elemC == $lastElemC);	
					}
					
					# Split values for one row, if necessary. Otherwise repeat value for row
					# as long as I have split elements.
					foreach my $key (keys(%{$hashR})) {
						if (exists $splits{$key}) {
							my @elems = split("; ", $hashR->{$key}->[$i], -1);
							die "ERROR: Split columns must have the same number of elements" if (@elems != $elemC);
					
							if (exists $tmp{$key}) {
								push(@{$tmp{$key}}, @elems)
							}
							else {
								$tmp{$key} = \@elems
							}
							
							# Add index column, if necessary
							if ($isIdx == 1) {
								if (exists $tmp{"idx"}) {
									push(@{$tmp{"idx"}}, (0..$elemC -1))
								}
								else {
									$tmp{"idx"} = [0..$elemC -1]
								}
							}
						}
						else {
							if (exists $tmp{$key}) {
								push(@{$tmp{$key}}, ($hashR->{$key}->[$i]) x $elemC)
							}
							else {
								$tmp{$key} = [($hashR->{$key}->[$i]) x $elemC]
							}
						}
					}

					# Reinitialize the split columns, in case one field was null
					%splits = map {$_ => undef} @tcols;
				}
				# Overwrite table hash
				%{$hashR} = %tmp;
			}
			# Make values unique. Dangerous with multiple columns (would loose
			# connections of fields). Only allowed, if hash only contains  one column.
			elsif ($exec eq 'unique') {
				if (keys(%{$hashR}) != 1) {
					die "ERROR: Unique only allowed, if table contains one column. Found " . keys(%{$hashR}) . ".\n"
				}
				else {
					my $column = (keys(%{$hashR}))[0];
					my %tmp = map {$_ => undef} @{$hashR->{$column}};
					$hashR->{$column} = [keys(%tmp)];
				}
			}
			# Remove empty fields. Dangerous with multiple columns (would loose
			# connections of fields). Only allowed, if hash only contains one column.
			elsif ($exec eq 'rm_empty') {
				if (keys(%{$hashR}) != 1) {
					die "ERROR: Remove empty only allowed, if table contains one column. Found " . keys(%{$hashR}) . ".\n"
				}
				else {
					my $column = (keys(%{$hashR}))[0];
					my @tmps = ();
					foreach my $element (@{$hashR->{$column}}) {
						if ($element !~ m/^\s*$/) {
							push(@tmps, $element)
						} 
					}
					$hashR->{$column} = \@tmps;
				}
			}
			else {
				die "ERROR: Unknown preprocessing mode $exec";
			}
		}
	}
	
	#-----------------------------------------------------#
	# Get foreign keys
	#-----------------------------------------------------#
	if (exists $tables{$table}{'foreign'}) {
		my @foreigns = @{$tables{$table}{'foreign'}};
		
		foreach my $foreignR (@foreigns) {
			
			my $foreignKey = (keys(%{$foreignR}))[0];
			
			# Get query columns from hashR
			my @queryCols = @{$foreignR->{$foreignKey}->{'qcol'}};
			my $tcol = $foreignR->{$foreignKey}->{'tcol'};
			my $ttable = $foreignR->{$foreignKey}->{'ttab'};
			
			# Keep undefined foreign keys
			my $keepUndef = 0;
			$keepUndef = 1 if (exists $foreignR->{$foreignKey}->{'keep_undef'});
			
			my %hash = %{$hashR};
			
			my $mode = "";
			$mode = $foreignR->{$foreignKey}->{'mode'} if (exists $foreignR->{$foreignKey}->{'mode'});
			my $modeOrR = $foreignR->{$foreignKey}->{'mode_or'} if (exists $foreignR->{$foreignKey}->{'mode_or'}) || undef;
			
			# Mark foreign keys with array references
			$refCols{$foreignKey} = undef if ($mode =~ m/all/);
			
			# Default insert into db must not be used if foreign mode has ever been 'all'
			$insertRow = 1 if ($mode ne 'all');
			
			print "INFO: Getting foreign key ->$foreignKey<- from column ->$tcol<- in table ->$ttable<- using mode ->$mode<-\n";
			
			if (defined $modeOrR) {
				print "INFO: At least one column of ".join(",", @{$modeOrR})." must be defined\n";
			}
			
			my %queries = %hash{@queryCols};
			
			# Remove query columns from table hash, if desired
			delete @hash{@queryCols} if (not exists $foreignR->{$foreignKey}->{'no_delete'});
			
			# Get foreign key using queries
			my ($keysR, $remIdxR) = Import::getKey (\%queries, $tcol, $ttable, $dbh, $mode, $modeOrR, $keepUndef);
			
			# Push foreign key to hashR
			$hash{$foreignKey} = $keysR;
			
			# Optionally delete rows by index from table hash.
			# These rows have a null foreign key.
			if (@{$remIdxR}) {
				print "INFO: Deleting foreign keys with null values\n";
				Import::deleteFromTable(\%hash, $remIdxR);
			}
			
			# No elements remaining
			if (not keys(%hash)) {
				print "WARNING: No data remaining. Skipping...\n";
				next TABLE;
			}
			
			$hashR = \%hash;
		}		
		
	}
	
	#-----------------------------------------------------#
	# Get metadata from insert and push to change table
	#-----------------------------------------------------#
	$changeh->execute($uorfdb::config::USER, uorfdb::utils::toSQL("127.0.0.1", "ip"), time);
	my $id_change = $changeh->fetchall_arrayref->[0]->[0];
	
	my @changes = ($id_change) x @{$hashR->{(keys(%{$hashR}))[0]}};
	$hashR -> {'id_change'} = \@changes;
	
	
	#-----------------------------------------------------#
	# Modify/delete existing keys or create a new one
	#-----------------------------------------------------#
	if (exists $tables{$table}{'mod'}) {
		
		my @mods = @{$tables{$table}{'mod'}};
		
		foreach my $modR (@mods) {
			
			my $modKey = (keys(%{$modR}))[0];
			
			print "INFO: Modifying column ->$modKey<-\n";
			
			# Input columns
			my @tcols = @{$modR->{$modKey}->{'tcol'}};
			
			# Operations to be performed, cp if empty
			my @execs = @{$modR->{$modKey}->{'exec'}};
			
			# Convert values to upper case
			if (@execs == 1 and $execs[0] eq 'upper') {
				foreach my $tcol (@tcols) {
					print "INFO: Converting ->$tcol<- to upper case\n";
					
					foreach my $field (@{$hashR->{$tcol}}) {
						$field = uc($field);
					}
				}
			}
			# Copy to a new hash
			elsif (@execs == 1 and defined reftype ($execs[0])) {
				
				die "ERROR: Illegal reference type for mod command" if (reftype ($execs[0]) ne "HASH");
				
				foreach my $tcol (@tcols) {
					print "INFO: Copying ->$tcol<- to new hash\n";
					$execs[0]->{$tcol} = $hashR->{$tcol};
				}
			}
			# Delete a key
			elsif (@execs == 1 and $execs[0] eq 'delete') {
				
				foreach my $tcol (@tcols) {
					print "INFO: Deleting key ->$tcol<-\n";
					delete $hashR->{$tcol};
				}
			}
			# Initialize a key with a static value
			elsif (@execs == 1 and $execs[0] eq 'init_static') {
				
				print "INFO: Creating key ->$modKey<-\n";
				my @res = ($tcols[0]) x @{$hashR->{$tcols[1]}};
				$hashR -> {$modKey} = \@res;
			}
			# Split a sequence in chunks
			elsif (@execs == 1 and $execs[0] =~ m/^chunk_seqs/) {
				print "INFO: Splitting sequence in chunks\n";
				
				# Users can provide a chunk size with the exec command
				# Chunk sizes and target columns are matched by index
				my %targets = ();
				my @chunkSizes = split("_", $execs[0]);
				splice(@chunkSizes, 0, 2);
				
				die "ERROR: There must be a chunk size for every target column" if (@chunkSizes != @tcols);
				
				for (my $i = 0; $i <= $#tcols; $i++) {
					$targets{$tcols[$i]} = $chunkSizes[$i]
				}
					
				# Temporary hash to store the new data
				my %tmp = map { $_ => [] } keys(%{$hashR});
				
				die "ERROR: Special column ->part<- already exists" if (exists $tmp{'part'});
				$tmp{'part'} = [];
				
				# Keys unaffected by chunking
				my @staticKeys = ();
				foreach my $key (keys(%{$hashR})) {
					push(@staticKeys, $key) if (not exists $targets{$key});	
				}
				
				# Arbitrarily select the first target column and use it as the basis
				# for chunking
				my $tcol = $tcols[0];
				for (my $i = 0; $i <=$#{$hashR->{$tcol}}; $i ++) {
					my $j = 0;
					my $seq = $hashR->{$tcol}->[$i];
					my $chunkSize = $targets{$tcol};
					
					# Chunk the first target column, but keep other columns static
					while($seq) {
						foreach my $static (@staticKeys) {
							push(@{$tmp{$static}}, @{$hashR->{$static}}[$i]);
						}						
						push(@{$tmp{$tcol}}, substr($seq, 0, $chunkSize, ""));
						
						# Process the other target columns and match them with the first
						foreach my $otherTcol (@tcols[1..$#tcols]) {
							my $otherSeq = $hashR->{$otherTcol}->[$i];
							my $otherChunkSize = $targets{$otherTcol};
							
							# Add chunk to temporary hash.
							# Save rest of string to input hash.
							push(@{$tmp{$otherTcol}}, substr($otherSeq, 0, $otherChunkSize, ""));
							$hashR->{$otherTcol}->[$i] = $otherSeq;
							
							if (not $seq and $otherSeq) {
								die "ERROR: Target columns must produce equal numbers of chunks"
							}
							elsif ($seq and not $otherSeq) {
								die "ERROR: Target columns must produce equal numbers of chunks"
							}
						}
						
						push (@{$tmp{'part'}}, $j);
						$j++;
					}
				}				
				# Replace table hash with temporary hash
				%{$hashR} = %tmp;
				%tmp = ();
			}
			# Concatenate columns or add/subtract. "Off-label use": Copy
			# column to new column in the current hash, if no execs.
			elsif (@tcols == @execs+1) {
							
				# Order matters for interlace!
				my $mod = '';
	
				for (my $i=0; $i <= $#tcols; $i++) {
					if ($i == 0) {
						$mod .= '$hashR->{'.$tcols[$i].'}->[$_]'
					}
					else {
						$mod .= $execs[($i-1)].'$hashR->{'.$tcols[$i].'}->[$_]';
					}
		 
				}
					
				# Perform the modification
				my @res = map {eval($mod)} 0..$#{$hashR->{$tcols[0]}};
				$hashR -> {$modKey}= \@res;
			}
			else {
				die "ERROR: Cannot modify column. Illegal number of execute statements."
			}
		}
	}
	
	if (exists $tables{$table}{'post'}) {
		$insertRow = 1
	}
	
	my $tmpTime = Time::Seconds->new(time() - $^T)->pretty;
	print "INFO: Finished pre-processing in $tmpTime\n";
	
	
	#-----------------------------------------------------#
	# Insert into database
	#-----------------------------------------------------#
	# Find out, how to handle duplicates between db and table hash
	my $insertR = undef;
	my $conflColsR = [];
	my $conflBehav = undef;
	
	if (exists $tables{$table}{'insert'}) {
		$insertR = $tables{$table}{'insert'};
		if (exists $insertR->{'columns'} and $insertR->{'behaviour'}) {
			$conflColsR = $insertR->{'columns'};
			$conflBehav = $insertR->{'behaviour'};
		}
	}
	
	# Sometime multiple virtual tables are used to fill the same db table
	my $dbTable = $table;
	$dbTable = (split("_", $table))[0] if ($table =~ m/_/);
	
	if ($insertRow == 0) {
		# Default insert
		Import::insertHash($hashR, $dbTable, $dbh, $conflColsR, $conflBehav);
	}
	else {
		# Insert which can handle arrays of arrays as hash values
		Import::insertRowwise($hashR, $dbTable, $dbh, \%refCols, $conflColsR, $conflBehav);
	}
	
	$dbh->commit;
	
	
	#-----------------------------------------------------#
	# Post-processing
	#-----------------------------------------------------#
	if (exists $tables{$table}{'post'}) {
		
		my @posts = @{$tables{$table}{'post'}};
		
		foreach my $postR (@posts) {
			
			my $exec = (keys(%{$postR}))[0];
			
			die "ERROR: Undefined post-processing $exec" if ($exec ne "insert");
			
			my $tcol = $postR->{$exec}->{'tcol'};
			my @data = @{$postR->{$exec}->{'data'}};
			my $hashR = $postR->{$exec}->{'thash'};
		
			$hashR->{$tcol} = \@data;
		}
	}
	
	# Empty the current table hash
	# I cannot do %{$hashR} = (), because $hashR is pointing to %hash
	# which is a subset of my table hash. It does not point to the table
	# hash anymore.
	%{$tables{$table}{'hash'}} = ();
}

# Disconnect
$dbh->disconnect or die $DBI::errstr;
print "INFO: Closed connection to " . $uorfdb::config::DB . "\n";


# $^T is the start time variable
my $execTime = Time::Seconds->new(time() - $^T)->pretty;
print "INFO: Finished in $execTime\n";
print "\nDONE\n";
