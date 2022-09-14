/* Record changes to the db */
create table change (
	id serial primary key,
	username varchar(16),
	ip int,
	ts int
);

/* Stores common terms */
create table reference (
	id serial primary key,
	name varchar(48) not Null unique,
	id_change int references change(id)
);

/* Needed to translate tax id -> name */
create table taxonomy (
	id int primary key,
	name varchar(200),
	lineage varchar(512),
	id_change int references change(id)
);

/* Authors of publications, one per row */
create table author (
	id serial primary key,
	name varchar(64) unique,
	id_change int references change(id) 
);

/* Metadata for publications */
create table publication (
	id int primary key,				/* PubMed ID */
	pubtype varchar(8),
	title varchar(256),
	pubdate varchar(16),
	abstract varchar(3072),
	startpage varchar(64),			/* Sometimes string is used for startpage */
	endpage varchar(64),			/* Sometimes string is used for endpage */
	volume varchar(32),
	issue varchar(32),
	doi varchar(64),
	journal varchar(168),
	journal_abbr varchar(48),
	publisher varchar(168),
	country varchar(32),
	id_change int references change(id)
);

create table content (
	id serial primary key,
	gene_name_paper varchar(512),
	number char,					/* Number */
	length char,					/* Length */
	alt_promoters char,				/* Alternative_Promoters */
	alt_splicing char,				/* Alternative_Splicing */
	tissue_spec_uorfs char,			/* Tissue_Specific_uORFs */
	non_aug_uorfs char,				/* non-AUG_uORFs */
	dist_5_cap char,				/* Distance_from_5-CAP */
	dist_uorf_stop_cds char,		/* Distance_From_uORF-Stop_to_CDS */
	cds_overlap char,				/* CDS_Overlap */
	rna_sec_struct char,			/* RNA_Secondary_Structure */
	cds_repression char,			/* CDS_Repression */
	cds_induction char,				/* CDS_Induction */
	start_site_sel char,			/* Start_Site_Selection */
	nonsense_mediated_decay char,	/* Nonsense-Mediated_Decay */
	mrna_destab char,				/* mRNA_Destabilisation */
	ribosome_load char,				/* Ribosome_Load */
	ribosome_pause_stall char,		/* Ribosome_Pausing_/_Stalling */
	ribosome_shunting char,			/* Ribosome_Shunting */
	kozak_cons_seq char,			/* Kozak_Consensus_Sequence */
	transl_status char,				/* Translational_Status */
	term_context char,				/* Termination_(context) */
	uorf_rna_pept_seq char,			/* uORF_RNA/Peptide_Sequence */
	regul_seq_motif char,			/* Regulatory_Sequence_Motif */
	cof_ribosome_interact char,		/* Co-factor_/_Ribosome_Interaction */
	disease_rel_uorfs char,			/* Disease-related_uORFs */
	acq_mutations_snps char,		/* Acquired_Mutations_/SNPs */
	mouse_models char,				/* Mouse_Models */
	ribosome_prof char,				/* Ribosome_Profiling */
	bioinf_array_screens char,		/* Bioinformatics_/_Array_/_Screens */
	proteomics char,				/* Proteomics */
	methods char,					/* Methods */
	reviews char,					/* Reviews */
	id_publication int references publication(id),
	id_change int references change(id)
);

/* Link author and publication */
create table pubauthor (
	id_publication int references publication(id) ON DELETE CASCADE,
	id_author int references author(id),
	idx int,
	id_change int references change(id)
);

create index idx_pubauthor_id_publication on pubauthor(id_publication);
create index idx_pubauthor_id_author on pubauthor(id_author);

/* Link content and taxonomy table */
create table conttax (
    id_content int references content(id),
    id_taxonomy int references taxonomy(id),
    id_change int references change(id),
    primary key (id_content, id_taxonomy)
);

create index idx_conttax_id_content on conttax(id_content);
create index idx_conttax_id_taxonomy on conttax(id_taxonomy);

create table gene (
	id serial primary key,
	assembly varchar(64),
	geneid varchar(32),				/* ENTREZ geneid */
	genbank_id varchar(32),
	symbol varchar(32),
	names varchar(2600),
	aliases varchar(1024),
	chrom varchar(32),
	id_taxonomy int references taxonomy(id),
	id_change int references change(id),
	unique(id_taxonomy, assembly, geneid, chrom)
);

create index idx_g_geneid on gene(geneid);
create index idx_gene_id_taxonomy on gene(id_taxonomy);

/* Link publication and gene table */
create table contgene (
    id_content int references content(id),
    id_gene int references gene(id),
    id_change int references change(id),
    primary key (id_content, id_gene)
);

create index idx_contgene_id_content on contgene(id_content);
create index idx_contgene_id_gene on contgene(id_gene);

create table transcript (
	id serial primary key,
	ncbiid varchar(32),
	pos int,
	strand char,
	length_noint int,
	length int,
	tlslength int,
	cds_start int,
	cds_end int,
	kozakcontext varchar(10),
	id_kozakstrength int references reference(id),
	id_gene int references gene(id) ON DELETE CASCADE,
	id_change int references change(id),
	unique(ncbiid)
);

create index idx_t_ncbiid on transcript(ncbiid);
create index idx_t_id_gene on transcript(id_gene);

create table exon (
	id serial primary key,
	startpos int,
	endpos int,
	id_transcript int references transcript(id) ON DELETE CASCADE,
	id_change int references change(id),
	unique (startpos, endpos, id_transcript)
);

create table uorf (
	id serial primary key,
	startpos int,
	endpos int,
	startcodon char(3),
	stopcodon char(3),
	cdsdist int,
	cdsdist_noint int,
	length int,
	length_noint int,
	fivepdist int,
	kozakcontext varchar(10),
	id_kozakstrength int references reference(id),
	id_type int references reference(id),
	ucscparams varchar(64),
	rframe smallint,
	id_transcript int references transcript(id) ON DELETE CASCADE,
	id_change int references change(id),
	unique(id_transcript, startpos, endpos)
);

create index idx_u_startpos on uorf(startpos);
create index idx_uorf_id_transcript on uorf(id_transcript);

/* uORF sequences*/
create table seq (
    	id_uorf int references uorf(id),
	part int,
	nucseq varchar(1023),
	aminoseq varchar(341),
	id_change int references change(id),
 	primary key (id_uorf, part)
);

create table variant (
	id serial primary key,
	chrom varchar(6),
	pos int,
	length int,												/* Length of longest allele of the variant: ref or alt */
	allele varchar (100),									/* Allele from VCF: ref or alt */
	isref bool,												/* Reference allele? */
	dbsnpid varchar(16),									/* Only for alt, but not mandatory */
	clinvarid int,											/* Only for alt, but not mandatory */
	unique(chrom, pos, allele),
	id_change int references change(id)
);

create table uorfvar (
	id serial primary key,
	id_uorf int references uorf(id) ON DELETE CASCADE,
	id_variant int references variant(id) on DELETE CASCADE,
	alt_cds_dist int,						/* CDS distance of uORF with the alternative allele */
	alt_uorf_length int,					/* Length of uORF with alternative allele */
	id_location int references reference(id), /* In which region of the uORF is the variant located? */
	id_change int references change(id),
	unique (id_uorf, id_variant)
);

create index idx_uorfvar_id_uorf on uorfvar(id_uorf);

create table uorfvaraltseq (
	id serial primary key,
	id_uorfvar int references uorfvar(id) ON DELETE CASCADE,
	id_type int references reference(id),	/* start/stop codon, Kozak or sequence */
	id_effect int references reference(id), /* Kozak strenght or codon change or other description */
	id_change int references change(id),
	unique (id_uorfvar, id_type)
);

/* Alternate sequences */
create table altseq (
	id_uorfvaraltseq int references uorfvaraltseq(id) ON DELETE CASCADE,
	part int,
	nucseq varchar(1023),
	id_change int references change(id),
	unique(id_uorfvaraltseq, part)
);

create table frequency (
	id_variant int references variant(id) ON DELETE CASCADE,
	freq float,									/* Frequency in cancer or in reference db */
	id_reference int references reference(id), /* Cancer type or refdb name */
	id_change int references change(id),
	primary key (id_variant, id_reference)
);

create table patient (
	id serial primary key,
	accession varchar(12),
	id_cancer int references reference(id),
	unique(accession, id_cancer),
	id_change int references change(id)
);

create table patvariant (
	id_patient int references patient(id) ON DELETE CASCADE,
	id_variant int references variant(id) ON DELETE CASCADE,
	id_change int references change(id),
	primary key (id_patient, id_variant)
);


/* VIEWS */

create view v_authors as
	select *,
	(select count(*) from pubauthor pa where pa.id_author = a.id) as pubcount
	from author a;

create view v_contents as
	select c.*,
	(select p.pubdate from publication p where p.id = c.id_publication),
	(select p.title from publication p where p.id = c.id_publication),
	array_to_string(array(select a.name from author a, pubauthor pa, publication p where pa.id_publication = p.id and p.id = c.id_publication and pa.id_author = a.id order by pa.idx asc), ';') as authors,
	array_to_string(array(select t.name from taxonomy t, conttax ct where ct.id_content = c.id and ct.id_taxonomy = t.id), ',') as taxons,
	array_to_string(array(select g.symbol from gene g, contgene cg where c.id = cg.id_content and g.id = cg.id_gene), ',') as symbols
	from content c;

create view v_genes as
	select g.id, g.assembly, g.symbol, g.names, g.aliases, g.chrom, g.id_taxonomy, g.id_change,
	concat(g.geneid, '#', g.genbank_id) as geneacc,
	(select count(*) from transcript t where t.id_gene = g.id) as transcount,
	(select count(*) from contgene ct where ct.id_gene = g.id) as contcount,
	t.name as taxonomy
	from gene g
	left join taxonomy t on g.id_taxonomy = t.id;

create view v_transcripts as
	select t.*,
	(t.pos + t.length) as endpostrans,
	g.symbol,
	g.chrom,
	(select count(*) from uorf u where startcodon = 'ATG' and u.id_transcript = t.id) as atgcount,
	(select count(*) from uorf u where startcodon != 'ATG' and u.id_transcript = t.id) as non_atgcount,
	(select name from reference r where r.id = t.id_kozakstrength) as kozakstrength
	from transcript t, gene g
	where t.id_gene = g.id;

create view v_uorfs as
	select u.*, t.strand,
	concat(t.ncbiid, '_', u.startcodon, '.', (case
	when t.strand = '+' then (select count(*) + 1 from uorf u2 where u2.id_transcript = u.id_transcript and u2.startcodon = u.startcodon and u2.startpos < u.startpos)
	when t.strand = '-' then (select count(*) + 1 from uorf u2 where u2.id_transcript = u.id_transcript and u2.startcodon = u.startcodon and u2.endpos > u.endpos)
	end)) as name,
	(select g.chrom from gene g where g.id = t.id_gene) as chrom,
	(select count(distinct(uv.id_variant)) from uorfvar uv where uv.id_uorf = u.id) as varcount,
	(select left(nucseq, 513) as nucseq from seq s where s.id_uorf = u.id and s.part = 0),
	(select left(aminoseq, 171) as aminoseq from seq s where s.id_uorf = u.id and s.part = 0),
	concat((select id_taxonomy from gene g where id = t.id_gene), ',', u.id) as vardbsnp,
	concat((select id_taxonomy from gene g where id = t.id_gene), ',', u.id) as varclinvar,
	(select name from reference r where r.id = u.id_kozakstrength) as kozakstrength,
	(select name from reference r where r.id = u.id_type) as type
	from uorf u, transcript t
	where u.id_transcript = t.id;
	
create view v_uorfs_fullseq as
	select u.*, t.strand,
	concat(t.ncbiid, '_', u.startcodon, '.', (case
	when t.strand = '+' then (select count(*) + 1 from uorf u2 where u2.id_transcript = u.id_transcript and u2.startcodon = u.startcodon and u2.startpos < u.startpos)
	when t.strand = '-' then (select count(*) + 1 from uorf u2 where u2.id_transcript = u.id_transcript and u2.startcodon = u.startcodon and u2.endpos > u.endpos)
	end)) as name,
	(select g.chrom from gene g where g.id = t.id_gene) as chrom,
	(select count(distinct(uv.id_variant)) from uorfvar uv where uv.id_uorf = u.id) as varcount,
	array_to_string(array(select nucseq from seq s where s.id_uorf = u.id order by part asc), '') as nucseq,
	array_to_string(array(select aminoseq from seq s where s.id_uorf = u.id order by part asc), '') as aminoseq,
	concat((select id_taxonomy from gene g where id = t.id_gene), ',', u.id) as vardbsnp,
	concat((select id_taxonomy from gene g where id = t.id_gene), ',', u.id) as varclinvar,
	(select name from reference r where r.id = u.id_kozakstrength) as kozakstrength,
	(select name from reference r where r.id = u.id_type) as type
	from uorf u, transcript t
	where u.id_transcript = t.id;

create view v_variants as
    	select uv.id,
    	concat(t.ncbiid, '_', u.startcodon, '.', (case
    	when t.strand = '+' then (select count(*) + 1 from uorf u2 where u2.id_transcript = u.id_transcript and u2.startcodon = u.startcodon and u2.startpos < u.startpos)
    	when t.strand = '-' then (select count(*) + 1 from uorf u2 where u2.id_transcript = u.id_transcript and u2.startcodon = u.startcodon and u2.endpos > u.endpos)
    	end)) as name,
    	concat((select id_taxonomy from gene g where id = t.id_gene), ',', v.id) as vardbsnp,
    	concat((select id_taxonomy from gene g where id = t.id_gene), ',', v.id) as varclinvar,
    	u.startcodon, u.stopcodon, u.kozakcontext, v.pos, v.allele as ref_allele, uv.alt_cds_dist, uv.alt_uorf_length
    	from uorfvar uv, variant v, uorf u, transcript t
    	where t.id = u.id_transcript and v.id = uv.id_variant and u.id = uv.id_uorf and v.isref = 't';

/* Used for variant search in query */
create view v_query_variants as
	select uv.id, v.isref, v2.dbsnpid, v2.clinvarid
	from uorfvar uv, variant v, uorfvar uv2, variant v2
	where v.isref = 't' and v2.isref = 'f'
	and uv.id_variant = v.id and uv2.id_variant = v2.id
	and uv2.id_uorf = uv.id_uorf
	and v2.pos = v.pos;
