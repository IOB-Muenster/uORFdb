#!/usr/bin/env bash


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
# Dump database contents
#=================================================================================================#
#
# DESCRIPTION
#
# This script dumps the database contents for four major fields into separate files.
#	1) Publications and their annotations
#	2) Genes, transcripts, exons, uORFs (incl. full nucleotide and amino acid sequences)
#	3) Variants, annotations in the reference databases, alternate sequences, and mutation
#		effects
#	4) Variants, annotations in the reference databases, and their frequencies
#
#
# USAGE
#
# ./dumpRelations.sh OUTDIR
#
#
# OUTPUT
#
# The output will be generated in the output directory that was specified on the command line.
# For each of the four major fields of interest (see DESCRIPTION), one TSV file is created.
# For every export, there is an example file (example.*) that just contains the
# first 10 records of its associated full file. This is useful to explore the format of the output.
# Each output file is accompanied by a MD5 checksum file (*.md5) that can be used to verify the
# contents of the respective export.
#
#
# CAVEATS
#
# By definition, the cancer frequencies for the ref allele are NA. For the
# web interface we define the frequency of the ref allele as 1 - sum(f(alt1)..f(altn))
# Although the script runs on BSD and macOS, the checksum files will be formated for the use
# with GNU md5sum.
#=================================================================================================#


# Dump data from the database. Fix backslashes.
dumpData () {
	echo -e ${1} > ${2}
	psql -h ${3} -U ${4} -d ${5} -c "COPY (${6}) ${7};" | sed 's/\\\\/\\/g' >> ${2} || { printf "ERROR: Could not dump relation"; exit 1; }
}


# Calculate a check sum file for the provided input file.
getCheckSum () {
	${1} ${2} > "${2}.md5" || { printf "ERROR: Cannot calculate checksum"; exit 1; }
	${1} --quiet -c "${2}.md5" || { printf "ERROR: Calculated checksum does not match"; exit 1; }
}


# Extract first 11 lines from input and calculate checksum
dumpExample () {
	head -n 11 _${1} > "_example.${1}"
	getCheckSum ${MD5SUM} "_example.${1}"
}


#--------------------------------------------------------------------------#
# Parse CLI
#--------------------------------------------------------------------------#
OUTP=$(realpath ${1})

if [ -z "${OUTP}" ]; then
	printf "ERROR: You need to specify an output path\n"
	exit 1
else
	OUTP=$(realpath "${OUTP}")
	if [ ! -d "${OUTP}" ]; then
		printf "ERROR: The output directory ${OUTP} does not exist\n"
		exit 1
	elif [ ! -w "${OUTP}" ]; then
		printf "ERROR: The output directory ${OUTP} is not writable\n"
		exit 1
	fi
fi

cd ${OUTP}


#--------------------------------------------------------------------------#
# Database variables
#--------------------------------------------------------------------------#
HOST="FOO"
USER="BAR"
DB="FOOBAR"


#--------------------------------------------------------------------------#
# SQL queries and output headers
#--------------------------------------------------------------------------#
TSV="TO STDOUT WITH NULL ''"

HEADERCONTTSV="PubMedID\t
Authors\t
Title\tAbstract\tPublicationType\tPublDate\tStartPage\tEndPage\tVolume\tIssue\tDOI\tJournalName\tAbbrJournalName\tPublisher\tCountry\tAlternativePromotors\t
AlternativeSplicing\tTissue-specific_uORFs\tNon-AUGuORFs\tNumber\tLength\tDistanceFrom5'Cap\tDistanceFrom_uORF-STOPtoCDS\tCDSoverlap\tCDSrepression\tCDSinduction\t
StartSiteSelection\tmRNAdestabilization\tNonsense-mediatedDecay\tRNAsecondaryStructure\tRibosomeLoad\tRibosomePausingStalling\tRibosomeShunting\tKozakConsensusSequence\t
TranslationalStatus\tTerminationContext\tuORF_RNApeptideSequence\tRegulatorySequenceMotif\tCo-factorRibosomeInteraction\tDisease-related_uORFs\tAcquiredMutationsSNPs\t
MouseModels\tRibosomeProfiling\tBioinformaticsArraysScreens\tProteomics\tMethods\tReview\t
Taxa\t
GeneIDs\t
GenBankIDs\t
GeneSymbols\t
GeneNameInPaper"
HEADERUORFTSV="Taxon\tAssembly\tChr\t
Symbol\tGeneID\tGenBankID\tSymbolAliases\tGeneNames\t
NCBIID\tTranscrStart\tTranscrEnd\tStrand\tTranscrLength\tTLSlength\tCDSstart\tCDSend\tTranscrKozakContext\t
TranscrKozakStrength\t
ExonStarts\t
ExonEnds\t
uORF_ID\tuORFstart\tuORFend\tuORFstartCodon\tuORFstopCodon\tuORFlength\tuORFCDSdistance\t
uORF5'-capDistance\tuORFkozakContext\tuORFkozakStrength\t
uORFtype\tuORFreadingFrame\tuORFnucleotideSeq\tuORFaminoSeq\t
SharedStartCodon
"
HEADERVARMAPTSV="uORF_ID\tChr\tPosition\tRef\tAlt\tdbSNP_ID\tClinVar_ID\tLocation\tVAR_ID\tAltCDSdistance\tAltuORFlength\tuORFvar_ID"
HEADERVARFREQTSV="VAR_ID\tExAC(ref; alt)\tgnomAD(ref; alt)\tTopMed(ref; alt)\tBRCA_AF(ref; alt)\tCOAD_AF(ref; alt)\tLAML_AF(ref; alt)\tLUAD_AF(ref; alt)\tPRAD_AF(ref; alt)\tSKCM_AF(ref; alt)\tTotal_AF(ref; alt)\t"
HEADERVARSEQTSV="uORFvar_ID\tEffectStartCodon\tAltStartCodon\tEffectStopCodon\tAltStopCodon\tEffectKozakContext\tAltKozakContext\tEffectSequence\tAltSequence"

CMDCONT="select p.id as pubmedid, 
array_to_string(array(select a.name from author a inner join pubauthor pa on a.id = pa.id_author where pa.id_publication = p.id order by pa.idx), '; ') as authors,
p.title, p.abstract, p.pubtype, p.pubdate, p.startpage, p.endpage, p.volume, p.issue, p.doi, p.journal, p.journal_abbr, p.publisher, p.country, c.alt_promoters,
c.alt_splicing, c.tissue_spec_uorfs, c.non_aug_uorfs, c.number, c.length, c.dist_5_cap, c.dist_uorf_stop_cds, c.cds_overlap, c.cds_repression, c.cds_induction,
c.start_site_sel, c.mrna_destab, c.nonsense_mediated_decay, c.rna_sec_struct,c.ribosome_load, c.ribosome_pause_stall, c.ribosome_shunting, c.kozak_cons_seq,
c.transl_status, c.term_context, c.uorf_rna_pept_seq, c.regul_seq_motif, c.cof_ribosome_interact, c.disease_rel_uorfs, c.acq_mutations_snps,
c.mouse_models, c.ribosome_prof, c.bioinf_array_screens, c.proteomics, c.methods, c.reviews,
array_to_string(array(select t.name from taxonomy t inner join conttax ct on ct.id_taxonomy = t.id where ct.id_content = c.id), '; ') as taxa, 
array_to_string(array(select g.geneid from gene g inner join contgene cg on cg.id_gene = g.id where cg.id_content = c.id), '; ') as geneids, 
array_to_string(array(select g.genbank_id from gene g inner join contgene cg on cg.id_gene = g.id where cg.id_content = c.id), '; ') as genbanks, 
array_to_string(array(select g.symbol from gene g inner join contgene cg on cg.id_gene = g.id where cg.id_content = c.id), '; ') as symbols,
c.gene_name_paper from publication p inner join content c on c.id_publication = p.id"

CMDUORF="select (select t.name from taxonomy t where g.id_taxonomy = t.id) as taxon, g.assembly, g.chrom, 
g.symbol, g.geneid, g.genbank_id, replace(g.aliases, '#', '; ') as aliases, replace(g.names, '#', '; ') as names,  
t.ncbiid, t.pos, (t.pos + t.length) as endpos, t.strand, t.length_noint, t.tlslength, t.cds_start, t.cds_end, t.kozakcontext, 
(select r.name from reference r where t.id_kozakstrength = r.id) as kozakstrength, 
string_agg(e.startpos::TEXT, '; ' order by e.startpos asc) as ex_starts, 
string_agg(e.endpos::TEXT, '; ' order by e.endpos asc) as ex_ends, 
vu.name, vu.startpos, vu.endpos, vu.startcodon, vu.stopcodon, vu.length_noint, vu.cdsdist_noint,  
vu.fivepdist, vu.kozakcontext, (select r.name from reference r where vu.id_kozakstrength = r.id) as kozakstrength_uorf, 
(select r.name from reference r where vu.id_type = r.id) as type, vu.rframe, vu.nucseq, vu.aminoseq,
replace(array_to_string(array(select x.ncbiid from transcript x where x.id_gene = t.id_gene and x.id != t.id and
((t.strand = '+' and x.strand = '+' and x.cds_start = vu.startpos and substr(x.kozakcontext,7,3) = vu.startcodon) or 
(t.strand = '+' and x.strand = '-' and x.cds_end = vu.startpos and substr(x.kozakcontext,7,3) = vu.startcodon) or 
(t.strand = '-' and x.strand = '-' and x.cds_end = vu.endpos and substr(x.kozakcontext,7,3) = vu.startcodon) or 
(t.strand = '-' and x.strand = '+' and x.cds_start = vu.endpos and substr(x.kozakcontext,7,3) = vu.startcodon)) order by x.ncbiid), 
','), ',', '; ') as sharedstarts
from gene g 
left outer join transcript t on t.id_gene = g.id left outer join exon e on e.id_transcript = t.id left outer join
v_uorfs_fullseq vu on vu.id_transcript = t.id group by g.id, g.id_taxonomy, t.id, vu.name, vu.startpos, vu.endpos,
vu.strand, vu.rframe, vu.startcodon, vu.stopcodon, vu.cdsdist_noint, vu.length_noint, vu.fivepdist, vu.kozakcontext,
vu.id_kozakstrength, vu.id_type, vu.nucseq, vu.aminoseq"

# Variant metadata to variant ID and uorfvar ID
CMDVARMAP="select vu.name, v.chrom, v.pos,
(select x.allele from variant x where x.chrom = v.chrom and x.pos = v.pos and x.length = v.length and x.isref = 't') as ref,
v.allele as alt, v.dbsnpid, v.clinvarid, (select r.name from reference r where r.id = uv.id_location) as location, v.id, uv.alt_cds_dist, uv.alt_uorf_length,
uv.id from v_uorfs vu inner join uorfvar uv on uv.id_uorf = vu.id inner join variant v on v.id = uv.id_variant where v.isref = 'f'"

# Variant ID to variant frequency
# Coalesce: Use "NA" as default for frequencies from reference dbs. Otherwise problems with multiallelics where one alt has no frequency in reference db.
# Use to_char, to get rid off scientific float representation, i.e. e^*
CMDVARFREQ="select v.id,
concat_ws('; ',
	coalesce(to_char(
		(select z.freq from variant x, frequency z, reference r where isref = 't' and v.pos = x.pos and v.chrom = x.chrom and v.length = x.length and x.id = z.id_variant and
			r.id = z.id_reference and r.name = 'ExAC'), 'FM0.000000'), 'NA'),
	coalesce(to_char(
		(select x.freq from frequency x, reference r where r.id = x.id_reference and x.id_variant = v.id and r.name = 'ExAC'), 'FM0.000000'), 'NA')
) as exacs,
concat_ws('; ',
	coalesce(to_char(
		(select z.freq from variant x, frequency z, reference r where isref = 't' and v.pos = x.pos and v.chrom = x.chrom and v.length = x.length and x.id = z.id_variant and
			r.id = z.id_reference and r.name = 'gnomAD'), 'FM0.000000'), 'NA'),
	coalesce(to_char(
		(select x.freq from frequency x, reference r where r.id = x.id_reference and x.id_variant = v.id and r.name = 'gnomAD'), 'FM0.000000'), 'NA')
) as gnomads,
concat_ws('; ',
	coalesce(to_char(
		(select z.freq from variant x, frequency z, reference r where isref = 't' and v.pos = x.pos and v.chrom = x.chrom and v.length = x.length and x.id = z.id_variant and
			r.id = z.id_reference and r.name = 'TopMed'), 'FM0.000000'), 'NA'),
	coalesce(to_char(
		(select x.freq from frequency x, reference r where r.id = x.id_reference and x.id_variant = v.id and r.name = 'TopMed'), 'FM0.000000'), 'NA')
) as topmeds,
(select concat_ws('; ', 'NA', to_char(x.freq, 'FM0.000000')) from frequency x, reference r where v.id = x.id_variant and x.id_reference = r.id and r.name = 'TCGA-BRCA') as brcas,
(select concat_ws('; ', 'NA', to_char(x.freq, 'FM0.000000')) from frequency x, reference r where v.id = x.id_variant and x.id_reference = r.id and r.name = 'TCGA-COAD') as coads,
(select concat_ws('; ', 'NA', to_char(x.freq, 'FM0.000000')) from frequency x, reference r where v.id = x.id_variant and x.id_reference = r.id and r.name = 'TCGA-LAML') as lamls,
(select concat_ws('; ', 'NA', to_char(x.freq, 'FM0.000000')) from frequency x, reference r where v.id = x.id_variant and x.id_reference = r.id and r.name = 'TCGA-LUAD') as luads,
(select concat_ws('; ', 'NA', to_char(x.freq, 'FM0.000000')) from frequency x, reference r where v.id = x.id_variant and x.id_reference = r.id and r.name = 'TCGA-PRAD') as prads,
(select concat_ws('; ', 'NA', to_char(x.freq, 'FM0.000000')) from frequency x, reference r where v.id = x.id_variant and x.id_reference = r.id and r.name = 'TCGA-SKCM') as skcms,
(select concat_ws('; ', 'NA', to_char(x.freq, 'FM0.000000')) from frequency x, reference r where v.id = x.id_variant and x.id_reference = r.id and r.name = 'Total') as totals
from variant v where v.isref = 'f'"

# uorfvar ID to altseq
CMDVARSEQ="select uv.id,
(select r.name from uorfvaraltseq uva, reference r where r.id = uva.id_effect and uva.id_uorfvar = uv.id and
	(select x.name from reference x where uva.id_type = x.id) = 'start codon') as starteff,
array_to_string(array(select a.nucseq from altseq a inner join uorfvaraltseq uva on uva.id = a.id_uorfvaraltseq where
	(select x.name from reference x where uva.id_type = x.id) = 'start codon' and uva.id_uorfvar = uv.id order by part asc), '') as altstart,
(select r.name from uorfvaraltseq uva, reference r where r.id = uva.id_effect and uva.id_uorfvar = uv.id and
	(select x.name from reference x where uva.id_type = x.id) = 'stop codon') as stopeff,
array_to_string(array(select a.nucseq from altseq a inner join uorfvaraltseq uva on uva.id = a.id_uorfvaraltseq where 
	(select x.name from reference x where uva.id_type = x.id) = 'stop codon' and uva.id_uorfvar = uv.id order by part asc), '') as altstop,
(select r.name from uorfvaraltseq uva, reference r where r.id = uva.id_effect and uva.id_uorfvar = uv.id and
	(select x.name from reference x where uva.id_type = x.id) = 'Kozak context') as kozakeff,
array_to_string(array(select a.nucseq from altseq a inner join uorfvaraltseq uva on uva.id = a.id_uorfvaraltseq where
	(select x.name from reference x where uva.id_type = x.id) = 'Kozak context' and uva.id_uorfvar = uv.id order by part asc), '') as altkozak,
(select r.name from uorfvaraltseq uva, reference r where r.id = uva.id_effect and uva.id_uorfvar = uv.id and
	(select x.name from reference x where uva.id_type = x.id) = 'sequence') as seqeff,
array_to_string(array(select a.nucseq from altseq a inner join uorfvaraltseq uva on uva.id = a.id_uorfvaraltseq where
	(select x.name from reference x where uva.id_type = x.id) = 'sequence' and uva.id_uorfvar = uv.id order by part asc), '') as altseq from uorfvar uv"


#--------------------------------------------------------------------------#
# Use Linux compatibility layer, as options are different on BSD
#--------------------------------------------------------------------------#
CUT="cut"
MD5SUM="md5sum"
OS=`uname -s`

if [ "$OS" == "FreeBSD" ]; then
	CUT="/compat/linux/bin/cut"
	MD5SUM="/compat/linux/bin/md5sum"
fi


#--------------------------------------------------------------------------#
# Dump publication and uORF TSVs
#--------------------------------------------------------------------------#
# Author, content, publication, contgene, conttax, gene symbol, gene ID, GenBank ID
printf "INFO: Dumping publications as TSVs\n"

OUTF="publication_dump_uORFdb.tsv"
dumpData "${HEADERCONTTSV}" "_${OUTF}" "${HOST}" "${USER}" "${DB}" "${CMDCONT}" "${TSV}"
getCheckSum ${MD5SUM} "_${OUTF}"
dumpExample ${OUTF}

# Gene, transcript, exon, uORF (incl. full nucleotide and amino acid sequences)
printf "INFO: Dumping uORFs as TSVs\n"

OUTF="uORF_dump_uORFdb.tsv"
dumpData "${HEADERUORFTSV}" "_${OUTF}" "${HOST}" "${USER}" "${DB}" "${CMDUORF}" "${TSV}"
getCheckSum ${MD5SUM} "_${OUTF}"
dumpExample ${OUTF}


#--------------------------------------------------------------------------#
# Dump TSVs variant maps
#--------------------------------------------------------------------------#
printf "INFO: Dumping variants as TSVs\n"

# Dump mapping and raw frequency file
OUTMAP="tmp.varmap.tsv"
dumpData "${HEADERVARMAPTSV}" "_${OUTMAP}" "${HOST}" "${USER}" "${DB}" "${CMDVARMAP}" "${TSV}"

# Loose id_uorfvar, and alt CDS distance and alt uORF length for mapping to frequency
OUTMAPFREQ="tmp.varmapfreq.tsv"
${CUT} -f10-12 --complement "_${OUTMAP}" > "_${OUTMAPFREQ}"

# Loose id_variant for mapping to altseq
OUTMAPSEQ="tmp.varmapseq.tsv"
${CUT} -f9 --complement "_${OUTMAP}" > "_${OUTMAPSEQ}"

rm "_${OUTMAP}" || { printf "ERROR: Could not remove temporary map"; exit 1; }


#--------------------------------------------------------------------------#
# Dump variant TSVs
#--------------------------------------------------------------------------#
OUTFREQ="tmp.varfreq.tsv"
dumpData "${HEADERVARFREQTSV}" "_${OUTFREQ}" "${HOST}" "${USER}" "${DB}" "${CMDVARFREQ}" "${TSV}"

# Join map and frequency file
OUTF="variant_frequency_dump_uORFdb.tsv"
join -1 9 -2 1 -t $'\t' \
	<(head -n1 "_${OUTMAPFREQ}") \
	<(head -n1 "_${OUTFREQ}") | \
	${CUT} -f1 --complement > "_${OUTF}" || { printf "ERROR: Could not join map and frequency"; exit 1; }
	
join -1 9 -2 1 -t $'\t' \
	<(tail -n+2 "_${OUTMAPFREQ}" | sort -t $'\t' -k 9,9) \
	<(tail -n+2 "_${OUTFREQ}" | sort -t $'\t' -k 1,1) | \
	${CUT} -f1 --complement | sort -k2,3V >> "_${OUTF}" || { printf "ERROR: Could not join map and frequency"; exit 1; }
getCheckSum ${MD5SUM} "_${OUTF}"
dumpExample ${OUTF}
rm "_${OUTMAPFREQ}" "_${OUTFREQ}" || { printf "ERROR: Could not remove temporary map and temporary frequency files"; exit 1; }

OUTSEQ="tmp.varfseq.tsv"
dumpData "${HEADERVARSEQTSV}" "_${OUTSEQ}" "${HOST}" "${USER}" "${DB}" "${CMDVARSEQ}" "${TSV}"

# Join map and sequence file
OUTF="variant_sequence_dump_uORFdb.tsv"
join -1 11 -2 1 -t $'\t' \
	<(head -n1 "_${OUTMAPSEQ}") \
	<(head -n1 "_${OUTSEQ}") | \
	${CUT} -f1 --complement > "_${OUTF}" || { printf "ERROR: Could not join map and frequency"; exit 1; }
	
join -1 11 -2 1 -t $'\t' \
	<(tail -n+2 "_${OUTMAPSEQ}" | sort -t $'\t' -k 11,11) \
	<(tail -n+2 "_${OUTSEQ}" | sort -t $'\t' -k 1,1) | \
	${CUT} -f1 --complement | sort -k2,3V >> "_${OUTF}" || { printf "ERROR: Could not join map and frequency"; exit 1; }
getCheckSum ${MD5SUM} "_${OUTF}"
dumpExample ${OUTF}
rm "_${OUTMAPSEQ}" "_${OUTSEQ}" || { printf "ERROR: Could not remove temporary map and temporary frequency files"; exit 1; }


#--------------------------------------------------------------------------#
# Rename database dumps
#--------------------------------------------------------------------------#
for FILE in $(find ./ -maxdepth 1 -name '_*dump_uORFdb.*sv'); do
	NEWFILE=$(basename ${FILE} | sed 's/^_//')
	mv ${FILE} ${NEWFILE} || { printf "ERROR: Cannot move final dump ${FILE}"; exit 1; }
done

for FILE in $(find ./ -maxdepth 1 -name '_*dump_uORFdb.*sv.md5'); do
	NEWFILE=$(basename ${FILE} | sed 's/^_//')
	sed 's/  _/  /' ${FILE} > ${NEWFILE} || { printf "ERROR: Could not move final checksum ${FILE}"; exit 1; }
	rm ${FILE} || { printf "ERROR: Could not remove temporary checksum ${FILE}"; exit 1; }
done


printf "\nDONE\n"
