import datetime
from decimal import * 
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import sqlite3 as sql
import os

# The interface of BLAST would be basically a textbox filling the entire screen, and one button for doing the blast

getcontext().prec = 2

def blast(input_file):
	e_default = 1e-50

	sequence = SeqIO.read(input_file,'fasta',generic_dna)
	output_name = 'temp-blast' + str(datetime.datetime.now() )+ '.xml'
	output_file = open(output_name,'w')
	blast_outcome = NCBIWWW.qblast('blastn', 'nr', sequence.seq, expect=e_default, service='dmegablast', alignments=100, megablast='TRUE')
	output_file.write(blast_outcome.read())
	output_file.close()
	return output_name

def parse_blast(file_name):
	File = open(file_name, 'r')
	record = NCBIXML.read(File)
	return record

def search_similarity(ID, blast_object):
	match_list = []
	super_colonizer = ['brassicacearum', 'CP002585', 'F113', 'CP003150', 'Q8r196']
	non_colonizer = ['Pf-5', 'CP000076', 'Pf0-1', 'CP000094', 'SBW25', 'AM181176', 'Q287']

	# If non_colonizer is e_value

	for align in blast_object.alignments:
		for hit in align.hsps:
			genome = align.accession
			organism = align.hit_def
			query_start = hit.query_start
			query_end = hit.query_end
			subject_start = hit.sbjct_start
			subject_end = hit.sbjct_end
			identity = float(Decimal(hit.identities) / Decimal(hit.align_length))
			match_list.append((ID, genome, organism, query_start, query_end, subject_start, subject_end, identity))

	for match in match_list:	# For demo
		print match		# For demo
	
	for match in match_list:
		if match[1] in non_colonizer:
			pursue = 0
			break
		else:
			pursue = 1

	selected = list(match_list[0])
	selected.append(pursue)
	return tuple(selected)

def run_blast(database):
	con = sql.Connection(database)
	with con:
		cur = con.cursor()
		cur.execute("select * from toBlast;")
		while True:
			try:
				ID, trace_file = cur.fetchone()
			except:
				trace_file = None
			if trace_file == None:
				break
			trace_file = 'toblast.fasta'
			#blast_result = blast(trace_file)
			blast_result = 'result.xml'
			result = parse_blast(blast_result)
			entry = search_similarity(ID, result)
			print '\n\nThis is the entry:', entry 		# For result demo
#			cur.execute("insert into Blast(ID, Genome, Organism, Query_from, Query_to, Subject_from, Subject_to, Identity, Pursue) values(?,?,?,?,?,?,?,?,?);", entry)
			# Clean up
#			cur.execute("delete from toBlast where ID=?",int(ID))
#			os.remove(blast_result)

run_blast('Seq.db')
