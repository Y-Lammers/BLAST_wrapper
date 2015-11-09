#!/usr/bin/python

# these are all part of the standard library
import argparse, logging, os, sys, StringIO, csv, time, collections, multiprocessing
from subprocess import call, Popen, PIPE

# The following modules are imported from BioPython
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
from Bio import Entrez


# get the commandline arguments that specify the input fastafile and the output file
parser = argparse.ArgumentParser(description = ('Identify a set of sequences and check if there are CITES species present'))
parser.add_argument('-i', '--input_file', metavar='fasta file', dest='i', type=str, 
			help='input data in FASTA format, the number of sequences is limited to a 100 when running online BLAST searches.', default='', required=True)
parser.add_argument('-o', '--output_file', metavar='output file', dest='o', type=str, 
			help='results file in TSV format. if "-" is provided, output is to STDOUT', default='-')
parser.add_argument('-ba', '--BLAST_algorithm', metavar='algorithm', dest='ba', type=str, 
			help='BLAST algorithm to use (optional, default=blastn)', default='blastn')
parser.add_argument('-task', metavar='task', dest='tk', type=str,
			help='BLAST task to use blastn/megablast', default='megablast')
parser.add_argument('-bd', '--BLAST_database', metavar='database', dest='bd', type=str,
			help = 'BLAST database to use (optional, default=nt)', default = 'nt')
parser.add_argument('-lb', '--local_blast', dest='lb', action='store_true', 
			help = 'blast using a local database (uses the ncbi-blast+ tool, this needs to installed separately)')
parser.add_argument('-tf', '--taxon_file', metavar='taxon file', dest='tf', type=str,
			help = 'Taxon file containing the taxonid\'s + matching scientific names', default = '')
parser.add_argument('-hs', '--hitlist_size', dest='hs', type=int,
			help = 'number of results BLAST will return (optional, default=10), the number of hits is limited to 20 when running online BLAST searches.', default = 10)
parser.add_argument('-mb', '--megablast', dest='mb', action='store_true', 
			help = 'use megablast, can only be used in combination with blastn (optional)')
parser.add_argument('-mi', '--min_identity', dest='mi', type=int, 
			help = 'lowest percentage identity for BLAST results to consider (optional, default=97)', default = 97)
parser.add_argument('-mc', '--min_coverage', dest='mc', type=int, 
			help = 'minimal coverage for BLAST results in number of bases (optional, default=100)', default = 100)
parser.add_argument('-me', '--max_evalue', dest='me', type=float, 
			help = 'threshold E-value for BLAST results (optional, default=0.05)', default = 0.05)
parser.add_argument('-l', '--logging', metavar='log level', dest='l', type=str,
			help = 'set log level to: debug, info, warning (default) or critical', default='warning')
parser.add_argument('-lf', '--log_file', metavar='log file', dest='lf', type=str,
			help = 'specifies a file to log to. if unset, logging is to STDERR', default='')
args = parser.parse_args()


def blast_version ():
	
	# this function checks if the blast version available at the
	# path is version 2.2.28 or higher
	pipe = Popen(['blastn', '-version'], stdout=PIPE, stderr=PIPE)
	version = pipe.communicate()[0]
	try:
		version = int(version.split('\n')[0].split(' ')[1].replace('+','').replace('.',''))
		if version >= 2228:
			return 'pass'
		else:
			return 'low'
	except:
		return 'error'


def parse_taxon ():
	
	# open the taxonid file and return
	# a dictionary with the taxonid's as key
	# and the scientific names as values
	taxon_dic = [{}, {}]
	for taxon in open(args.tf):
		taxon = taxon.strip().split('\t')
		taxon_dic[0][taxon[0]] = taxon[1]
		taxon_dic[1][taxon[1]] = taxon[0]
	return taxon_dic

	
def obtain_tax (code):
	
	# try to obtain the taxon id
	taxon, count = ['unknown ID', 'unknown species'], 0
	while count < 3:
		try:
			# based on the genbank id the taxon id is retrieved from genbank
			Entrez.email = "HTS-barcode-checker@gmail.com"
			handle = Entrez.efetch(db="nucleotide", id= code, rettype="gb",retmode="text")
			record = SeqIO.read(handle, "gb")
	
			# parse through the features and grap the taxon_id
			sub = record.features
			db_xref = [code.split(':')[1] for code in sub[0].qualifiers['db_xref'] if 'taxon' in code][0]
			organism = sub[0].qualifiers['organism'][0]
			return [db_xref, organism]
		except:
			count += 1
	return taxon


def write_results (result, mode):

	# ignore empty strings!
	if result != '':
	
		try:
			csvfile = sys.stdout
			if args.o != '-':
				csvfile = open(args.o, mode)
			writer = csv.writer(csvfile, delimiter = '\t',  quoting=csv.QUOTE_NONE, lineterminator='\n')
			writer.writerow(result)
		except:
			pass
	

def local_blast ():

	# BLAST the sequences against a local database using the ncbi-blast+ sofware package
	logging.debug('Starting local BLAST')
	temp_output = '{0}_temp-blast.tsv'.format(os.path.splitext(args.i)[0])
	cores = multiprocessing.cpu_count()

	# Use the selected BLAST algorithm
	if args.ba == 'blastn':
		BLAST_handle = Popen(['blastn', '-query', args.i, '-out', temp_output, '-db', args.bd, '-task', args.tk,
					'-max_target_seqs', str(args.hs), '-num_threads', str(cores), '-outfmt',
					'6 qseqid sseqid stitle sgi sacc pident length qlen evalue bitscore staxids'],
					stdout=PIPE, stderr=PIPE)
	elif args.ba == 'blastp':
		BLAST_handle = Popen(['blastp', '-query', args.i, '-out', temp_output, '-db', args.bd,
					'-max_target_seqs', str(args.hs), '-num_threads', str(cores), '-outfmt',
					'6 qseqid sseqid stitle sgi sacc pident length qlen evalue bitscore staxids'],
					stdout=PIPE, stderr=PIPE)
	elif args.ba == 'blastx':
		BLAST_handle = Popen(['blastx', '-query', args.i, '-out', temp_output, '-db', args.bd,
					'-max_target_seqs', str(args.hs), '-num_threads', str(cores), '-outfmt',
					'6 qseqid sseqid stitle sgi sacc pident length qlen evalue bitscore staxids'],
					stdout=PIPE, stderr=PIPE)
	elif args.ba == 'tblastn':
		BLAST_handle = Popen(['tblastn', '-query', args.i, '-out', temp_output, '-db', args.bd,
					'-max_target_seqs', str(args.hs), '-num_threads', str(cores), '-outfmt',
					'6 qseqid sseqid stitle sgi sacc pident length qlen evalue bitscore staxids'],
					stdout=PIPE, stderr=PIPE)
	elif args.ba == 'tblastx':
		BLAST_handle = Popen(['tblastx', '-query', args.i, '-out', temp_output, '-db', args.bd,
					'-max_target_seqs', str(args.hs), '-num_threads', str(cores), '-outfmt',
					'6 qseqid sseqid stitle sgi sacc pident length qlen evalue bitscore staxids'],
					stdout=PIPE, stderr=PIPE)
	else:
		logging.critical('Invalid blast algorithm selected')
		sys.exit()

	logging.debug('Waiting for the local blast run to finish')
	logging.debug('BLAST message: {0}'.format(BLAST_handle.communicate()))

	# return the blast output file
	return temp_output
	

def parse_local_blast():

	# obtain the blast results
	temp_output = local_blast()
	
	# obtain the taxon dic if provided
	if args.tf != '':
		taxon_dic = parse_taxon()

	# parse through the output and format the data for filtering
	for line in open(temp_output):
		if '\"' in line: line = line.replace('\"','')
		line = line.strip().split('\t')
		if line[2] != 'N/A' and line[2] != line[1]:
			line[1:3] = [' '.join(line[1:3])]
		else:
			line[1:3] = [line[1]]
		if ';' in line[-1]: line[-1] = line[-1].split(';')[0]
		
		# BOLD species name exception
		if 'BoLD' in args.bd and line[-1] == 'N/A':
			name = line[1].split('|')[1].replace('_', ' ')
			try:
				if name in taxon_dic[1]:
					line[-1] = taxon_dic[1][name]
				else:
					line[-1] = taxon_dic[1][name.split(' ')[0]]
				line += [name]
			except:
				line += [name]
		elif line[-1] == 'N/A':
			line += ['N/A']
		else:	
			# try to add the taxon id
			try:
				line += [taxon_dic[0][line[-1]]]
			except:
				line += ['N/A']

		# send line away for filtering
		filter_hits(line)

	# remove the temp file
	os.remove(temp_output)


def online_blast (seq_list):

	# convert the sequences to a sequence file (stored in the
	# working memory)
	temp = StringIO.StringIO()
	SeqIO.write(seq_list, temp, 'fasta')	
	temp.seek(0,0)

	# BLAST the sequences online against a NCBI database
	logging.debug('BLASTING sequences agaist NCBI')
	result_handle = NCBIWWW.qblast(args.ba, args.bd, temp.read(), megablast=args.mb, hitlist_size=args.hs)

	# return the results handle with the blast results
	return result_handle


def parse_online_blast (seq_list):

	# get the result handle and set the taxon dic
	blast_handle, taxon_dic = online_blast(seq_list), {}

	# use the biopython xml parse module to get the results
	logging.debug('Parsing blast result XML file.')
	blast_list = [item for item in NCBIXML.parse(blast_handle)]

	# walk through the blast results and prepare them for filtering
	for blast_result in blast_list:
		for alignment in blast_result.alignments:
			for hsp in alignment.hsps:
				            		
				# calculate the %identity
				identity = float(hsp.identities/(len(hsp.match)*0.01))


				# grab the genbank number
				gb_num = alignment.title.split('|')[1:4:2]
				gb_num[1] = gb_num[1].split('.')[0]

				# get the taxon id based on the genbank identifier
				if gb_num[0] not in taxon_dic:
					taxon = obtain_tax(gb_num[0])
					taxon_dic[gb_num[0]] = taxon
				else:
					taxon = taxon_dic[gb_num[0]]

				# pull all the results together and sent them to the filter function
				filter_hits([str(blast_result.query), str(alignment.title), str(gb_num[0]), str(gb_num[1]),
						str(identity), str(len(hsp.query)), str(blast_result.query_length),
						str(hsp.expect), str(hsp.bits), taxon[0], taxon[1]])


def filter_hits (blast):
	
	# filter the blast hits, based on the minimum
	# identity, minimum coverage, e-value and the user blacklist
	if float(blast[4]) >= args.mi and int(blast[5]) >= args.mc and float(blast[7]) <= args.me:
		blast[6] = round((float(blast[5])/float(blast[6]))*100, 2)
		results = blast
		del results[3]

		# write the results
		write_results(results, 'a')

def parse_seq_file ():

	# start the local or online blast
	if args.lb == True:
		logging.debug('Starting local blast functions.')
		parse_local_blast()
	else:
		# parse the fasta file
		logging.info('Reading user provided fasta file: %s.' % args.i)
		seq_list = [seq for seq in SeqIO.parse(args.i, 'fasta')]
	
		# check # sequences exceeds the 100 limit for online blasting
		if len(seq_list) > 100:
			logging.critical('To many sequences for online blasting, try using a sequence set'+
					' containing a maximum of a 100 sequences.')
			return

		logging.debug('Starting online blast functions.')
		parse_online_blast(seq_list)


def main ():

	# configure logging
	log_level  = getattr(logging, args.l.upper(), None)
	log_format = '%(funcName)s [%(lineno)d]: %(levelname)s: %(message)s'
	if not isinstance(log_level, int):
		raise ValueError('Invalid log level: %s' % loglevel)
		return	
	if args.lf == '':
		logging.basicConfig(format=log_format, level=log_level)
	else:
		logging.basicConfig(filename=args.lf, filemode='a', format=log_format, level=log_level)

	# Write input commands to log
	logging.debug(' '.join(sys.argv))

	# if the users wants to run a local blast, check if the local blast version is high enough.
	if args.lb == True:
		bv = blast_version()
		if bv == 'pass': logging.info('Local blast version is 2.2.28 or higher')
		elif bv == 'low':
			logging.critical('The local blast version is 2.2.27 or lower, update the blast+ software to a more' +
			' recent version or try using the online blast option.')
			return
		else:
			logging.critical('No local blast version is detected, make sure blast+ is installed correctly and' + 
			' available from the commandline.')
			return

	# check if the maximum number of blast hits exceeds 20 if the user wants
	# to blast the sequences online
	if args.lb == False and args.hs > 20:
		logging.critical('Maximum of online blast hits exceeds the limit of 20. Please lower the number blast hits.')
		return
	
	# Check if input fasta file and output file are provided
	if args.i == '':
		logging.critical('No fasta file or output file provided, see -help for details.')
		return

	# create a blank result file and write the header
	header = [ 
		'Query ID',
		'Hit description',
		'GI',
		'% identity',
		'Hit length',
		'Hit % length',
		'E-value',
		'Bit score',
		'Taxon ID',
		'Species name' 
	]
	write_results(header, 'w')

	# parse through the sequence file, blast all sequences and write the blast hits + CITES info
	logging.info('Processing the sequence file.')
	parse_seq_file()

	logging.critical('Done. Results are written to the: ' + args.o + ' output file')


if __name__ == "__main__":
    main()

