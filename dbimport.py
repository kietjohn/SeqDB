import os, shutil
import datetime
import sqlite3 as sql
import subprocess, sys
from Bio import SeqIO

present = datetime.datetime.now()
time = str(present)
date = str(present.year) + '-' + str(present.month) + '-' + str(present.day) + '_' + str(present.hour) + '-' + str(present.minute)

def checkQuality():
	file_source = os.path.abspath('Input/') + '/'
	file_list = os.listdir(file_source)
	filtered_list, exceptions = abi(file_source)
	good_count = 0
	bad_count = 0
	uncertain_count = 0

	for trace_file in filtered_list:
		read = SeqIO.read(file_source + trace_file, 'abi')
		phred_list = read._per_letter_annotations['phred_quality']
		count_good = [1 for x in phred_list if x > 30]
		count_ok = [1 for x in phred_list if x > 20]
	
		if sum(count_good) > 300:
			shutil.copy2(file_source + trace_file, 'Input/Good/')
			os.remove(file_source + trace_file)
			good_count += 1		
		elif sum(count_ok) > 300:
			shutil.copy2(file_source + trace_file, 'Input/Uncertain/')
			os.remove(file_source + trace_file)
			uncertain_count += 1
		else:
			shutil.copy2(file_source + trace_file, 'Input/Bad/')
			os.remove(file_source + trace_file)
			bad_count += 1
	if exceptions:
		clean_up(file_source, exceptions, 'Exceptions/' + date)
	return '\nDate:' + date + '\nQuality check result:\n\t- Number of good sequences: %d\n\t- Number of bad sequences: %d\n\t- Number of uncertain sequences: %d\n\t- Files in wrong format:%d' % (good_count, bad_count, uncertain_count,len(exceptions))

# Filter a list and return two lists, one is abi files, the other is not abi file.
def abi(source):
	file_list = os.listdir(source)
	abi_files = []
	not_abi = []
	folder = [ 'Good', 'Bad', 'Uncertain']
	file_list = list(set(file_list) - set(folder))

	if not file_list:
		print 'There is no file in this folder.'
		return abi_files, not_abi

	for item in file_list:
		name, ext = os.path.splitext(item)
		if ext == '.abi' or ext == '.ab1':
			abi_files.append(item)
		else:
			not_abi.append(item)
	return abi_files, not_abi

def import_files(source,plate,instructor_name,institution_name,quality,database):
	'''This function takes a list of files, extract information from those files, move them to their organized folders, and save all the operations in a log.'''

	entries = []
	exceptions = []

	black_list = open('black_list.txt','r').read().split('\n')

	# Initiate lists for successful and failed file processing.
	reportSuccess = []
	reportFailure = []

	# First the file list is filtered for files ending with either abi or ab1. The rest are added to the exceptions list.
	filtered_list, exceptions = abi(source)
	if not filtered_list:
		return entries, exceptions
	for entry in exceptions:
		reportFailure.append(entry + ' is not in the abi file format.')

	# Loop through the filtered file list, generate an entry for each file. If there is any failure point, the file will be added to the exceptions list.
	for trace_file in filtered_list:
		try:
			read = SeqIO.read(source + trace_file, 'abi') 
		except:
			exceptions.append(trace_file)
			reportFailure.append(trace_file + ' failed to parse.')
			continue

		# Filter out the files that are not a sample (e.g. water, control).
		if read.id in black_list:
			exceptions.append(trace_file)
			reportFailure.append(trace_file + ' is not a sample.')
			continue

		# Extract clone name, primer name, and sequencing date from the file. Value for plate has to be initiated manually. If the file's primer is not either F or R, then it will be added to the exceptions list.

		sample_info = read.__str__().split('\n')
		for line in sample_info:
			if 'run_finish' in line:
				trunk, info = line.split('=')
				run_date, run_time = info.split(' ')

		try:
			sample_name, author = read.id.split('-')
		except :
			sample_name = read.id
			author = 'none'

		instructor = str(instructor_name)
		institution = str(institution_name)

		#Put the ones with no plate alphabet into exceptions list.
		try:
			check = int(sample_name[0])
			exceptions.append(trace_file)
			reportFailure.append(run_date + ' ' + trace_file + ' does not have plate alphabet.')
			continue
		except ValueError:
			pass

		clone, primer = sample_name[:-1], sample_name[-1]
		if not (primer == 'F' or primer == 'R'):
			exceptions.append(trace_file)
			reportFailure.append(run_date + ' ' + trace_file + ' does not have primer record.')
			continue

		file_name = plate + '_' + sample_name + '_' + run_date

		#Output an ABI and a FASTA file
		abi_dest = 'Sequence Files/ABI/' + plate + '/' + clone[0] + '/'
		abi_path = abi_dest + file_name + '.abi'
		fasta_dest = 'Sequence Files/FASTA/' + plate + '/' + clone[0] + '/'
		fasta_path = fasta_dest + file_name + '.fasta'

		if not os.path.exists(abi_dest): os.makedirs(abi_dest)
		shutil.copyfile(source + trace_file, abi_path)

		if not os.path.exists(fasta_dest): os.makedirs(fasta_dest)
		SeqIO.write(read, fasta_path,'fasta')

		# Add the file information to entries list
		entries.append( (plate[-1], clone, primer, run_date, author, instructor, institution, quality, fasta_path, abi_path) )
		os.remove(source + trace_file)
		reportSuccess.append(trace_file)

	# After looping through the file list, add entries to database, and move all exceptions to the "Exceptions" folder.
	add_to_database(entries,database)
	if len(exceptions) > 0:
		clean_up(source, exceptions, 'Exceptions/' + date)

	print 'Import complete.'
	return reportSuccess, reportFailure

def add_to_database(entries,database):
	'''This function takes a list of tuples to be entered into the database named by the variable "database". '''
	con = sql.connect(database)

# Sequence table:  ID INTEGER PRIMARY KEY AUTOINCREMENT,
#	Plate INT, Clone TEXT, Primer TEXT (F orR), Run_date TEXT,
#	Student TEXT, Instructor TEXT, Institution TEXT, Quality INT (0 or 1),
#	FASTA TEXT NOT NULL, ABI TEXT NOT NULL
#	Pursue INT (0,1, or null), 'Comment' TEXT

	with con:
		cur = con.cursor()
		cur.executemany('INSERT INTO Sequence(Plate, Clone, Primer, Run_date, Student, Instructor, Institution, Quality, FASTA, ABI, Pursue, "Comment") VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, null, null)',entries)

def clean_up(source, file_list,dest):
	'''This function takes in a list of files, and a destination folder. It will copy the files to the destination and then delete the files in the original folder.'''

	if not os.path.exists(dest):	os.makedirs(dest)
	for item in file_list:
		shutil.copy2(source + '/' + item, dest + '/' + item)
		os.remove(source + '/' + item)

def writeToLog(msg):
	try:
		log = open('log.txt','a')
	except:
		log = open('log.txt','w')
	log.write(msg)
	log.close()

def log_process(reportSuccess, reportFailure):
	if not (reportSuccess or reportFailure):
		return ''
	log_message = '\n----------------------------------------------------------------\nDate:' + time + '\n\nFiles that have been successfully imported into the database:\n'
	for success in reportSuccess:
		log_message += '\t- ' + success + '\n'
	log_message += '\nExceptions:\n'
	for failure in reportFailure:
		log_message += '\t- ' + failure + '\n'

	return log_message
