from Bio import SeqIO
import os

file_list = os.listdir('Exceptions/2012-11-19_22-59')

for trace_file in file_list:
	try:
		read = SeqIO.read('Exceptions/2012-11-19_22-59/' + trace_file,'abi')
		print read.id
	except:
		pass
