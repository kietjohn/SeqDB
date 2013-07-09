import os

from Bio import SeqIO
# from Bio.Blast import NCBIWWW
# from Blast import BlastRecord

import parameters
import tools

class Sequence(object):
    """The Sequence class is used to hold data from one single trace file. Upon
    initiation, the object will try to extract the information required to
    enter the sequence into the record database as designed in the schema file.
    The object can be called as:
        trace_file = Sequence(file_path, plate, instructor, institution[, clone, quality, student, comment])
    (The parameters can be implicit or explicitly stated.)

    This object contains a record() method that returns a list containing the
    validity of the trace file, current file path, and a tuple that agrees with
    the database schema. The prior two attributes can be used to manipulate
    the trace item files. Files that cannot be loaded or does not follow the
    standard naming convention will be listed as non-valid (valid=False), and
    they can be eliminated via other commands.

    """

    def __init__(self, file_path, plate, instructor, institution,
                clone=None, primer=None, student=None,
                comment=None, ID=None):
        self.valid = True
        self.plate = plate
        self.author = student
        self.instructor = instructor
        self.quality = None
        self.institution = institution
        self.path = file_path
        self.comment = comment
        self.db_index = ID

        # Parse the ABI file
        try:
            self.seq = SeqIO.read(self.path, 'abi')
            self.read = self.seq.seq
            tools.generate_single_fasta(self.path)
        except:
            print "The sequence cannot be read!"
            self.valid = False

        # Filter out black listed names
        if self.seq.id in parameters.BLACK_LIST:
            self.valid = False

        # If the author is not supplied, then get sample name and author
        # from ID of trace file.
        if self.author:
            try:
                self.sample_name = self.seq.id.split('-')[0]
            except:
                self.sample_name = self.seq.id
        else:
            try:
                self.sample_name, self.author = self.seq.id.split('-')
            except:
                self.sample_name = self.seq.id
                self.author = None

        # Get start run date from trace file
        for line in str(self.seq).split('\n'):
            if 'run_start' in line:
                self.iden, self.run = line.split('=')
                self.run_date, self.run_time = self.run.split(' ')
        del self.iden, self.run

        # Put the ones with no plate alphabet into exceptions list.
        try:
            int(self.seq.id[0])
            self.valid = False
        except ValueError:
            pass

        if clone:
            self.clone = clone
            self.primer = primer
        else:
            self.clone, self.primer = self.sample_name[:-1], self.sample_name[-1]

        if self.primer != 'F' and self.primer != 'R':
            self.valid = False

        # Quality check, if not already provided
        if not self.quality:
            self.phred_list = self.seq._per_letter_annotations['phred_quality']
            self.good_count = sum([1 for x in self.phred_list
                                    if x > parameters.GOOD_THRESHOLD])

            if self.good_count > parameters.MIN_NUCLEOTIDES:
                self.quality = 1
            else:
                self.quality = 0
        del self.phred_list, self.good_count

        # Assemble sequence storage path
        if self.plate > 100:
            self.plat_string = 'PF' + str(self.plate)
        elif 10 < self.plate < 100:
            self.plate_string = 'PF0' + str(self.plate)
        elif self.plate < 10:
            self.plate_string = 'PF00' + str(self.plate)

        self.file_name = ''.join([self.plate_string, '_',
                                    self.clone, self.primer, '_',
                                    self.run_date])
        self.abi_path = os.path.join('Sequence Files/ABI/', self.plate_string, self.clone[0], self.file_name + '.ab1')
        self.fasta_path = os.path.join('Sequence Files/FASTA/', self.plate_string, self.clone[0], self.file_name + '.fasta')

        # Entry record for this sequence
        if self.valid is True and self.quality != 2:
            self.record = [
                True,
                self.path,
                (self.plate, self.clone, self.primer,
                    self.run_date, self.author,
                    self.instructor, self.institution,
                    self.quality,
                    self.fasta_path, self.abi_path,
                    self.comment)]
        elif self.valid is False:
            self.record = [False, self.path, self.seq.id + '.ab1']

    def __str__(self):
        self.form = '\n'.join(["Plate: " + self.plate_string,
                "Clone: " + ''.join([self.clone, self.primer]),
                "Run date: " + str(self.run_date),
                "Student: " + str(self.author),
                "Instructor: " + self.instructor,
                "Institution: " + self.institution,
                "Sequence: " + str(self.seq.seq),
                "Quality: " + str(self.quality),
                "FASTA file: " + self.fasta_path,
                "Trace file: " + self.abi_path,
                "Comment: " + str(self.comment)
        ])
        return self.form

    def __repr__(self):
        if self.valid is True:
            return "Sequence(file_path= %s, plate= %d, clone= %s, primer= %s, run_date= %s, author= %s, instructor= %s, institution= %s, quality= %d, fasta_path= %s, abi_path= %s, comment= %s)" % (self.record[0], self.record[2])
        if self.valid is False:
            return "Sequence(file_path = %s, id = %s)" % self.record[1:2]

    def __nonzero__(self):
        return self.valid

    def __bool__(self):
        return self.non__zero__()

'''
    def blast(self):
        """Perform BLASTn on the sequence file of this Sequence object,
        uses the E_THRESHOLD from parameters.py."""
        if self.valid is True:
            self.temp_file = 'temp/blast/' + self.file_name + '_blast.xml'
            with open(self.temp_file, 'w') as temp:
                self.blast = NCBIWWW.qblast('blastn', 'nr', self.seq.seq,
                        expect=parameters.E_THRESHOLD, service='dmegablast',
                        alignments=parameters.alignments, megablast=True)
                temp.write(self.blast.read())
            return BlastRecord(self.db_index, raw_blast_result=self.temp_file)
        else:
            print 'This sequence is not suitable for BLAST.'
'''