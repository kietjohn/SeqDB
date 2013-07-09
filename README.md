SequenceRecord
==============

**Using the database**

1. Copy/Move the sequence files into the Input folder.
2. Click the "Check Quality" button to separate out the good and bad quality ones. Then click "View Uncertains" to see the files whose quality is inbetween good and bad. Judge if the sequence is *good enough* for giving meaningful BLASTn result. If yes, move it to the Good folder; if not, Bad folder.
3. Input the information about the sequence files to the fields in Step 2, and check that the information are correct.
4. Click the "Import Files" button. This will import the sequence files into the SQL database, create FASTA version of each sequence, and organize all the files into folders according to their Plate/Clone Alphabet. This will also produce a log file, which is equivalent to the transaction log on the interface.
5. If there is any exception, click the "View Exceptions" button to view those files. If any exception file is a valid trace file, then use the manual input interface to add it into the database. Otherwise, the user should remove the file.

Note: If there is any error during the import process, the error message is usually displayed in the command line window that accompanies the software window.

**The database**

The database is an elementary SQLite database containing two tables, along with three views.

The first table is the Sequence table, which contains most of the information about the sequence files. The other table is the Blast table, which contains *only* the first result of the BLASTn search on the sequence.

The three views are: LowQuality, Rerun, and Pursue.

Sequence table
------------------------
<table>
	<tr>
	<th>ID</th>
	<th>Plate</th>
	<th>Clone</th>
	<th>Primer</th>
	<th>Run_date</th>
	<th>Student</th>
	<th>Instructor</th>
	<th>Institution</th>
	<th>Quality</th>
	<th>FASTA</th>
	<th>ABI</th>
	<th>Pursue</th>
	<th>Comment</th>
	</tr>
	<tr>
	<td>2384</td>
	<td>2</td>
	<td>H5</td>
	<td>F</td>
	<td>2006-20-01</td>
	<td>NLA</td>
	<td>Bangera</td>
	<td>Bellevue College</td>
	<td>1</td>
	<td>Sequence Files/FASTA/PF002/H/Pf002_H5F_2006-20-01.fasta</td>
	<td>Sequence Files/ABI/PF002/H/Pf002_H5F_2006-20-01.abi</td>
	<td>0</td>
	<td></td>
	</tr>
</table>

<dl>
	<dt>ID</dt>
	<dd>Unique identifier of the sequence in the database. This information is used to correlate the sequence with its BLASTn search result in the Blast table.</dd>
	<dt>Plate</dt>
	<dd>The plate number of the bacteria clone.</dd>
	<dt>Clone</dt>
	<dd>The clone number of the bacteria clone.</dd>
	<dt>Primer</dt>
	<dd>The primer used for the sequencing reaction. F stands for SL1, R stands for SR2.</dd>
	<dt>Run_date</dt>
	<dd>The specific date when the sequence was read on the sequencing machine. When the sequencing machine is left to run overnight, it is possible that two different sequencing in that same batch will have different run date. This is because the run date extracted from the sequence file is the date of *that* individual file. For instance, both sample A and B are loaded into the sequencing machine at the same time, but sample A is sequenced on 11:30 PM of 2010-10-23, sample B is sequenced on 1:05 AM of 2010-10-24, then they will have different run date.</dd>
	<dt>Student</dt>
	<dd>The student who performed the procedures to sequence the clone.</dd>
	<dt>Instructor</dt>
	<dd>The instructor who oversee the student.</dd>
	<dt>Institution</dt>
	<dd>The place where the sequencing was done.</dd>
	<dt>Quality</dt>
	<dd>How good the sequence is. There are two possible values: 0 and 1. "1" means that the sequence is good enough for BLAST search, "0" means otherwise. Even if a record has a quality of "1", it does not mean that it will be perfect, but only that it might still give meaningful result in BLAST searches.</dd>
	<dt>FASTA</dt>
	<dd>Path to FASTA file of the sequence.</dd>
	<dt>ABI</dt>
	<dd>Path to ABI trace file of the sequence.</dd>
	<dt>Pursue</dt>
	<dd>Whether the sequence is selected for further investigation. "1" means yes, and "0" means no.</dd>
	<dt>Comment</dt>
	<dd>Any additional comment.</dd>
</dl>

Blast table
----------------
<table>
	<tr>
	<th>ID</th>
	<th>Genome</th>
	<th>Organism</th>
<<<<<<< HEAD
=======
	<th>E_value</th>
>>>>>>> develop
	<th>Query_from</th>
	<th>Query_to</th>
	<th>Subject_from</th>
	<th>Subject_to</th>
	<th>Identity</th>
	<th>Similarity</th>
	</tr>
	<tr>
	<td>2384</td>
	<td>CP002585.1</td>
	<td>Pseudomonas brassicacearum subsp. brassicacearum NFM421</td>
<<<<<<< HEAD
=======
	<td>1e-145</td>
>>>>>>> develop
	<td>100</td>
	<td>745</td>
	<td>40555456</td>
	<td>40556101</td>
	<td>99%</td>
	<td>F113</td>
	</tr>
</table>

<dl>
	<dt>ID</dt>
	<dd>Unique identifier for the sequence in this database. This value is used to link the BLASTn result with a sequence in the Sequence table.</dd>
	<dt>Genome</dt>
	<dd>Genome ID of the matching organism. The value used here is the Genebank ID, not the accession ID.</dd>
	<dt>Organism</dt>
	<dd>Name of the organism that this gene is found in, with the highest certainty by the BLASTn search.</dd>
<<<<<<< HEAD
=======
	<dt>E_value</dt>
	<dd>Expected matching sequences in the database just by chance.</dd>
>>>>>>> develop
	<dt>Query_from and Query_to</dt>
	<dd>Location of the matching hit on the query.</dd>
	<dt>Subject_from and Subject_to</dt>
	<dd>Location of the hit on the subject.</dd>
	<dt>Identity</dt>
	<dd>The percentage similar of the query to the match.</dd>
	<dt>Similarity</dt>
	<dd>Other organism that this sequence is also found in.</dd>
</dl>

**The program**
The program is written in Python. It requires BioPython to parse and handle the sequence files, as well as perform and parse the BLAST searches. The program uses Tkinter as the interface package.

<<<<<<< HEAD
In short summary, the program will load the list of files inside the Good and Bad folder, iterate through the files to extract the necessary information, and then add those information into the record database. Then it move the files into their designated folders. One pitfall of this approach is that the sequence files might be entered into the record database, only to have it being deleted at the target site.
=======
In short summary, the program will load the list of files inside the Good and Bad folder, iterate through the files to extract the necessary information, and then add those information into the record database. Then it move the files into their designated folders. One pitfall of this approach is that the sequence files might be entered into the record database, only to have it being deleted at the target site.
>>>>>>> develop
