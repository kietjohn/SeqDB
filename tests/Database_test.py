from nose.tools import *
from SeqDB.Database import *
from SeqDB.parameters import schema

import os


def seq_insert_test():
    '''Insert list of sequences'''
    db_path = 'tests/test_data/database/insert.db'
    db = create_database(db_path)

    entries = [
        [True, 'tests/test_data/database/trace_file/file1.ab1', (3, 'C5', 'R', '2012-03-09', 'IW', 'Bangera', 'Bellevue College', 0, 'Sequence Files/FASTA/PF003/C/PF003_C5R_2012-03-09.fasta', 'Sequence Files/ABI/PF003/C/PF003_C5R_2012-03-09.abi', 'Test cases')],
        [True, 'tests/test_data/database/trace_file/file2.ab1', (1, 'C12', 'R', '2012-03-30', 'MB', 'Bangera', 'Bellevue College', 0, 'Sequence Files/FASTA/PF001/C/PF001_C12R_2012-03-30.fasta', 'Sequence Files/ABI/PF001/C/PF001_C12R_2012-03-30.abi', None)],
        [True, 'tests/test_data/database/trace_file/file3.ab1', (1, 'A1', 'R', '2012-03-30', 'MB', 'Bangera', 'Bellevue College', 1, 'Sequence Files/FASTA/PF001/A/PF001_A1R_2012-03-30.fasta', 'Sequence Files/ABI/PF001/A/PF001_A1R_2012-03-30.abi', 'Test')],
        [True, 'tests/test_data/database/trace_file/file4.ab1', (3, 'C8', 'R', '2012-03-30', 'SAS', 'Bangera', 'Bellevue College', 0, 'Sequence Files/FASTA/PF003/C/PF003_C8R_2012-03-30.fasta', 'Sequence Files/ABI/PF003/C/PF003_C8R_2012-03-30.abi', 'Corrected Insert')],
        [True, 'tests/test_data/database/trace_file/file5.ab1', (1, 'A1', 'R', '2012-03-30', 'MB', 'Bangera', 'Bellevue College', 1, 'Sequence Files/FASTA/PF001/A/PF001_A1R_2012-03-30.fasta', 'Sequence Files/ABI/PF001/A/PF001_A1R_2012-03-30.abi', 'Test')]
        ]
    result = [
        (3, 'C5', 'R', '2012-03-09', 'IW', 'Bangera', 'Bellevue College', 0, 'Sequence Files/FASTA/PF003/C/PF003_C5R_2012-03-09.fasta', 'Sequence Files/ABI/PF003/C/PF003_C5R_2012-03-09.abi', None, 'Test cases'),
        (1, 'C12', 'R', '2012-03-30', 'MB', 'Bangera', 'Bellevue College', 0, 'Sequence Files/FASTA/PF001/C/PF001_C12R_2012-03-30.fasta', 'Sequence Files/ABI/PF001/C/PF001_C12R_2012-03-30.abi', None, None),
        (1, 'A1', 'R', '2012-03-30', 'MB', 'Bangera', 'Bellevue College', 1, 'Sequence Files/FASTA/PF001/A/PF001_A1R_2012-03-30.fasta', 'Sequence Files/ABI/PF001/A/PF001_A1R_2012-03-30.abi', None, 'Test'),
        (3, 'C8', 'R', '2012-03-30', 'SAS', 'Bangera', 'Bellevue College', 0, 'Sequence Files/FASTA/PF003/C/PF003_C8R_2012-03-30.fasta', 'Sequence Files/ABI/PF003/C/PF003_C8R_2012-03-30.abi', None, 'Corrected Insert'),
        (1, 'A1', 'R', '2012-03-30', 'MB', 'Bangera', 'Bellevue College', 1, 'Sequence Files/FASTA/PF001/A/PF001_A1R_2012-03-30.fasta', 'Sequence Files/ABI/PF001/A/PF001_A1R_2012-03-30.abi', None, 'Test')
        ]

    db.insert_sequence(entries, dev=True)

    # Obtain results inserted into db, and turns into form of template
    db_entered = []
    with db.con:
        cur = db.con.cursor()
        cur.execute("SELECT * FROM Sequence;")
        holder = cur.fetchall()
        for x in holder:
            a = list(x)
            del a[0]
            db_entered.append(tuple(a))

    assert_equal(db_entered, result)

    db.con.close()
    os.remove(db_path)


def seq_insert_fail_test():
    '''Reject weird items'''
    db_path = 'tests/test_data/database/fail.db'
    db = create_database(db_path)

    entries = [
        [True, 'file_path', ('Missing', 'entries', 'Test fail 1')],
        [True, 'file_path', ('A1', 'R', 'date', 'MB', 'Bangera', 'BC', 1, 'Fasta_path', 'Abi_path', 'Test failure2')],
        [True, 'file_path', (1, 'A1', 'R', 'date', 'MB', 'Bangera', 'BC', 1, 'Fasta_path', 'Abi_path', 'Extra', 'Test failure3')],
        [True, 'file_path', ('')]
        ]

    db.insert_sequence(entries, dev=True)

    assert_equal(db.fail, entries)
    os.remove(db_path)


def create_database(db_path):

    try:
        os.remove(db_path)
    except:
        pass

    db = Database(db_path)
    with db.con as con:
        cur = con.cursor()
        cur.executescript(schema)
    return db
