from nose.tools import *
from SeqDB.Sequence import *

import os
import shutil


def sequence_object_test():
    '''Sequence objects parsing

    Given several sequences, load them into the object, then compare the output
    record with manual inputs.
    '''

    src = 'tests/test_data/data_src/parse'
    dst = 'tests/test_data/parsing'
    test_files = sorted(prep_data(src, dst))
    print test_files

    bad1 = Sequence(test_files[0], 3, 'Bangera', 'Bellevue College', comment='Test cases')
    assert_equal(bad1.record, [True, test_files[0], (3, 'C5', 'R', '2012-03-09', 'IW', 'Bangera', 'Bellevue College', 0, 'Sequence Files/FASTA/PF003/C/PF003_C5R_2012-03-09.fasta', 'Sequence Files/ABI/PF003/C/PF003_C5R_2012-03-09.ab1', 'Test cases')])

    bad2= Sequence(test_files[1], 1, 'Bangera', 'Bellevue College')
    assert_equal(bad2.record, [True, test_files[1], (1, 'C12', 'R', '2012-03-30', 'MB', 'Bangera', 'Bellevue College', 0, 'Sequence Files/FASTA/PF001/C/PF001_C12R_2012-03-30.fasta', 'Sequence Files/ABI/PF001/C/PF001_C12R_2012-03-30.ab1', None)])

    fail1 = Sequence(test_files[2], 3, 'Bangera', 'Bellevue College', comment='Test cases')
    assert_equal(fail1.record, [False, test_files[2], 'POSITIVE.ab1'])

    good = Sequence(test_files[3], 1, 'Bangera', 'Bellevue College', comment='Test')
    assert_equal(good.record, [True, test_files[3], (1, 'A1', 'R', '2012-03-30', 'MB', 'Bangera', 'Bellevue College', 1, 'Sequence Files/FASTA/PF001/A/PF001_A1R_2012-03-30.fasta', 'Sequence Files/ABI/PF001/A/PF001_A1R_2012-03-30.ab1', 'Test')])

    fail2 = Sequence(test_files[4], 6, 'Bangera', 'Bellevue College')
    assert_equal(fail2.record, [False, test_files[4], 'C8B-SAS.ab1'])

    manual = Sequence(test_files[4], 3, 'Bangera', 'Bellevue College', clone='C8', primer='R', student='SAS', comment='Corrected Insert')
    assert_equal(manual.record, [True, test_files[4], (3, 'C8', 'R', '2012-03-30', 'SAS', 'Bangera', 'Bellevue College', 0, 'Sequence Files/FASTA/PF003/C/PF003_C8R_2012-03-30.fasta', 'Sequence Files/ABI/PF003/C/PF003_C8R_2012-03-30.ab1', 'Corrected Insert')])

    shutil.rmtree(dst)

def prep_data(src, dst):
    if not os.path.isdir(dst):
        os.mkdir(dst)
    for item in [os.path.join(src, x) for x in os.listdir(src)]:
        shutil.copy2(item, dst)
    return [os.path.join(dst, x) for x in os.listdir(dst)]