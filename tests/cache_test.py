# -*- coding: utf-8 *-*
from nose.tools import *
from SeqDB.TraceFile import TraceFile
import os
import shutil


def prep_data(src, dst):
    if not os.path.isdir(dst):
        os.mkdir(dst)
    for item in [os.path.join(src, x) for x in os.listdir(src)]:
        shutil.copy2(item, dst)
    return [os.path.join(dst, x) for x in os.listdir(dst)]


def temp_file_test():
    '''Initate temp cache

    Given a set of data, load them into TraceFile object, and check if the files
    have been copied into a temp folder.
    '''

    src = 'tests/test_data/data_src/raw'
    dst = 'tests/test_data/cache'
    cache_path = 'temp/trace_files'
    test_files = prep_data(src, dst)

    cache = TraceFile(test_files, 1, 'Bangera', 'Bellevue College')

    cache_file_list = os.listdir(cache_path)
    temp_file_list = [os.path.join(cache_path, x) for x in cache_file_list]
    temp_file_list = [x for x in temp_file_list if 'fasta' not in x]

    print 'cache.file_list: ', cache.file_list

    assert_equal(sorted(cache.file_list), sorted(temp_file_list))

    shutil.rmtree(dst)
    shutil.rmtree(cache_path)


def temp_record_test():
    '''Temp file entries

    Given a set of data, load them into TraceFile object, and comapre the record
    generated with manual inputs.
    '''

    src = 'tests/test_data/data_src/parse'
    dst = 'tests/test_data/cache_parse'
    test_files = prep_data(src, dst)

    entries_list = sorted([
        [True, 'temp/trace_files/file1.ab1', (1, 'C5', 'R', '2012-03-09', 'IW', 'Bangera', 'Bellevue College', 0, 'Sequence Files/FASTA/PF001/C/PF001_C5R_2012-03-09.fasta', 'Sequence Files/ABI/PF001/C/PF001_C5R_2012-03-09.ab1', None)],
        [True, 'temp/trace_files/file2.ab1', (1, 'C12', 'R', '2012-03-30', 'MB', 'Bangera', 'Bellevue College', 0, 'Sequence Files/FASTA/PF001/C/PF001_C12R_2012-03-30.fasta', 'Sequence Files/ABI/PF001/C/PF001_C12R_2012-03-30.ab1', None)],
        [False, 'temp/trace_files/file3.ab1', 'POSITIVE.ab1'],
        [True, 'temp/trace_files/file4.ab1', (1, 'A1', 'R', '2012-03-30', 'MB', 'Bangera', 'Bellevue College', 1, 'Sequence Files/FASTA/PF001/A/PF001_A1R_2012-03-30.fasta', 'Sequence Files/ABI/PF001/A/PF001_A1R_2012-03-30.ab1', None)],
        [False, 'temp/trace_files/file5.ab1', 'C8B-SAS.ab1']
        ])

    good_list = [
        [True, 'temp/trace_files/file1.ab1', (1, 'C5', 'R', '2012-03-09', 'IW', 'Bangera', 'Bellevue College', 0, 'Sequence Files/FASTA/PF001/C/PF001_C5R_2012-03-09.fasta', 'Sequence Files/ABI/PF001/C/PF001_C5R_2012-03-09.ab1', None)],
        [True, 'temp/trace_files/file2.ab1', (1, 'C12', 'R', '2012-03-30', 'MB', 'Bangera', 'Bellevue College', 0, 'Sequence Files/FASTA/PF001/C/PF001_C12R_2012-03-30.fasta', 'Sequence Files/ABI/PF001/C/PF001_C12R_2012-03-30.ab1', None)],
        [True, 'temp/trace_files/file4.ab1', (1, 'A1', 'R', '2012-03-30', 'MB', 'Bangera', 'Bellevue College', 1, 'Sequence Files/FASTA/PF001/A/PF001_A1R_2012-03-30.fasta', 'Sequence Files/ABI/PF001/A/PF001_A1R_2012-03-30.ab1', None)]
        ]

    test_cache = TraceFile(test_files, 1, 'Bangera', 'Bellevue College')
    generated_entries_list = sorted([seq.record for seq in test_cache.seq_list])

    assert_equal(generated_entries_list, entries_list)
    assert_equal(sorted(test_cache.entries), sorted(good_list))

    shutil.rmtree(dst)
    shutil.rmtree(test_cache.temp_dir)


def cache_exceptions_handling_test():
    '''Exception files handling

    Prepare a mix of files, and generate a TraceFile object consisting of those
    files. Then check if the non-abi files are filtered into the non_abi variable,

    Perform clean_up() of TraceFile, Assert that the files are moved to the
    Exceptions folder (set within the test folder).

    Assert that the file lists match, then move the files back into original folder.
    '''
    src = 'tests/test_data/data_src/mixed'
    dst = 'tests/test_data/mixed_files'
    file_list = prep_data(src, dst)

    non_abi = []
    for x in file_list:
        if 'abi' not in x:
            if 'ab1' not in x:
                non_abi.append(x)
    non_abi_file_names = sorted([os.path.basename(x) for x in non_abi])

    test_cache = TraceFile(file_list, 1, 'Bangera', 'Bellevue College')
    cache_non_abi = sorted([os.path.basename(x) for x in test_cache.non_abi])
    assert_equal(cache_non_abi, non_abi_file_names)

    test_cache.fail_folder = 'tests/test_data/mixed_files/Exceptions'
    test_cache.clean_up()
    fail_list = [os.path.join(test_cache.fail_folder,item) for item in os.listdir(test_cache.fail_folder)]
    cleaned_up_list = [os.path.join(test_cache.fail_folder, x) for x in non_abi_file_names]
    assert_equal(sorted(fail_list), sorted(cleaned_up_list))

    shutil.rmtree(test_cache.fail_folder)
    shutil.rmtree(dst)