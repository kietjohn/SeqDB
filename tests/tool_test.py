# -*- coding: utf-8 *-*
from nose.tools import *
from SeqDB.tools import *
import os


def multiple_files_move_test():
    '''Move a list of files

    Prepare a list of file for testing. Move the files to a temp folder,
    assert that the files are in that folder, then move them out, assert again
    that the file has been moved. Delete the folder.
    '''

    start = 'tests/test_data/tool_test/move_files'
    end = 'tests/test_data/tool_test/move_files/end'
    test_files = prep_data(start)
    file_list = [os.path.basename(item) for item in test_files]

    move_file(test_files, end)
    result = listdir_joined(end)
    result_file_list = [os.path.basename(item) for item in result]
    assert_equal(file_list, result_file_list)
    assert_equal(listdir_joined(start), [])
    print "Files are moved successfully."

    move_file(result, start)
    assert_equal(os.listdir(end), [])
    print "Files are put back successfully."

    shutil.rmtree(start)


def copy_to_file_test():
    '''Copy and rename files

    Prepare a list of files, take the first one for test. Move the file into
    rename folder, check content of that folder. Move file back out, delete
    folder.
    '''

    start = 'tests/test_data/tool_test/rename'
    end = 'tests/test_data/tool_test/rename/end'
    os.makedirs(end)
    test_files = prep_data(start)[0:10]
    names = [str(x) + '.ab1' for x in xrange(0, 10)]
    end_names = [os.path.join(end, x) for x in names]

    file_list_and_names = zip(test_files, end_names)

    for item in file_list_and_names:
        copy_file(*item)
    assert_equal(sorted(os.listdir(end)), names)

    shutil.rmtree(start)


def generate_fasta_test():
    '''Generate fasta files

    Prepare a list of files. Generate a fasta file for each, assert that the
    fasta files have been generated. Then deleted any new file.
    '''

    test_dir = 'tests/test_data/tool_test/generate_fasta'

    test_files = prep_data(test_dir)

    for item in test_files:
        generate_single_fasta(item)

    abi_list = [os.path.basename(item) for item in test_files]
    fasta_list = [item.replace('.ab1', '.fasta') for item in abi_list]
    combined_list = sorted(abi_list + fasta_list)
    generated_files = sorted(os.listdir(test_dir))

    assert_equal(generated_files, combined_list)

    print listdir_joined(test_dir)

    shutil.rmtree(test_dir)


def prep_data(dst):
    src = 'tests/test_data/data_src/raw'
    src_list = listdir_joined(src)

    if os.path.isdir(dst):
        shutil.rmtree(dst)
        os.makedirs(dst)
    else:
        os.makedirs(dst)
    for item in src_list:
        shutil.copy2(item, dst)
    return sorted(listdir_joined(dst))
