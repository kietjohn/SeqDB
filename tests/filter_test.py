from nose.tools import *
from SeqDB.tools import abi_filter

def filter_test():
    '''Filter abi files'''
    file_list = ['test_data/file1.ab1',
                'test_data/file2.abi',
                'dummy.ab1',
                'test_data/blast/single_blast.xml',
                'no_extension',
                'test_data/file1.ab1.dummy']
    abi, non_abi = abi_filter(file_list)

    assert_equal(abi, file_list[0:3])
    assert_equal(non_abi, file_list[3:])
