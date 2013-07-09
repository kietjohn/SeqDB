# -*- coding: utf-8 *-*
import os
import shutil


def abi_filter(file_list):
    abi = []
    non_abi = []

    for item in file_list:
        file_name, file_ext = os.path.splitext(item)
        if file_ext == '.ab1' or file_ext == '.abi':
            abi.append(item)
        else:
            non_abi.append(item)

    return abi, non_abi


def listdir_joined(path):
    '''Generate a list of file paths for files in the specified folder.'''
    file_list = []
    for item in os.listdir(path):
        file_path = os.path.join(path, item)
        if bool(os.path.isfile(file_path)):
            file_list.append(file_path)
    return sorted(file_list)


def generate_single_fasta(src):
    '''Generate a fasta file for the sequence file.'''
    file_name, ext = os.path.splitext(src)
    from Bio import SeqIO
    SeqIO.convert(src, 'abi', file_name + '.fasta', 'fasta')


def copy_list_to_folder(src, dst):
    '''Copy files to a destination folder. If that folder does not exist,
    a new folder with that name will be created. The files should be supplied
    as a list of paths.'''
    # Check if the name supplied is a folder, then check its existence
    if not bool(os.path.splitext(dst)[1]):
        if not os.path.isdir(dst):
            os.makedirs(dst)
    else:
        par = os.path.dirname(dst)
        if not os.path.isdir(par):
            os.makedirs(par)

    for item in src:
        shutil.copy2(item, dst)


def move_file(src, dst):
    '''Move files to a destination folder. If that folder does not exist,
    a new folder with that name will be created. The files should be supplied
    as a list of paths.'''
    # Check if the name supplied is a folder, then check its existence
    copy_list_to_folder(src, dst)
    for item in src:
        os.remove(item)


def copy_file(src, dst):
    '''Copy a file to a new destination, and rename it according
    to the path provided.'''
    # Check if the name supplied is a folder, then check its existence
    par = os.path.dirname(dst)
    if not os.path.isdir(par):
        os.makedirs(par)

    shutil.copy2(src, dst)


def move_single_file(src, dst):
    '''Move a single file to a destination file. This is essentially the
    copy_file function followed by removing the src file.'''
    copy_file(src, dst)
    os.remove(src)
