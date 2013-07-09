# -*- coding: utf-8 -*-
import os
import shutil
import datetime
import tools

from Sequence import Sequence


class TraceFile(object):
    '''An object to handle working with all the ab1 files. Upon initation,
    it will generate a list of Sequence object corresponding to the list of
    files, and a list of entries that can be inserted into the database.
    There are also list variables to hold the files that don't get entered
    into the database.'''

    def __init__(self, file_list, *args, **kwargs):
        self.temp_dir = 'temp/trace_files'

        # Copy the files provided into the temp directory
        tools.copy_list_to_folder(file_list, self.temp_dir)
        self.file_list = tools.listdir_joined(self.temp_dir)

        self.trace_file, self.non_abi = tools.abi_filter(self.file_list)

        self.seq_list = []
        self.entries = []

        self.parse_failed = []
        self.exceptions = []
        self.fail_folder = os.path.join('Exceptions', str(datetime.date.today()))

        for item in self.trace_file:
            try:
                holder = Sequence(item, *args, **kwargs)
                self.seq_list.append(holder)
            except:
                self.parse_failed.append(item)
                print "Failed to parse ", item

        # Select the valid files for entry
        for seq in self.seq_list:
            if seq.record[0] is True:
                self.entries.append(seq.record)
            elif seq.record[0] is False:
                self.exceptions.append(seq.record)
                print seq.record[1], ' is not a valid file.'

    def clean_up(self):
        '''Move exception files into the exception folder.'''
        dont_rename = self.parse_failed + self.non_abi
        tools.move_file(dont_rename, self.fail_folder)

        for entry in self.exceptions:
            src = entry[1]
            dst = os.path.join(self.fail_folder, entry[2])
            tools.copy_file(src, dst)
            os.remove(entry[1])
        assert(os.listdir(self.temp_dir), [])
        shutil.rmtree(self.temp_dir)

    def report(self):
        failed_to_parse = [x + ' cannot be read.' for x in self.parse_failed]
        not_abi = [x + ' is not an ABI file.' for x in self.non_abi]

        invalid = [x[1] + ' is does not follow the format.' for x in self.exceptions]
        valid = [x[1] + ' is read successfully.' for x in self.entries]

        msg = [str(datetime.datetime.now())] + failed_to_parse + not_abi + invalid + valid
        return '\n'.join(msg)