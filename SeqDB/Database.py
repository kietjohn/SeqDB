import parameters
import tools
import sqlite3 as sql
import os
import datetime


class Database(object):
    """The Database object holds the connection to the Sequence database,
    and can be initiated with a path to the database:
    >> database = Database('../Sequence Files/Seq.db')

    """

    def __init__(self, database_path):
        self.path = database_path
        self.con = sql.connect(self.path)

    def __enter__(self):
        self.con = sql.connect(self.path)
        return self.con

    def __exit__(self, type, value, tb):
        self.con.close()

    def insert_sequence(self, entries, dev=False):
        """Insert a list of entries into the given table in the  database.
        The entries will be inserted one-by-one. After completing the
        interation, the function will return the table name, the  number of
        successfull and failed insertion, as well as a list of the failed
        entries."""

        self.success = []
        self.fail = []
        self.table = 'Sequence'

        with self.con:
            self.con.text_factory = str
            self.cur = self.con.cursor()
            for x in entries:
                try:
                    self.cur.execute(parameters.SQL_COMMAND[self.table], x[2])
                    self.success.append(x)
                except:
                    self.fail.append(x)

        # Move the files inserted to their places, and failed ones to exception
        # Then check if everything has been moved
        if dev is False:
            report = self.process_files(self.success, self.fail)

        print "Import complete. %d inserted, %d failed." % (len(self.success),len(self.fail))
        return self.success, self.fail

    def process_files(self, success, fail):
        '''Moves the files entered into the database into their respective
        places. The successfully entered will go to its designated folder,
        while the failed ones will go to the exceptions folder. All operations
        will be double-checked with the list supplied by database.'''

        logging = {'success': [], 'fail': []}

        # Handle the success
        for x in success:
            fasta = x[1].replace('ab1', 'fasta')
            abi = x[1]
            tools.move_single_file(fasta, x[2][8])
            tools.move_single_file(abi, x[2][9])
            logging['success'].append(x)

        # Handle the fail
        self.fail_folder = os.path.join('Exceptions', str(datetime.date.today()))
        for x in fail:
            fasta = x[1].replace('ab1', 'fasta')
            abi = x[1]
            tools.move_single_file(abi, self.fail_folder)
            os.remove(fasta)
            logging['fail'].append(x)

        return logging

    def get_to_blast(self):
        with self.con:
            self.con.text_factory = str
            self.cur = self.con.cursor()
            self.cur.execute('SELECT * FROM toBlast')
            toBlast_list = self.cur.fetchall()
        return toBlast_list

    def insert_blast(self, entries, dev=False):
        """Insert a list of entries into the given table in the  database.
        The entries will be inserted one-by-one. After completing the
        interation, the function will return the table name, the  number of
        successfull and failed insertion, as well as a list of the failed
        entries."""

        self.success = []
        self.fail = []
        self.table = 'Blast'

        with self.con:
            self.con.text_factory = str
            self.cur = self.con.cursor()
            for x in entries:
                try:
                    self.cur.execute(parameters.SQL_COMMAND[self.table], x)
                    self.success.append(x)
                except:
                    self.fail.append(x)

        print "Import complete. %d inserted, %d failed." % (len(self.success),len(self.fail))
        return self.success, self.fail
