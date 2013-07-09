from SeqDB import *
from SeqDB.TraceFile import *
import os
import sys
import subprocess
import datetime
from PyQt4 import QtCore, QtGui


class MainWindow(QtGui.QMainWindow):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        self.ui = gui.Ui_MainWindow()
        self.ui.setupUi(self)
        self.ui.statusbar.showMessage('Greetings!')

        # Manage the actions of the widgets
        QtCore.QObject.connect(self.ui.file_input, QtCore.SIGNAL("dropped"), self.main_add_file)
        QtCore.QObject.connect(self.ui.main_import, QtCore.SIGNAL("clicked()"), self.import_sequences)
        QtCore.QObject.connect(self.ui.view_exceptions, QtCore.SIGNAL("clicked()"), self.view_exceptions)
        QtCore.QObject.connect(self.ui.view_log, QtCore.SIGNAL("clicked()"), self.show_log)
        QtCore.QObject.connect(self.ui.main_clear_entry, QtCore.SIGNAL("clicked()"), self.main_clear_entry)

    # Main input interface, automatically grabs info from the files
    def main_add_file(self):
        '''Handles adding files to the list window, and then display the path
        to the files.'''
        self.ui.file_input.clear()
        for url in self.ui.file_input.file_list:
            if os.path.exists(url):
                print(url)
                self.ui.file_input.addItem(url)

    def main_clear_entry(self):
        '''Delete all entries from the file list'''
        self.ui.file_input.file_list[:] = []
        self.ui.file_input.clear()
        self.ui.statusbar.showMessage('File list cleared.')

    def import_sequences(self):
        '''Get information from the text fields.'''
        self.file_list = self.ui.file_input.file_list
        self.plate = int(self.ui.plate_entry.text())
        self.instructor = str(self.ui.instructor_entry.text())
        self.institution = str(self.ui.institution_entry.text())

        # Initiate the sequence cache, then insert the entries obtained to db
        print self.plate, self.instructor, self.institution
        self.cache = TraceFile(self.file_list,
            self.plate,
            self.instructor,
            self.institution)
        cache_report = self.cache.report()
        success, fail = db.insert_sequence(self.cache.entries)
        self.logging(cache_report, success, fail)
        self.cache.clean_up()

        self.ui.file_input.clear()
        self.ui.file_input.file_list[:] = []
        del(self.cache)

        self.ui.statusbar.showMessage('Import complete.')

    def view_exceptions(self):
        try:
            subprocess.check_call(['explorer', 'Exceptions'])
        except:
            subprocess.check_call(['gnome-open', 'Exceptions'])

        self.ui.statusbar.showMessage('Open Exceptions folder.')

    # Manual input interface. Use has to input some info for the files
    def get_file_path(self):
        self.ui.man_file_path_entry.setText(QtGui.QFileDialog.getOpenFileName())
        self.ui.statusbar.showMessage('File path obtained.')

    def logging(self, cache_report, success, fail):
        msg = [cache_report,
            str(datetime.datetime.now()),
            'Files that have been successfully imported:'
            ]

        for item in success:
            path = item[1]
            try:
                target = '_'.join([str(i) for i in item[2][0:5]])
            except:
                pass
            line = "\t- %s imported as %s" % (path, target)
            msg.append(line)
        msg.append('Files that have failed:')
        for item in fail:
            line = "\t- %s failed to import." % path
            msg.append(line)
        msg.append('=====================')

        text = '\n'.join(msg)
        self.ui.log.append(text)

        if not os.path.isfile('log.txt'):
            with open('log.txt', 'w') as handle:
                handle.write('Log file for importing program')

        with open('log.txt', 'a') as handle:
            handle.write(text)

    def show_log(self):
        os.startfile('log.txt')

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    db = Database.Database(parameters.DB_PATH)
    myapp = MainWindow()
    myapp.show()
    sys.exit(app.exec_())
