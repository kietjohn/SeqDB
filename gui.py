from dbimport import *
from manual import *
import subprocess
import Tkinter
import tkFileDialog
import tkFont

class interface(Tkinter.Tk):


	def __init__(self,parent):
		Tkinter.Tk.__init__(self,parent)
		self.parent = parent
		self.customFont = tkFont.Font(family='Arial', size=12)
		self.initialize()


	def initialize(self):
		self.grid()

		self.label_xsep = 15
		self.label_ysep = 2
		self.button_xsep = 10
		self.button_ysep = 2
		self.button_bottom = 40

		# Main frame: contains input fields and functional buttons
	
		self.main_frame = Tkinter.Frame()
		self.main_frame.grid(column=0, row=0, sticky='EWNS',padx=10,pady=10)

		# Separator
		self.separator1 = Tkinter.Frame(self.main_frame,height=10)
		self.separator1.grid()

		# Step 1
		self.step1 = Tkinter.Frame(self.main_frame)
		self.step1.grid(sticky='EWNS')

		self.step1_label = Tkinter.Label(self.step1, text='Step 1: Copy sequence files into Input folder and check quality', anchor='w', font=self.customFont)
		self.step1_label.grid(column=0,row=0,columnspan=6,sticky='EW',pady=20)
		
		self.quality = Tkinter.Button(self.step1, text='1. Check Quality',command=self.check_procedure,font=self.customFont)
		self.quality.grid(column=1, row=1, sticky='EW')
		
		self.uncertain = Tkinter.Button(self.step1, text='2. View uncertains', command=self.OpenUncertain, font=self.customFont)
		self.uncertain.grid(column=3, row=1,sticky='EW')

		# Separator
		self.separator2 = Tkinter.Frame(self.main_frame,height=20)
		self.separator2.grid()

		# Step 2
		self.step2 = Tkinter.Frame(self.main_frame)
		self.step2.grid(sticky='EWNS')

		self.step2_label = Tkinter.Label(self.step2, text='Step 2: Input information for the files and import into database', anchor='w', font=self.customFont)
		self.step2_label.grid(column=0,row=0,columnspan=6,sticky='W', pady=20)

		self.plateLabel = Tkinter.Label(self.step2, text="Plate number (ex: PF003):", anchor='w', font=self.customFont)
		self.plateLabel.grid(column=0, row=1, columnspan=2, sticky='EW', padx=self.label_xsep,pady=self.label_ysep)

		self.instructorLabel = Tkinter.Label(self.step2,text='Instructor name:',anchor='w',font=self.customFont)
		self.instructorLabel.grid(column=0,row=2,columnspan=2,sticky='EW',padx=self.label_xsep,pady=self.label_ysep)

		self.instructorLabel = Tkinter.Label(self.step2,text='Sequencing Institution:',anchor='w',font=self.customFont)
		self.instructorLabel.grid(column=0,row=3,columnspan=2,sticky='EW',padx=self.label_xsep,pady=self.label_ysep)

		self.good_label = Tkinter.Label(self.step2, text='Good Sequence Folder:', anchor='w', font=self.customFont)
		self.good_label.grid(column=0, row=4, columnspan=2, sticky='EW', padx=self.label_xsep,pady=self.label_ysep)

		self.bad_label = Tkinter.Label(self.step2, text='Bad Sequence Folder:', anchor='w', font=self.customFont)
		self.bad_label.grid(column=0, row=5, columnspan=2, sticky='EW', padx=self.label_xsep,pady=self.label_ysep)

		self.database_path_label = Tkinter.Label(self.step2, text='Database path:', anchor='w', font=self.customFont)
		self.database_path_label.grid(column=0, row=6, columnspan=2, sticky='EW', padx=self.label_xsep,pady=self.label_ysep)

		# Text fields
		self.plate_number = Tkinter.StringVar()
		self.plate_entry = Tkinter.Entry(self.step2, textvariable=self.plate_number, font=self.customFont)
		self.plate_entry.grid(column=2, row=1, columnspan=4, sticky='EW', padx=self.button_xsep,pady=self.button_ysep)
		self.plate_entry.focus_set()
		self.plate_entry.selection_range(0, Tkinter.END)

		self.instructor = Tkinter.StringVar()
		self.instructor_entry = Tkinter.Entry(self.step2, textvariable = self.instructor, font=self.customFont)
		self.instructor_entry.grid(column=2, row=2, columnspan=4, sticky='EW', padx=self.button_xsep,pady=self.button_ysep)
		self.instructor.set('Bangera')

		self.institution = Tkinter.StringVar()
		self.institution_entry = Tkinter.Entry(self.step2, textvariable = self.institution, font=self.customFont)
		self.institution_entry.grid(column=2, row=3, columnspan=4, sticky='EW', padx=self.button_xsep,pady=self.button_ysep)
		self.institution.set('Bellevue College')

		self.good_folder_path = Tkinter.StringVar()
		self.good_folder = Tkinter.Entry(self.step2, textvariable = self.good_folder_path, font=self.customFont)
		self.good_folder.grid(column=2, row=4, columnspan=4, sticky='EW', padx=self.button_xsep,pady=self.button_ysep)
		self.good_folder_path.set('Input/Good')

		self.bad_folder_path = Tkinter.StringVar()
		self.bad_folder = Tkinter.Entry(self.step2,textvariable = self.bad_folder_path, font=self.customFont)
		self.bad_folder.grid(column=2, row=5, columnspan=4, sticky='EW', padx=self.button_xsep,pady=self.button_ysep)
		self.bad_folder_path.set('Input/Bad')

		self.database_path = Tkinter.StringVar()
		self.database = Tkinter.Entry(self.step2, textvariable = self.database_path, font=self.customFont)
		self.database.grid(column=2, row=6, columnspan=4, sticky='EW', padx=self.button_xsep,pady=self.button_ysep)
		self.database_path.set('Sequence Files/Seq.db')

		# Buttons

		self.good_folder_choose = Tkinter.Button(self.step2, text='...', command=self.browse_good_folder, font=self.customFont)
		self.good_folder_choose.grid(column=6, row=4, sticky='EW')

		self.bad_folder_choose = Tkinter.Button(self.step2, text='...', command=self.browse_bad_folder, font=self.customFont)
		self.bad_folder_choose.grid(column=6, row=5, sticky='EW')

		self.database_path_button = Tkinter.Button(self.step2, text='...', command=self.browse_database,width=2, font=self.customFont)
		self.database_path_button.grid(column=6, row=6, sticky='EW')

		self.open_log = Tkinter.Button(self.step2, text="Open log", command=self.OpenLog, font=self.customFont)
		self.open_log.grid(column=1, row=7, sticky='EW', pady=self.button_bottom)

		self.import_file_button = Tkinter.Button(self.step2, text="Import Files", command=self.import_procedure, font=self.customFont)
		self.import_file_button.grid(column=3, row=7, sticky='EW', pady=self.button_bottom)

		#Step 3
		self.step3 = Tkinter.Frame(self.main_frame)
		self.step3.grid(sticky='EWNS')

		self.step3_label = Tkinter.Label(self.step3, text="Step 3: View the exception files.", anchor='w', font=self.customFont)
		self.step3_label.grid(column=0, row=0, columnspan=6, sticky='W')

		self.open_exception = Tkinter.Button(self.step3, text="View Exceptions", command=self.OpenException, font=self.customFont)
		self.open_exception.grid(column=5, sticky='EW', pady=self.button_bottom)

		# Log frame
		self.log_frame = Tkinter.Frame(bd=2)
		self.log_frame.grid(column=1, row=0, sticky='EWNS', padx=10, pady=10)

		self.log_label = Tkinter.Label(self.log_frame, text='Transaction Log', anchor='w', font=self.customFont)
		self.log_label.grid(column=0, row=0, sticky='EW', padx=self.label_xsep,pady=self.label_ysep)

		self.scrollbar = Tkinter.Scrollbar(self.log_frame)
		self.scrollbar.grid(column=6, row=1, rowspan=3, sticky='ns')

		self.log_display = Tkinter.Text(self.log_frame,width=55, state='disabled', font=self.customFont, yscrollcommand=self.scrollbar.set)
		self.log_display.grid(column=0, row=1, columnspan=6, rowspan=3, sticky='EWNS', padx=10, pady=5)
		self.scrollbar.config(command=self.log_display.yview)

		self.close = Tkinter.Button(self.log_frame, text='Close', command=self.destroy, font=self.customFont)
		self.close.grid(column=5, row=5, sticky='EW', ipadx=5,pady=self.button_ysep)

		# Settings commands
		self.grid_columnconfigure(0,weight=1,minsize=150)
		self.grid_columnconfigure(1,weight=1,minsize=150)
		self.resizable(False,False)
		self.update()
		self.geometry(self.geometry())

	# Button commands
	def check_procedure(self):
		self.log_message = checkQuality()
		writeToLog(self.log_message)
		self.displayLog(self.log_message)

	def OpenUncertain(self):
		try:
			subprocess.check_call(['explorer', 'Input/Uncertain'])
		except:
			subprocess.check_call(['gnome-open', 'Input/Uncertain'])
	
	def OpenLog(self):
		os.startfile('log.txt')

	def browse_good_folder(self):
		self.dirname = tkFileDialog.askdirectory(parent = self.parent, initialdir = '/', title = 'Select Folder with Good sequences')
		if len(self.dirname) > 0:
			self.good_folder_path.set(self.dirname)

	def browse_bad_folder(self):
		self.dirname = tkFileDialog.askdirectory(parent=self.parent, title = 'Select Folder with Bad sequences')
		if len(self.dirname) > 0:
			self.bad_folder_path.set(self.dirname)

	def browse_database(self):
		self.database_file = tkFileDialog.askopenfilename(filetypes = ( ("Database File", "*.db") , ("SQL database","*.sql") ))
		self.database_path.set(self.database_file)
	
	def import_procedure(self):
		self.plate = self.plate_number.get()
		self.gd_source = self.good_folder.get()
		self.bd_source = self.bad_folder.get()
		self.db_path = self.database.get()
		self.instructor_name = self.instructor.get()
		self.institution_name = self.institution.get()
		
		self.good_source = self.get_directory(self.gd_source)
		self.bad_source= self.get_directory(self.bd_source)
		self.database_path = os.path.abspath(self.db_path)
		
		self.good_success, self.good_fail = import_files(self.good_source,self.plate,self.instructor_name,self.institution_name,1,self.database_path)
		self.bad_success, self.bad_fail = import_files(self.bad_source,self.plate,self.instructor_name,self.institution_name,0,self.database_path)
		self.log_message = log_process(self.good_success + self.bad_success, self.good_fail + self.bad_fail)
		writeToLog(self.log_message)
		self.displayLog(self.log_message)

	def displayLog(self, msg):
		self.log_display['state'] = 'normal'
		self.log_display.insert('end', msg)
		self.log_display['state'] = 'disabled'

	def OpenException(self):
		try:
			subprocess.check_call(['explorer', 'Exceptions'])
		except:
			subprocess.check_call(['gnome-open', 'Exceptions'])

	def get_directory(self, directory):
		self.path = os.path.abspath(directory)
		return os.path.join(self.path, '')

def manual_interface():
	man = manual_input(None)
	man.title('Manual input')
	man.mainloop()

# Initiates the application interface.
if __name__ == "__main__":
	app = interface(None)
	app.title('Database entry')
	app.mainloop()
