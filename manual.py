import Tkinter, tkFileDialog, tkFont
from dbimport import *
import os

class manual_input(Tkinter.Tk):


	def __init__(self,parent):
		Tkinter.Tk.__init__(self,parent)
		self.parent = parent
		self.customFont = tkFont.Font(family='Arial', size=12)
		self.initialize()


	def initialize(self):
		self.grid()

		self.entries = []

		self.label_xsep = 15
		self.label_ysep = 2
		self.button_xsep = 10
		self.button_ysep = 2
		self.button_bottom = 40
		self.customFont = tkFont.Font(family='Arial', size=12)

		self.main_frame = Tkinter.Frame()
		self.main_frame.grid()

		self.path = Tkinter.Label(self.main_frame, text="File:", anchor='w', font=self.customFont)
		self.path.grid(column=0, row=0, columnspan=2, sticky='EW', padx=self.label_xsep,pady=self.label_ysep)
		
		self.plateLabel = Tkinter.Label(self.main_frame, text="Plate number:", anchor='w', font=self.customFont)
		self.plateLabel.grid(column=0, row=1, columnspan=2, sticky='EW', padx=self.label_xsep,pady=self.label_ysep)

		self.cloneLabel = Tkinter.Label(self.main_frame, text="Clone number:", anchor='w', font=self.customFont)
		self.cloneLabel.grid(column=0, row=2, columnspan=2, sticky='EW', padx=self.label_xsep,pady=self.label_ysep)

		self.primerLabel = Tkinter.Label(self.main_frame, text="Primer:", anchor='w', font=self.customFont)
		self.primerLabel.grid(column=0, row=3, columnspan=2, sticky='EW', padx=self.label_xsep,pady=self.label_ysep)

		self.qualityLabel = Tkinter.Label(self.main_frame, text="Quality:", anchor='w', font=self.customFont)
		self.qualityLabel.grid(column=0, row=4, columnspan=2, sticky='EW', padx=self.label_xsep,pady=self.label_ysep)

		self.instructorLabel = Tkinter.Label(self.main_frame,text='Instructor name:',anchor='w',font=self.customFont)
		self.instructorLabel.grid(column=0,row=5,columnspan=2,sticky='EW',padx=self.label_xsep,pady=self.label_ysep)

		self.instructorLabel = Tkinter.Label(self.main_frame,text='Sequencing Institution:',anchor='w',font=self.customFont)
		self.instructorLabel.grid(column=0,row=6,columnspan=2,sticky='EW',padx=self.label_xsep,pady=self.label_ysep)

		#Entry fields
		self.file_path = Tkinter.StringVar()
		self.file_entry = Tkinter.Entry(self.main_frame, textvariable=self.file_path, font=self.customFont)
		self.file_entry.grid(column=2, row=0, columnspan=4, sticky='EW', padx=self.button_xsep,pady=self.button_ysep)

		self.plate_number = Tkinter.StringVar()
		self.plate_entry = Tkinter.Entry(self.main_frame, textvariable=self.plate_number, font=self.customFont)
		self.plate_entry.grid(column=2, row=1, columnspan=4, sticky='EW', padx=self.button_xsep,pady=self.button_ysep)

		self.clone = Tkinter.StringVar()
		self.clone_entry = Tkinter.Entry(self.main_frame, textvariable=self.clone, font=self.customFont)
		self.clone_entry.grid(column=2, row=2, columnspan=4, sticky='EW', padx=self.button_xsep,pady=self.button_ysep)

		self.primer = Tkinter.StringVar()
		self.primer_entry = Tkinter.Entry(self.main_frame, textvariable=self.primer, font=self.customFont)
		self.primer_entry.grid(column=2, row=3, columnspan=4, sticky='EW', padx=self.button_xsep,pady=self.button_ysep)

		self.quality = Tkinter.StringVar()
		self.quality_entry = Tkinter.Entry(self.main_frame, textvariable=self.quality, font=self.customFont)
		self.quality_entry.grid(column=2, row=4, columnspan=4, sticky='EW', padx=self.button_xsep, pady=self.button_ysep)

		self.instructor = Tkinter.StringVar()
		self.instructor_entry = Tkinter.Entry(self.main_frame, textvariable=self.instructor, font=self.customFont)
		self.instructor_entry.grid(column=2, row=5, columnspan=4, sticky='EW', padx=self.button_xsep,pady=self.button_ysep)
		self.instructor.set('Bangera')

		self.institution = Tkinter.StringVar()
		self.institution_entry = Tkinter.Entry(self.main_frame, textvariable=self.institution, font=self.customFont)
		self.institution_entry.grid(column=2, row=6, columnspan=4, sticky='EW', padx=self.button_xsep,pady=self.button_ysep)
		self.institution.set('Bellevue College')

		# Button
		self.file_path_button = Tkinter.Button(self.main_frame, text='...', command=self.browse_file,width=2, font=self.customFont)
		self.file_path_button.grid(column=6, row=0, sticky='EW')
		
		self.add_to_list = Tkinter.Button(self.main_frame, text='Add to list', command=self.add_list, width=2, font=self.customFont)
		self.add_to_list.grid(column=1, row=7, sticky='EW', pady=self.button_bottom)

		self.file_import_button = Tkinter.Button(self.main_frame, text='Import files', command=self.importing, width=2, font=self.customFont)
		self.file_import_button.grid(column=3, row=7, sticky='EW', pady=self.button_bottom)

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

		self.grid_columnconfigure(0,weight=1,minsize=150)
		self.grid_columnconfigure(1,weight=1,minsize=150)
		self.resizable(False,False)
		self.update()
		self.geometry(self.geometry())

	def browse_file(self):
		self.file_name = tkFileDialog.askopenfilename(filetypes = ( ("Trace file", "*.abi") , ("All Files","*.*") ))
		self.file_path.set(os.path.abspath(self.file_name))

	def add_list(self):
		self.plate = str(self.plate_number.get())
		self.clone_num = str(self.clone.get())
		self.primer = str(self.primer.get())
		self.quality = str(self.quality.get())
		self.instructor_name = str(self.instructor.get())
		self.institution_name = str(self.institution.get())
		
		self.author, self.run_date, self.fasta_path, self.abi_path = file_handling(str(self.file_path.get()),self.plate,self.clone_num)

		self.entries.append( (self.plate[-1], self.clone_num, self.primer, self.run_date, self.author, self.instructor, self.institution, self.quality, self.fasta_path, self.abi_path) )
		self.msg = """File %s
		Plate: %s
		Clone: %s
		Primer: %s
		Quality: %s
		Instructor: %s
		Institution: %s
		""" % self.file_path, self.plate, self.clone, self.primer, self.quality, self.instructor_name, self.institution_name
		self.displayLog(self.msg)

	def import_procedure(self):
		try:
			add_to_database(self.entries,interface.database_path.get())
			self.log_message = 'The files have been successfully imported.'
			writeToLog(self.log_message)
			self.displayLog(self.log_message)
		except:
			print "Failed to import files."
		
		self.entries = []

	def importing(self):
		pass

	def displayLog(self, msg):
		self.log_display['state'] = 'normal'
		self.log_display.insert('end', msg)
		self.log_display['state'] = 'disabled'

	def get_directory(self, directory):
		self.path = os.path.abspath(directory)
		return os.path.join(self.path, '')

def file_handling(file_path,plate,clone):
	handle = SeqIO.read(file_path,'abi')

	info = handle.__str__().split('\n')
	for line in info:
		if 'run_finish' in line:
			trunk, info = line.split('=')
			run_date, run_time = info.split(' ')

	try:
		sample_name, author = handle.id.split('-')
	except :
		sample_name = read.id
		author = 'none'

	#Output an ABI and a FASTA file
	file_name = os.path.basename(file_path)
	abi_dest = 'Sequence Files/ABI/' + plate + '/' + clone[0] + '/'
	abi_path = abi_dest + file_name + '.abi'
	fasta_dest = 'Sequence Files/FASTA/' + plate + '/' + clone[0] + '/'
	fasta_path = fasta_dest + file_name + '.fasta'
	
	return author, run_date, abi_path, fasta_path

def blast():
	pass

# Initiates the application interface.
if __name__ == "__main__":
	man = manual_input(None)
	man.title('Manual input')
	man.mainloop()
