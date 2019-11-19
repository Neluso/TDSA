from characterization import characterization
from imaging import imaging
from tkinter import *


char_text = 'Characterization: browses for the data for reference and sample and performs the refractive index '
char_text += 'extraction'
imag_text = 'Imaging: uses the data in "imaging_data" directory to perform an image reconstruction'


def characterize():
    characterization(master_window, log_label)


def image():
    imaging(master_window)


master_window = Tk()
master_window.title('Time Domain Spectroscopy Analyzer (TDSA)')
master_window.geometry('750x200')
btn_char = Button(master_window, text='Characterization', command=characterize)
btn_char.grid(column=0, row=0)
btn_imag = Button(master_window, text='Imaging', command=image)
btn_imag.grid(column=1, row=0)
btn_quit = Button(master_window, text='Quit', command=quit)
btn_quit.grid(column=2, row=0)
log_label = Label(master_window, text=char_text)
log_label.grid(row=1, columnspan=5000, sticky='w')
log_label = Label(master_window, text=imag_text)
log_label.grid(row=2, columnspan=5000, sticky='w')


master_window.mainloop()
