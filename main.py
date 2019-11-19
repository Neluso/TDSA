from characterization import characterization
from imaging import imaging
from tkinter import *
from tkinter import messagebox


def characterize():
    try:
        thick = float(thickness.get())
    except:
        messagebox.showerror('Error', 'Thickness mus be a real number')
        return 0
    characterization(bool(show_plots.get()), thick)


def image():
    imaging(bool(show_plots.get()))


char_text = 'Characterization: browses for the data for reference and sample and performs the refractive index '
char_text += 'extraction'
imag_text = 'Imaging: uses the data in "imaging_data" directory to perform an image reconstruction'


master_window = Tk()
master_window.title('Time Domain Spectroscopy Analyzer (TDSA)')
master_window.geometry('450x100')


show_plots = BooleanVar()
thickness = StringVar()


btn_char = Button(master_window, text='Characterization', command=characterize)
btn_char.grid(column=0, row=0)
btn_imag = Button(master_window, text='Imaging', command=image)
btn_imag.grid(column=1, row=0)
btn_quit = Button(master_window, text='Quit', command=quit)
btn_quit.grid(column=2, row=0, sticky='w')

thick_label = Label(master_window, text='Sample thickness (mm)')
thick_label.grid(column=0, row=1, sticky='w', columnspan=2)
thick_entry = Entry(master_window, textvariable=thickness)
thick_entry.grid(column=2, row=1)

rad_lbl = Label(master_window, text='Show plots?')
rad_lbl.grid(column=0, row=2, sticky='w')
rad_yes = Radiobutton(master_window, text='Yes', variable=show_plots, value=True, command='')
rad_yes.grid(column=1, row=2)
rad_no = Radiobutton(master_window, text='No', variable=show_plots, value=False, command='')
rad_no.grid(column=1, row=3)
# char_label = Label(master_window, text=char_text)
# char_label.grid(row=1, columnspan=5000, sticky='w')
# imag_label = Label(master_window, text=imag_text)
# imag_label.grid(row=2, columnspan=5000, sticky='w')


master_window.mainloop()
