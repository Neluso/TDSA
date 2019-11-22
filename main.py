from characterization import characterization
from imaging import imaging
from paintmeter import paint_meter
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
    try:
        hoff = float(hoffset.get())
        voff = float(voffset.get())
        resol = float(resolution.get())
    except:
        messagebox.showerror('Error', 'Offset values must be real numbers')
        return 0
    imaging(bool(show_plots.get()), hoff, voff, resol, master_window)


def paint_measure():
    paint_meter(bool(show_plots.get()))


char_text = 'Characterization: browses for the data for reference and sample and performs the refractive index '
char_text += 'extraction'
imag_text = 'Imaging: uses the data in "imaging_data" directory to perform an image reconstruction'


master_window = Tk()
master_window.title('Time Domain Spectroscopy Analyzer (TDSA)')
master_window.geometry('450x300')


show_plots = BooleanVar()
thickness = StringVar(master_window, value='1.95')  # mm
hoffset = StringVar(master_window, value='0')  # mm
voffset = StringVar(master_window, value='20')  # mm
resolution = StringVar(master_window, value='0.5')  # delta_mm


btn_char = Button(master_window, text='Characterization', command=characterize)
btn_char.grid(column=0, row=0)
btn_imag = Button(master_window, text='Imaging', command=image)
btn_imag.grid(column=1, row=0)
btn_imag = Button(master_window, text='Paintmeter', command=paint_measure)
btn_imag.grid(column=2, row=0)
btn_quit = Button(master_window, text='Quit', command=quit)
btn_quit.grid(column=3, row=0, sticky='w')

gen_label = Label(master_window, text='General tweaks', font='Helvetica 12 bold')
gen_label.grid(column=0, row=1, sticky='w', columnspan=20)
rad_lbl = Label(master_window, text='Show plots?')
rad_lbl.grid(column=0, row=2, sticky='w')
rad_yes = Radiobutton(master_window, text='Yes', variable=show_plots, value=True, command='')
rad_yes.grid(column=1, row=2)
rad_yes.select()
rad_no = Radiobutton(master_window, text='No', variable=show_plots, value=False, command='')
rad_no.grid(column=1, row=3)
rad_no.deselect()

spacer_1 = Label(master_window, text='')
spacer_1.grid(column=0, row=4, sticky='w')

char_label = Label(master_window, text='Characterization tweaks', font='Helvetica 12 bold')
char_label.grid(column=0, row=5, sticky='w', columnspan=20)
thick_label = Label(master_window, text='Sample thickness (mm)')
thick_label.grid(column=0, row=6, sticky='w', columnspan=20)
thick_entry = Entry(master_window, textvariable=thickness)
thick_entry.grid(column=2, row=6, columnspan=20)

spacer_1 = Label(master_window, text='')
spacer_1.grid(column=0, row=7, sticky='w')

imag_label = Label(master_window, text='Imaging tweaks', font='Helvetica 12 bold')
imag_label.grid(column=0, row=8, sticky='w', columnspan=20)
hoff_label = Label(master_window, text='Horizontal offset (mm)')
hoff_label.grid(column=0, row=9, sticky='w', columnspan=20)
hoff_entry = Entry(master_window, textvariable=hoffset)
hoff_entry.grid(column=2, row=9, columnspan=20)
voff_label = Label(master_window, text='Vertical offset (mm)')
voff_label.grid(column=0, row=10, sticky='w', columnspan=20)
voff_entry = Entry(master_window, textvariable=voffset)
voff_entry.grid(column=2, row=10, columnspan=20)
resol_label = Label(master_window, text='Resolution (mm)')
resol_label.grid(column=0, row=11, sticky='w', columnspan=20)
resol_entry = Entry(master_window, textvariable=resolution)
resol_entry.grid(column=2, row=11, columnspan=20)


master_window.mainloop()
