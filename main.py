from characterization import characterization
from imaging import imaging
from paintmeter import paint_meter
from tkinter import *
from tkinter import messagebox
from spectra_plot import spectra_plot


def characterize():
    try:
        thick = float(thickness.get())
        n_f_f = float(noise_floor_freq.get())
    except:
        messagebox.showerror('Error', 'Thickness and noise floor frequency must be a real number')
        return 0
    characterization(bool(show_plots.get()), thick, temporal_window.get(), n_f_f)  #, bool(dispersive_media.get()))
    return 0


def image():
    try:
        hoff = float(hoffset.get())
        voff = float(voffset.get())
        resol = float(resolution.get())
    except:
        messagebox.showerror('Error', 'Offset values must be real numbers')
        return 0
    imaging(bool(show_plots.get()), hoff, voff, resol, temporal_window.get())
    return 0


def spectrum():
    spectra_plot()
    return 0


def paint_measure():
    paint_meter(bool(show_plots.get()))
    return 0


char_text = 'Characterization: browses for the data for reference and sample and performs the refractive index '
char_text += 'extraction'
imag_text = 'Imaging: uses the data in "imaging_data" directory to perform an image reconstruction'


master_window = Tk()
master_window.title('Time Domain Spectroscopy Analyzer (TDSA)')


show_plots = BooleanVar()
thickness = StringVar(master_window, value='0.2')  # mm
hoffset = StringVar(master_window, value='0')  # mm
voffset = StringVar(master_window, value='20')  # mm
resolution = StringVar(master_window, value='0.5')  # delta_mm
temporal_window = StringVar(master_window, value='tukey')
noise_floor_freq = StringVar(master_window, value='4')  # THz
dispersive_media = BooleanVar()

at_row = 0

btn_char = Button(master_window, text='Characterization', command=characterize)
btn_char.grid(column=0, row=at_row)
btn_imag = Button(master_window, text='Imaging', command=image)
btn_imag.grid(column=1, row=at_row)
btn_imag = Button(master_window, text='Paintmeter', command=paint_measure)
btn_imag.grid(column=2, row=at_row)
btn_imag = Button(master_window, text='Spectra', command=spectrum)
btn_imag.grid(column=3, row=at_row)
btn_quit = Button(master_window, text='Quit', command=quit)
btn_quit.grid(column=4, row=at_row, sticky='w')

at_row += 1  # row 1
gen_label = Label(master_window, text='General tweaks', font='Helvetica 12 bold')
gen_label.grid(column=0, row=at_row, sticky='w', columnspan=20)
at_row += 1  # row 2
rad_lbl = Label(master_window, text='Show plots?')
rad_lbl.grid(column=0, row=at_row, sticky='w')
rad_yes = Radiobutton(master_window, text='Yes', variable=show_plots, value=True, command='')
rad_yes.grid(column=1, row=at_row)
rad_yes.select()
at_row += 1  # row 3
rad_no = Radiobutton(master_window, text='No', variable=show_plots, value=False, command='')
rad_no.grid(column=1, row=at_row)
rad_no.deselect()

at_row -= 1  # row 2
rad_lbl2 = Label(master_window, text='Temporal window: ')
rad_lbl2.grid(column=3, row=at_row, sticky='w')
at_row += 1  #row 3
rad_tuk = Radiobutton(master_window, text='Tukey', variable=temporal_window, value='tukey', command='')
rad_tuk.grid(column=3, row=at_row, sticky='w')
rad_tuk.select()
rad_cheb = Radiobutton(master_window, text='Chebishev', variable=temporal_window, value='chebisehev', command='')
rad_cheb.grid(column=4, row=at_row, sticky='w')
rad_cheb.deselect()
at_row += 1  # row 4
rad_hann = Radiobutton(master_window, text='Hann', variable=temporal_window, value='hann', command='')
rad_hann.grid(column=3, row=at_row, sticky='w')
rad_hann.deselect()
# at_row += 1  # row 5
rad_fexp = Radiobutton(master_window, text='Force-Exp', variable=temporal_window, value='force_exp', command='')
rad_fexp.grid(column=4, row=at_row, sticky='w')
rad_fexp.deselect()
at_row += 1  # row 5
rad_bh = Radiobutton(master_window, text='Blackman-Harris', variable=temporal_window, value='blackman_harris', command='')
rad_bh.grid(column=4, row=at_row, sticky='w')
rad_bh.deselect()
# at_row += 1  # row 5
rad_none = Radiobutton(master_window, text='None', variable=temporal_window, value='None', command='')
rad_none.grid(column=3, row=at_row, sticky='w')
rad_none.deselect()

at_row += 1  # row 6
spacer_1 = Label(master_window, text='')
spacer_1.grid(column=0, row=at_row, sticky='w')

at_row += 1  # row 7
char_label = Label(master_window, text='Characterization tweaks', font='Helvetica 12 bold')
char_label.grid(column=0, row=at_row, sticky='w', columnspan=20)
at_row += 1  # row 8
thick_label = Label(master_window, text=' Sample thickness (mm)')
thick_label.grid(column=0, row=at_row, sticky='w', columnspan=20)
thick_entry = Entry(master_window, textvariable=thickness)
thick_entry.grid(column=2, row=at_row, columnspan=20)
at_row += 1  # row 2
rad_disp = Label(master_window, text='Dispersive media?')
rad_disp.grid(column=0, row=at_row, sticky='w')
disp_yes = Radiobutton(master_window, text='Yes', variable=dispersive_media, value=True, command='')
disp_yes.grid(column=1, row=at_row)
disp_yes.deselect()
disp_no = Radiobutton(master_window, text='No', variable=dispersive_media, value=False, command='')
disp_no.grid(column=2, row=at_row)
disp_no.select()
at_row += 1
noise_label = Label(master_window, text=' Noise floor (THz)')
noise_label.grid(column=0, row=at_row, sticky='w', columnspan=20)
noise_entry = Entry(master_window, textvariable=noise_floor_freq)
noise_entry.grid(column=2, row=at_row, columnspan=20)

at_row += 1  # row 9
spacer_1 = Label(master_window, text='')
spacer_1.grid(column=0, row=at_row, sticky='w')

at_row += 1  # row 10
imag_label = Label(master_window, text='Imaging tweaks', font='Helvetica 12 bold')
imag_label.grid(column=0, row=at_row, sticky='w', columnspan=20)
at_row += 1  # row 11
hoff_label = Label(master_window, text='Horizontal offset (mm)')
hoff_label.grid(column=0, row=at_row, sticky='w', columnspan=20)
hoff_entry = Entry(master_window, textvariable=hoffset)
hoff_entry.grid(column=2, row=at_row, columnspan=20)
at_row += 1  # row 12
voff_label = Label(master_window, text='Vertical offset (mm)')
voff_label.grid(column=0, row=at_row, sticky='w', columnspan=20)
voff_entry = Entry(master_window, textvariable=voffset)
voff_entry.grid(column=2, row=at_row, columnspan=20)
at_row += 1  # row 13
resol_label = Label(master_window, text='Resolution (mm)')
resol_label.grid(column=0, row=at_row, sticky='w', columnspan=20)
resol_entry = Entry(master_window, textvariable=resolution)
resol_entry.grid(column=2, row=at_row, columnspan=20)

dims = '500x' + str(at_row * 25 + 25)
master_window.geometry(dims)

master_window.mainloop()
