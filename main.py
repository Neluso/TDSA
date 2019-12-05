from characterization import characterization
from imaging import imaging
from multilayer_meter import paint_meter
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


btn_char = Button(master_window, text='Characterization', command=characterize)
btn_char.place(width=100)
btn_imag = Button(master_window, text='Imaging', command=image)
btn_imag.place(width=70, x=100)
btn_imag = Button(master_window, text='Paintmeter', command=paint_measure)
btn_imag.place(width=90, x=170)
btn_imag = Button(master_window, text='Spectra', command=spectrum)
btn_imag.place(width=70, x=260)
btn_quit = Button(master_window, text='Quit', command=quit)
btn_quit.place(width=70, x=330)


# General Tweaks block
block_x = 0
block_y = 50
gen_label = Label(master_window, text='General tweaks', font='Helvetica 12 bold')
gen_label.place(x=block_x, y=block_y, anchor='w')
# Show plots? sub-block
sub_block_x = block_x
sub_block_y = block_y+25
rad_lbl = Label(master_window, text='Show plots?')
rad_lbl.place(x=sub_block_x, y=sub_block_y, anchor='w')
rad_yes = Radiobutton(master_window, text='Yes', variable=show_plots, value=True, command='')
rad_yes.place(x=sub_block_x, y=sub_block_y+25, anchor='w')
rad_yes.select()
rad_no = Radiobutton(master_window, text='No', variable=show_plots, value=False, command='')
rad_no.place(x=sub_block_x, y=sub_block_y+50, anchor='w')
rad_no.deselect()
# Temporal window sub-block
sub_block_x = block_x+100
sub_block_y = block_y+25
rad_lbl2 = Label(master_window, text='Temporal window: ')
rad_lbl2.place(x=sub_block_x, y=sub_block_y, anchor='w')
sub_block_x -= 125
sub_block_y += 25
rad_tuk = Radiobutton(master_window, text='Tukey', variable=temporal_window, value='tukey', command='')
rad_tuk.place(x=sub_block_x+125, y=sub_block_y, anchor='w')
rad_tuk.select()
rad_none = Radiobutton(master_window, text='None', variable=temporal_window, value='None', command='')
rad_none.place(x=sub_block_x+125, y=sub_block_y+25, anchor='w')
rad_none.deselect()
rad_fexp = Radiobutton(master_window, text='Force-Exp', variable=temporal_window, value='force_exp', command='')
rad_fexp.place(x=sub_block_x+195, y=sub_block_y, anchor='w')
rad_fexp.deselect()
rad_hann = Radiobutton(master_window, text='Hann', variable=temporal_window, value='hann', command='')
rad_hann.place(x=sub_block_x+195, y=sub_block_y+25, anchor='w')
rad_hann.deselect()
rad_bh = Radiobutton(master_window, text='Blackman-Harris', variable=temporal_window, value='blackman_harris', command='')
rad_bh.place(x=sub_block_x+275, y=sub_block_y+25, anchor='w')
rad_bh.deselect()
rad_cheb = Radiobutton(master_window, text='Chebishev', variable=temporal_window, value='chebisehev', command='')
rad_cheb.place(x=sub_block_x+275, y=sub_block_y, anchor='w')
rad_cheb.deselect()

# Characterization tweaks block
block_x = 0
block_y = 150
char_label = Label(master_window, text='Characterization tweaks', font='Helvetica 12 bold')
char_label.place(x=block_x, y=block_y)
thick_label = Label(master_window, text=' Sample thickness (mm)')
thick_label.place(x=block_x, y=block_y+25)
thick_entry = Entry(master_window, textvariable=thickness)
thick_entry.place(x=block_x+200, y=block_y+25)
noise_label = Label(master_window, text=' Noise floor (THz)')
noise_label.place(x=block_x, y=block_y+50)
noise_entry = Entry(master_window, textvariable=noise_floor_freq)
noise_entry.place(x=block_x+200, y=block_y+50)

# Imaging tweaks block
block_x = 0
block_y = 250
imag_label = Label(master_window, text='Imaging tweaks', font='Helvetica 12 bold')
imag_label.place(x=block_x, y=block_y)
hoff_label = Label(master_window, text='Horizontal offset (mm)')
hoff_label.place(x=block_x, y=block_y+25)
hoff_entry = Entry(master_window, textvariable=hoffset)
hoff_entry.place(x=block_x+200, y=block_y+25)
voff_label = Label(master_window, text='Vertical offset (mm)')
voff_label.place(x=block_x, y=block_y+50)
voff_entry = Entry(master_window, textvariable=voffset)
voff_entry.place(x=block_x+200, y=block_y+50)
resol_label = Label(master_window, text='Resolution (mm)')
resol_label.place(x=block_x, y=block_y+75)
resol_entry = Entry(master_window, textvariable=resolution)
resol_entry.place(x=block_x+200, y=block_y+75)


dims = '400x375'
master_window.geometry(dims)

master_window.mainloop()
