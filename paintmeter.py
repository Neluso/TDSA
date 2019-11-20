from tkinter import *
from tkinter.filedialog import askopenfilename
from tkinter import messagebox
from read_data import read_data


def paint_meter(show_plots):
    messagebox.showinfo('Building', 'This function is still being built')
    Tk().withdraw()
    ref_file_path = askopenfilename(initialdir='./data', title='Select refractive index file')
    f_ref_rind, n_ref_ind, is_error = read_data(ref_file_path)
    if is_error:
        return 0
    
    return 0
