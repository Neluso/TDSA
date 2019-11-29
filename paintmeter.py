from tkinter import *
from tkinter.filedialog import askopenfilenames
from tkinter import messagebox
from read_data import read_data
from numpy import *
from DPS_functions import *
from TDS_constants import *
from aux_functions import *


def paint_meter(show_plots):
    # messagebox.showinfo('Building', 'This function is still being built')
    Tk().withdraw()
    ref_file_path = askopenfilenames(initialdir='./data', title='Select reference file/s')
    t_refs, E_refs, is_error = read_data(ref_file_path)
    if is_error:
        return 0
    if len(t_refs.shape) > 1:
        t_ref = mean(t_refs, 0)
        E_ref = mean(E_refs, 0)
    else:
        t_ref = t_refs[0]
        E_ref = E_refs[0]
    t_ref_aux = t_ref
    E_ref_aux = E_ref
    t_junk = t_ref
    E_junk = E_ref
    sam_file_path = askopenfilenames(initialdir=ref_file_path, title='Select sample file/s')
    t_sams, E_sams, is_error = read_data(sam_file_path)
    if is_error:
        return 0
    t_ref *= 1e-12  # seconds
    t_sams *= 1e-12  # seconds

    nSamp = E_ref.size
    nSamp_pow = nextpow2(nSamp)
    f_ref, E_ref_w = fourier_analysis(t_ref, E_ref, nSamp_pow)
    f_min_idx, f_max_idx = f_min_max_idx(f_ref)
    E_ref_w = abs(E_ref_w)[f_min_idx:f_max_idx*2]
    f_ref = f_ref[f_min_idx:f_max_idx*2]

    E_ref_w_filt = rect_low_filter(E_ref_w, 0.5)
    plot(f_ref, E_ref_w)
    plot(f_ref, E_ref_w_filt)
    show()
    
    return 0
