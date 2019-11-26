from tkinter import *
from tkinter.filedialog import askopenfilename
from tkinter import messagebox
from matplotlib.pyplot import plot, show, legend, xlabel, ylabel, figure, title, savefig, close, xlim, axvline
from read_data import read_data
from write_data import *
from jepsen_index import jepsen_index
from DPS_functions import *
from aux_functions import *
from numpy import amax, real
from shutil import copy


Tk().withdraw()
ref_file_path = askopenfilename(initialdir='./data', title='Select reference data')
# ref_file_path = './data/demo_data/test_ref.txt'
t_ref, E_ref, is_error = read_data(ref_file_path)
if is_error:
    quit()
Tk().withdraw()
sam_file_path = askopenfilename(initialdir=ref_file_path, title='Select sample data')
# sam_file_path = './data/demo_data/test_sam.txt'
t_sam, E_sam, is_error = read_data(sam_file_path)
if is_error:
    quit()
t_ref *= 1e-12  # seconds
t_sam *= 1e-12  # seconds

t_ref_bh, E_ref_bh, t_aux_j, E_aux_j = bh_windowing(t_ref, E_ref, t_sam, E_sam)
t_ref_ch, E_ref_ch, t_aux_j, E_aux_j = cheb_windowing(t_ref, E_ref, t_sam, E_sam)
t_ref_hn, E_ref_hn, t_aux_j, E_aux_j = hann_windowing(t_ref, E_ref, t_sam, E_sam)
t_ref_fx, E_ref_fx, t_aux_j, E_aux_j = force_exp_windowing(t_ref, E_ref, t_sam, E_sam)
t_ref_tk, E_ref_tk, t_aux_j, E_aux_j = force_exp_windowing(t_ref, E_ref, t_sam, E_sam)

nSamp = E_ref.size
nSamp_pow = nextpow2(nSamp)

f_ref, E_ref_w = fourier_analysis(t_ref, E_ref, nSamp_pow)
f_ref_bh, E_ref_w_bh = fourier_analysis(t_ref_bh, E_ref_bh, nSamp_pow)
f_ref_ch, E_ref_w_ch = fourier_analysis(t_ref_ch, E_ref_ch, nSamp_pow)
f_ref_hn, E_ref_w_hn = fourier_analysis(t_ref_hn, E_ref_hn, nSamp_pow)
f_ref_fx, E_ref_w_fx = fourier_analysis(t_ref_fx, E_ref_fx, nSamp_pow)
f_ref_tk, E_ref_w_tk = fourier_analysis(t_ref_tk, E_ref_tk, nSamp_pow)


f_min_idx, f_max_idx = f_min_max_idx(f_ref)
E_max = 1  # amax(E_ref_w)
plot_length = 5
plot(f_ref[f_min_idx:plot_length * f_max_idx], prettyfy(E_ref_w, E_max)[f_min_idx:plot_length * f_max_idx], label='no window', lw=1)
plot(f_ref_bh[f_min_idx:plot_length * f_max_idx], prettyfy(E_ref_w_bh[f_min_idx:plot_length * f_max_idx], E_max), label='b-harr', lw=1)
plot(f_ref_ch[f_min_idx:plot_length * f_max_idx], prettyfy(E_ref_w_ch[f_min_idx:plot_length * f_max_idx], E_max), label='cheb', lw=1)
plot(f_ref_hn[f_min_idx:plot_length * f_max_idx], prettyfy(E_ref_w_hn[f_min_idx:plot_length * f_max_idx], E_max), label='hann', lw=1)
plot(f_ref_fx[f_min_idx:plot_length * f_max_idx], prettyfy(E_ref_w_fx[f_min_idx:plot_length * f_max_idx], E_max), label='fexp', lw=1)
plot(f_ref_tk[f_min_idx:plot_length * f_max_idx], prettyfy(E_ref_w_tk[f_min_idx:plot_length * f_max_idx], E_max), label='tuk', lw=1)
xlim([f_ref[f_min_idx], f_ref[plot_length * f_max_idx]])
legend()
show()
