from tkinter import *
from tkinter.filedialog import askopenfilenames
from tkinter import messagebox
from matplotlib.pyplot import *
from read_data import read_data
from write_data import *
from jepsen_index import jepsen_index
from DPS_functions import *
from aux_functions import *
from numpy import *
from shutil import *


Tk().withdraw()
ref_file_path = askopenfilenames(initialdir='./data', title='Select reference data')
t_ref, E_ref, is_error = read_data(ref_file_path)
if is_error:
    quit()
Tk().withdraw()
sam_file_path = askopenfilenames(initialdir=ref_file_path, title='Select sample data')
t_sam, E_sam, is_error = read_data(sam_file_path)
if is_error:
    quit()
Tk().withdraw()
edge_file_path = askopenfilenames(initialdir=ref_file_path, title='Select sample data')
t_edge, E_edge, is_error = read_data(edge_file_path)
if is_error:
    quit()

t_ref = mean(t_ref, 0)
E_ref = mean(E_ref, 0)
t_sam = mean(t_sam, 0)
E_sam = mean(E_sam, 0)
t_edge = mean(t_edge, 0)
E_edge = mean(E_edge, 0)

t_ref *= 1e-12  # seconds
t_sam *= 1e-12  # seconds
t_edge *= 1e-12  # seconds

t_ref_aux = t_ref
E_ref_aux = E_ref

t_ref, E_ref, t_sam, E_sam = bh_windowing(t_ref, E_ref, t_sam, E_sam)
t_ref_aux, E_ref_aux, t_edge, E_edge = bh_windowing(t_ref_aux, E_ref_aux, t_edge, E_edge)

nSamp = E_ref.size
nSamp_pow = nextpow2(nSamp)

f_ref, E_ref_w = fourier_analysis(t_ref, E_ref, nSamp_pow)
f_sam, E_sam_w = fourier_analysis(t_sam, E_sam, nSamp_pow)
f_edge, E_edge_w = fourier_analysis(t_edge, E_edge, nSamp_pow)


E_ref_w = abs(E_ref_w)
E_sam_w = abs(E_sam_w)
E_edge_w = abs(E_edge_w)

f_min_idx, f_max_idx = f_min_max_idx(f_ref)

n_sam, alpha_sam = jepsen_index(t_ref, E_ref, t_sam, E_sam, 1.0)

# figure(1)
# plot(t_ref, E_ref, label='Ref', lw=1)
# plot(t_sam, E_sam, label='Sam', lw=1)
# plot(t_edge, E_edge, label='Edge', lw=1)
# legend()

E_ref_w_max = amax(E_ref_w)
plot_length = 3
figure(2)
plot(f_ref[f_min_idx:plot_length * f_max_idx], prettyfy(E_ref_w, E_ref_w_max)[f_min_idx:plot_length * f_max_idx], label='Ref', lw=1)
plot(f_sam[f_min_idx:plot_length * f_max_idx], prettyfy(E_sam_w, E_ref_w_max)[f_min_idx:plot_length * f_max_idx], label='Sam', lw=1)
plot(f_edge[f_min_idx:plot_length * f_max_idx], prettyfy(E_edge_w, E_ref_w_max)[f_min_idx:plot_length * f_max_idx], label='Edge', lw=1)
legend()


show()

