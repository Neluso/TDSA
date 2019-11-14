from tkinter import Tk
from tkinter.filedialog import askopenfilename
from matplotlib.pyplot import plot, show, legend, xlabel, ylabel, figure, title
from read_data import read_data
from write_data import write_data
from jepsen_index import jepsen_index
from DPS_functions import *
from aux_functions import *

Tk().withdraw()
ref_file = askopenfilename(initialdir='./data', title='Select reference data')
Tk().withdraw()
sam_file = askopenfilename(initialdir='./data', title='Select sample data')


t_ref, E_ref = read_data(ref_file)
t_sam, E_sam = read_data(sam_file)
t_ref *= 1e-12  # seconds
t_sam *= 1e-12  # seconds
thickness = 1.95e-3  # m


nSamp = E_ref.size
nSamp_pow = nextpow2(nSamp)


n, alpha_f = jepsen_index(t_ref, E_ref, t_sam, E_sam, thickness)
f_ref, E_ref_w = fourier_analysis(t_ref, E_ref, nSamp)
f_sam, E_sam_w = fourier_analysis(t_sam, E_sam, nSamp)
f_min_idx, f_max_idx = f_min_max_idx(f_ref)


print('Plotting data')
figure(1)
plot(f_ref[f_min_idx:2*f_max_idx], 0.01 * alpha_f[f_min_idx:2*f_max_idx])
title('Absorption')
xlabel('f (THz)')
ylabel(r'$\alpha (cm^{-1})$')
legend()


show()


print('Writting data')
write_data(f_ref, E_ref_w, E_sam_w)
print('Closing')