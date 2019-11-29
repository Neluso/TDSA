from DPS_functions import *
from aux_functions import *
from read_data import read_data
from tkinter import *
from tkinter.filedialog import askopenfilenames
from matplotlib.pyplot import *


def spectra_plot():
    Tk().withdraw()
    filenames = askopenfilenames(initialdir='./data', title='Select reference data')
    t_s, E_s, is_error = read_data(filenames)
    for i in range(E_s.shape[0]):
        t_s_i = t_s[i] * 1e-12
        E_s_i = E_s[i]
        nSamp = E_s_i.size
        nSampPow = nextpow2(nSamp)
        f_i, E_w_i = fourier_analysis(t_s_i, E_s_i, nSampPow)
        f_min_idx, f_max_idx = f_min_max_idx(f_i)
        plot_length = 3
        E_w_i = toDb(E_w_i)
        t_s_i *= 1e12
        f_i *= 1e-12
        
        figure(1)
        xlabel(r'$\Delta t \  (ps)$')
        xlim([0, t_s_i[-1]])
        plot(t_s_i, E_s_i, lw=1)
        figure(2)
        plot(f_i[f_min_idx:plot_length * f_max_idx], E_w_i[f_min_idx:plot_length * f_max_idx], lw=1)
        xlabel(r'$f \ (THz)$')
        ylabel(r'$E_{\omega} \ (dB)$')
        xlim([f_i[f_min_idx], f_i[plot_length * f_max_idx]])
    show()
    return 0