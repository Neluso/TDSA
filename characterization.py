# This script uses  the time domain mapping from a Reference pulse and the Sample pulse to perform spectroscopic
# analysis of the sample.


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


def characterization(show_plots, thickness, temp_window):  # , temp_window):
    
    # Tk().withdraw()
    # ref_file_path = askopenfilename(initialdir='./data', title='Select reference data')
    ref_file_path = './data/demo_data/test_ref.txt'
    t_ref, E_ref, is_error = read_data(ref_file_path)
    if is_error:
        return 0
    # Tk().withdraw()
    # sam_file_path = askopenfilename(initialdir=ref_file_path, title='Select sample data')
    sam_file_path = './data/demo_data/test_sam.txt'
    t_sam, E_sam, is_error = read_data(sam_file_path)
    if is_error:
        return 0

    t_ref *= 1e-12  # seconds
    t_sam *= 1e-12  # seconds
    thickness *= 1e-3  # 1.95e-3  # m
    
    if True:  #temp_window is 'blackman_harris':
        t_ref, E_ref, t_sam, E_sam= force_exp_windowing(t_ref, E_ref, t_sam, E_sam)
        # t_sam, E_sam = bh_windowing(t_sam, E_sam)
    
    
    nSamp = E_ref.size
    nSamp_pow = nextpow2(nSamp)


    n, alpha_f = jepsen_index(t_ref, E_ref, t_sam, E_sam, thickness)
    f_ref, E_ref_w = fourier_analysis(t_ref, E_ref, nSamp_pow)
    f_sam, E_sam_w = fourier_analysis(t_sam, E_sam, nSamp_pow)
    f_min_idx, f_max_idx = f_min_max_idx(f_ref)
    noisefloor = noise_floor(f_ref, prettyfy(E_ref_w, amax(E_ref_w)), 4e12)
    
    sam_file = sam_file_path.split('/')[-1]
    sam_file = sam_file.replace('.', '_')
    sam_file = sam_file.replace(' ', '_')
    save_path = write_spectra(f_ref, E_ref_w, f_sam, E_sam_w, sam_file)
    write_data(f_ref, n, 'refractive_index', save_path)
    copy(ref_file_path, save_path)
    copy(sam_file_path, save_path)

    figure(1)  # THz time domain pulses
    fig_name = 'Pulses'
    plot(t_ref, E_ref, label='Reference', linewidth=1)
    plot(t_sam, E_sam, label='Sample', linewidth=1)
    xlabel(r'$\Delta t \  (ps)$')
    xlim([t_ref[0], t_ref[-1]])
    title(fig_name)
    legend()
    savefig(save_path + fig_name + '_' + sam_file + '.svg', format='svg')

    figure(2)  # Spectra
    fig_name = "Spectra"
    plot_length = 6
    plot(f_ref[f_min_idx:plot_length * f_max_idx],
         prettyfy(E_ref_w, amax(E_ref_w))[f_min_idx:plot_length * f_max_idx], label='Reference', linewidth=1)
    plot(f_sam[f_min_idx:plot_length * f_max_idx],
         prettyfy(E_sam_w, amax(E_ref_w))[f_min_idx:plot_length * f_max_idx], label='Sample', linewidth=1)
    plot(f_ref[f_min_idx:plot_length * f_max_idx], noisefloor[f_min_idx:plot_length * f_max_idx], 'r--',
         label='Noise Floor', linewidth=1)
    xlabel(r'$f \ (THz)$')
    ylabel(r'$E_{\omega} \ (dB)$')
    xlim([f_ref[f_min_idx], f_ref[plot_length * f_max_idx]])
    title(fig_name)
    legend()
    savefig(save_path + fig_name + '_' + sam_file + '.svg', format='svg')

    figure(3)  # Alpha_f
    fig_name = 'Absorption'
    plot_length = 2
    plot(f_ref[f_min_idx:plot_length * f_max_idx], 0.01 * real(alpha_f[f_min_idx:plot_length * f_max_idx]), linewidth=1)
    xlabel(r'$f \ (THz)$')
    ylabel(r'$\alpha \ (cm^{-1})$')
    xlim([f_ref[f_min_idx], f_ref[plot_length * f_max_idx]])
    title(fig_name + ' coefficient')
    savefig(save_path + fig_name + '_' + sam_file + '.svg', format='svg')

    figure(4)  # Refractive Index
    fig_name = 'Index'
    plot_length = 2
    plot(f_ref[f_min_idx:plot_length * f_max_idx], n[f_min_idx:plot_length * f_max_idx], linewidth=1)
    xlabel(r'$f \ (THz)$')
    xlim([f_ref[f_min_idx], f_ref[plot_length * f_max_idx]])
    title(fig_name)
    savefig(save_path + fig_name + '_' + sam_file + '.svg', format='svg')
    
    if show_plots:
        show()
    
    close('all')
    # messagebox.showinfo('Process ended correctly', 'Output saved in ' + save_path) TODO uncomment
    
    return 0
