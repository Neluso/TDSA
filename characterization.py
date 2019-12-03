# This script uses  the time domain mapping from a Reference pulse and the Sample pulse to perform spectroscopic
# analysis of the sample.


from tkinter import *
from tkinter.filedialog import askopenfilenames
from tkinter import messagebox
from read_data import read_data
from jepsen_index import jepsen_index
from DSP_functions import *
from aux_functions import *
from numpy import *
from matplotlib.pyplot import *
from TDS_constants import *


def characterization(show_plots, thickness, temp_window, noise_floor_freq):
    
    Tk().withdraw()
    ref_file_path = askopenfilenames(initialdir='./data', title='Select reference data')
    t_ref, E_ref, is_error = read_data(ref_file_path)
    if len(t_ref.shape) > 1:
        t_ref = mean(t_ref, 0)
        E_ref = mean(E_ref, 0)
    else:
        t_ref = t_ref[0]
        E_ref = E_ref[0]
    t_ref_aux = t_ref
    E_ref_aux = E_ref
    t_junk = t_ref
    E_junk = E_ref
    if is_error:
        return 0
    Tk().withdraw()
    sam_file_path = askopenfilenames(initialdir=ref_file_path, title='Select sample data')
    t_sams, E_sams, is_error = read_data(sam_file_path)
    if is_error:
        return 0

    t_ref *= 1e-12  # seconds
    t_sams *= 1e-12  # seconds
    thickness *= 1e-3  # 1.95e-3  # m
    noise_floor_freq *= 1e12  # to THz
    
    if temp_window == 'blackman_harris':
        t_ref, E_ref, t_junk, E_junk = bh_windowing(t_ref_aux, E_ref_aux, t_junk, E_junk)
    if temp_window == 'chebisehev':
        t_ref, E_ref, t_junk, E_junk = cheb_windowing(t_ref_aux, E_ref_aux, t_junk, E_junk)
    if temp_window == 'hann':
        t_ref, E_ref, t_junk, E_junk = hann_windowing(t_ref_aux, E_ref_aux, t_junk, E_junk)
    if temp_window == 'force_exp':
        t_ref, E_ref, t_junk, E_junk = force_exp_windowing(t_ref_aux, E_ref_aux, t_junk, E_junk)
    if temp_window == 'tukey':
        t_ref, E_ref, t_junk, E_junk = tukey_windowing(t_ref_aux, E_ref_aux, t_junk, E_junk)

    
    nSamp = E_ref.size
    nSamp_pow = nextpow2(nSamp)
    f_ref, E_ref_w = fourier_analysis(t_ref, E_ref, nSamp_pow)
    
    n = list()
    alphas = list()
    t_sam_mean = list()
    E_sam_mean = list()
    f_sam_mean = list()
    E_sam_w_mean = list()
    
    for i in range(t_sams.shape[0]):
        
        t_sam = t_sams[i]
        E_sam = E_sams[i]
        
        if temp_window == 'blackman_harris':
            t_junk, E_junk, t_sam, E_sam = bh_windowing(t_ref_aux, E_ref_aux, t_sam, E_sam)
        if temp_window == 'chebisehev':
            t_junk, E_junk, t_sam, E_sam = cheb_windowing(t_ref_aux, E_ref_aux, t_sam, E_sam)
        if temp_window == 'hann':
            t_junk, E_junk, t_sam, E_sam = hann_windowing(t_ref_aux, E_ref_aux, t_sam, E_sam)
        if temp_window == 'force_exp':
            t_junk, E_junk, t_sam, E_sam = force_exp_windowing(t_ref_aux, E_ref_aux, t_sam, E_sam)
        if temp_window == 'tukey':
            t_junk, E_junk, t_sam, E_sam = tukey_windowing(t_ref_aux, E_ref_aux, t_sam, E_sam)
        
        n_i, alpha_f_i = jepsen_index(t_ref, E_ref, t_sam, E_sam, thickness)
        f_sam, E_sam_w = fourier_analysis(t_sam, E_sam, nSamp_pow)

        n.append(n_i)
        alphas.append(alpha_f_i)
        t_sam_mean.append(t_sam)
        E_sam_mean.append(E_sam)
        f_sam_mean.append(f_sam)
        E_sam_w_mean.append(E_sam_w)
    
    t_sam_mean = array(t_sam_mean)
    E_sam_mean = array(E_sam_mean)
    f_sam_mean = array(f_sam_mean)
    E_sam_w_mean = array(E_sam_w_mean)
    n = array(n)
    alphas = array(alphas)
    t_sam = mean(t_sam_mean, 0)
    E_sam = mean(E_sam_mean, 0)
    f_sam = mean(f_sam_mean, 0)
    E_sam_w = mean(E_sam_w_mean, 0)
    sigma_n = std(n, 0)
    n = mean(n, 0)
    alpha_f = mean(alphas, 0)
    sigma_alpha_f = std(alphas, 0)

    f_min_idx, f_max_idx = f_min_max_idx(f_ref)
    noisefloor = noise_floor(f_ref, prettyfy(E_ref_w, amax(E_ref_w)), noise_floor_freq)
    
    t_ref *= 1e12
    t_sam *= 1e12
    f_ref *= 1e-12
    f_sam *= 1e-12
    
    figure(1)  # THz time domain pulses
    fig_name = 'Pulses'
    plot(t_ref, E_ref, label='Reference', linewidth=1)
    plot(t_sam, E_sam, label='Sample', linewidth=1)
    xlabel(r'$\Delta t \  (ps)$')
    xlim([0, t_ref[-1]])
    title(fig_name)
    legend()
    savefig('./output/' + fig_name + '_' + '.svg', format='svg')

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
    savefig('./output/' + fig_name + '_' + '.svg', format='svg')

    figure(3)  # Alpha_f
    fig_name = 'Absorption'
    plot_length = 1
    # plot(f_ref[f_min_idx:plot_length * f_max_idx], 0.01 * real(alpha_f[f_min_idx:plot_length * f_max_idx]), linewidth=1)
    errorbar(f_ref[f_min_idx:plot_length * f_max_idx], 0.01 * real(alpha_f[f_min_idx:plot_length * f_max_idx]),
             yerr=0.01 * real(sigma_alpha_f[f_min_idx:plot_length * f_max_idx]), lw=1)
    xlabel(r'$f \ (THz)$')
    ylabel(r'$\alpha \ (cm^{-1})$')
    xlim([f_ref[f_min_idx], f_ref[plot_length * f_max_idx]])
    title(fig_name + ' coefficient')
    savefig('./output/' + fig_name + '_' + '.svg', format='svg')
    
    figure(4)  # Refractive Index
    fig_name = 'Index'
    plot_length = 2
    # plot(f_ref[f_min_idx:plot_length * f_max_idx], n[f_min_idx:plot_length * f_max_idx], linewidth=1)
    errorbar(f_ref[f_min_idx:plot_length * f_max_idx], n[f_min_idx:plot_length * f_max_idx],
             yerr=sigma_n[f_min_idx:plot_length * f_max_idx], lw=1)
    xlabel(r'$f \ (THz)$')
    xlim([f_ref[f_min_idx], f_ref[plot_length * f_max_idx]])
    title(fig_name)
    savefig('./output/' + fig_name + '_' + '.svg', format='svg')
    
    if show_plots:
        show()
    
    close('all')
    messagebox.showinfo('Process ended correctly', 'Output saved in directory')
    
    t0_ref = t_ref[centre_loc(E_ref)]
    t0_sam = t_sam[centre_loc(E_sam)]
    n_mean = 1 + (t0_sam - t0_ref) * c_0 / thickness
    
    return 0
