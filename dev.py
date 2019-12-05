# Developing script


from read_data import read_data
from numpy import *
from numpy.fft import *
from DSP_functions import *
from aux_functions import *
from matplotlib.pyplot import *
from tkinter import *
from tkinter.filedialog import askopenfilenames
from scipy import signal
from jepsen_index import jepsen_index
from tqdm import *


def bh_windowing_dev(t_val, E_val, t_sub, E_sub):
    center_val_idx = centre_loc(E_val)
    center_sub_idx = centre_loc(E_sub)
    left_padd_idxs = E_val.size - 2 * center_val_idx
    left_padd_idxs_sub = E_sub.size - 2 * center_sub_idx
    E_val = zero_padding(E_val, left_padd_idxs, 0)
    E_sub = zero_padding(E_sub, left_padd_idxs_sub, 0)
    t_val_rev = - flip(t_val[1:left_padd_idxs + 1])
    t_val = concatenate((t_val_rev, t_val))
    t_sub = concatenate((t_val_rev, t_sub))
    E_val *= signal.blackmanharris(E_val.size)
    win = signal.blackmanharris(E_sub.size)
    E_sub *= win
    E_sub = zero_padding(E_sub, left_padd_idxs - left_padd_idxs_sub, 0)
    win = zero_padding(win, left_padd_idxs - left_padd_idxs_sub, 0)
    return t_val, E_val, t_sub, E_sub, win


fh = open('./data/marca_autodestructiva/ref.txt')
data = fh.read()
data = data.split('\n')
t_ref = list()
E_ref = list()
for item in data:
    item = item.split(',')
    if item == ['']:
        break
    t_ref.append(float(item[0]))
    E_ref.append(float(item[1]))
t_ref = array(t_ref) * 1e-12
E_ref = array(E_ref)

fh = open('./data/marca_autodestructiva/sam.txt')
data = fh.read()
data = data.split('\n')
t_sam = list()
E_sam = list()
for item in data:
    item = item.split(',')
    if item == ['']:
        break
    t_sam.append(float(item[0]))
    E_sam.append(float(item[1]))
t_sam = array(t_sam) * 1e-12
E_sam = array(E_sam)


nSamp = E_ref.size
nSamp_pow = nextpow2(nSamp)


E_ref_centre = centre_loc(E_ref)
E_sam_centre = centre_loc(E_sam)
delta_time = 10  # 10 ps (peak span) / 2 = 5 ps
delta_time *= 1e-12
delta_time = t_ref[E_ref_centre]
delta_peak_ref_min = where(t_ref >= t_ref[E_ref_centre]-delta_time)[0][0]
delta_peak_ref_max = where(t_ref >= t_ref[E_ref_centre]+delta_time)[0][0]
delta_peak_sam_min = where(t_sam >= t_sam[E_sam_centre]-delta_time)[0][0]
delta_peak_sam_max = where(t_sam >= t_sam[E_sam_centre]+delta_time)[0][0]


win_ref = signal.tukey(delta_peak_ref_max - delta_peak_ref_min)
win_ref = zero_padding(win_ref, delta_peak_ref_min, E_ref.size - delta_peak_ref_max)
win_sam = signal.tukey(delta_peak_sam_max - delta_peak_sam_min)
win_sam = zero_padding(win_sam, delta_peak_sam_min, E_sam.size - delta_peak_sam_max)


E_ref_win = E_ref * win_ref
E_sam_win = E_sam * win_sam

f_ref, E_ref_w = fourier_analysis(t_ref, E_ref, nSamp_pow)
f_ref_win, E_ref_w_win = fourier_analysis(t_ref, E_ref_win, nSamp_pow)
f_sam, E_sam_w = fourier_analysis(t_sam, E_sam, nSamp_pow)
f_sam_win, E_sam_w_win = fourier_analysis(t_sam, E_sam_win, nSamp_pow)

# E_ref_w = abs(E_ref_w)
# E_sam_w = abs(E_sam_w)
# E_ref_w_win = abs(E_ref_w_win)
# E_sam_w_win = abs(E_sam_w_win)
#
# n, alpha_f, = jepsen_index(t_ref, E_ref, t_sam, E_sam, 0.97e-3)
# n_win, alpha_f_win = jepsen_index(t_ref, E_ref_win, t_sam, E_sam_win, 0.97e-3)


f_min_idx, f_max_idx = f_min_max_idx(f_ref)

plot_length = 1  # from 0.1 to plot_length THz

# figure(1)
# plot(f_ref[f_min_idx:plot_length * f_max_idx], n[f_min_idx:plot_length * f_max_idx], '-', label='no_win')
# plot(f_ref[f_min_idx:plot_length * f_max_idx], n_win[f_min_idx:plot_length * f_max_idx], '-', label='win')
# title('n')
# legend()
#
# figure(2)
# plot(f_ref[f_min_idx:plot_length * f_max_idx], 0.01 * alpha_f[f_min_idx:plot_length * f_max_idx], label='no_win')
# plot(f_ref[f_min_idx:plot_length * f_max_idx], 0.01 * alpha_f_win[f_min_idx:plot_length * f_max_idx], label='win')
# title(r'$\alpha$')
# legend()
#
# figure(3)
# plot(t_ref, E_ref)
# plot(t_ref, E_ref_win)
# plot(t_sam, E_sam)
# plot(t_sam, E_sam_win)
#
# figure(4)
# plot(f_ref[f_min_idx:plot_length * f_max_idx], E_ref_w[f_min_idx:plot_length * f_max_idx])
# plot(f_ref[f_min_idx:plot_length * f_max_idx], E_ref_w_win[f_min_idx:plot_length * f_max_idx])
# plot(f_sam[f_min_idx:plot_length * f_max_idx], E_sam_w[f_min_idx:plot_length * f_max_idx])
# plot(f_sam[f_min_idx:plot_length * f_max_idx], E_sam_w_win[f_min_idx:plot_length * f_max_idx])
# plot(f_ref[f_min_idx:plot_length * f_max_idx], abs(rfft(signal.unit_impulse(E_ref.size, 'mid'), n=nSamp_pow))[f_min_idx:plot_length * f_max_idx])

figure(5)
H_w = E_sam_w / E_ref_w

E_recon = 15 * irfft(H_w, n=t_ref.size)
max_idx = E_recon.size - argmax(abs(E_ref))
E_recon = concatenate((E_recon[max_idx:], E_recon[:max_idx]))
plot(t_ref, E_ref, lw=1)
E_delta = zeros(E_ref.size)
E_delta += signal.unit_impulse(E_ref.size, int(E_ref.size / 2))
E_delta += signal.unit_impulse(E_ref.size, int(1.05 * E_ref.size / 2.05)) / 2
plot(t_ref, E_delta, lw=1)
plot(t_ref, convolve(E_ref, E_delta, mode='same'), lw=1)


figure(6)
H_w = E_sam_w / E_ref_w
f_ref = f_ref[f_min_idx:plot_length * f_max_idx]

max_idx = E_recon.size - argmax(abs(E_ref))
E_recon = concatenate((E_recon[max_idx:], E_recon[:max_idx]))
plot(f_ref, abs(E_ref_w[f_min_idx:plot_length * f_max_idx]), lw=1)
E_delta = zeros(E_ref.size)
E_delta += signal.unit_impulse(E_ref.size, int(E_ref.size / 2))
E_delta += signal.unit_impulse(E_ref.size, int(1.05 * E_ref.size / 2.05)) / 2
plot(f_ref, abs(rfft(E_delta, n=nSamp_pow)[f_min_idx:plot_length * f_max_idx]), lw=1)
plot(f_ref, abs(rfft(convolve(E_ref, E_delta, mode='same'), n=nSamp_pow)[f_min_idx:plot_length * f_max_idx]), lw=1)

show()
