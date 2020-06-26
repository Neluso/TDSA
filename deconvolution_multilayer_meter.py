from TDSA import *
from TDS_constants import *
from DSP_functions import *
from pyswarm import pso
from read_data import read_1file
from tqdm import trange
from tkinter.ttk import *
from tkinter import *
from aux_functions import *
from numpy.fft import *
from scipy import signal
from matplotlib.pyplot import *


def objective_function(k, *args):
    k0, k1, k2, k3 = k
    E_val, E_meas, t_val = args
    idx_k1 = centre_loc(E_val) - int(round(k1 / mean(diff(t_val))))
    idx_k3 = centre_loc(E_val) - int(round(k3 / mean(diff(t_val))))
    return k0 * roll(E_val, - idx_k1) + k2 * roll(E_val, - idx_k3)


def min_function(k, *args):
    E_val, E_meas, t_val = args
    return sum((objective_function(k, *args) - E_meas)**2)


def constraints(k, *args):
    k0, k1, k2, k3 = k
    return [k0, abs(k0) - abs(k2), - k0 - k2 + 1, k3 - k1, abs(k2)]


# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_w_coat/ref metal wcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_w_coat/sam metal wcoat1_avg_f.txt')
t_ref, E_ref = read_1file(
    'data/muestras_airbus_boleto_176054_fecha_15_06_2018/cork_w_coat/ref metal cork wcoat_avg_f.txt')
t_sam, E_sam = read_1file(
    'data/muestras_airbus_boleto_176054_fecha_15_06_2018/cork_w_coat/sam cork wcoat3_avg_f.txt')

delta_t_ref = mean(diff(t_ref))
enlargement = (5 * 2**7 - 1) * E_ref.size
E_ref = zero_padding(E_ref, 0, enlargement)
t_ref = concatenate((t_ref, t_ref[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
E_sam = zero_padding(E_sam, 0, enlargement)
t_sam = concatenate((t_sam, t_sam[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
E_sam_max = amax(abs(E_sam))
E_sam_max_idx = where(abs(E_sam) == E_sam_max)[0][0]


f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
delta_f_ref = mean(diff(f_ref))



H_w = E_sam_w / E_ref_w
g_filt = gauss_low_filter(f_ref, f_ref[argmax(abs(E_sam_w))], 1.2)
wiener_filt = wiener_filter(E_ref_w, beta=fromDb(noise_floor(f_ref, E_ref_w, 20)))

H_w_filt = H_w * wiener_filt


irf_filt = irfft(H_w_filt, t_ref.size)
irf_filt_max = amax(irf_filt)
irf_filt_max_idx = where(irf_filt == irf_filt_max)[0][0]
t_irf = arange(irf_filt.size) / (irf_filt.size * delta_f_ref)

irf = roll(irf_filt / amax(abs(irf_filt)), E_sam_max_idx - irf_filt_max_idx)  # align deconvolution to t_ref
irf_dno = SWT_denoising(irf, 7, 0.1)  # apply de-noising

irf_peaks = signal.find_peaks(abs(irf), 0.125 * amax(abs(irf)))
print(irf_peaks)


figure(1)
plot(t_sam, E_sam / amax(E_sam), lw=1)
plot(t_ref, irf, lw=1)
plot(t_ref, irf_dno, lw=1)
for idx in irf_peaks[0]:  # mark peaks
    t_peak = t_ref[idx]
    E_peak = irf[idx]
    text(t_peak, E_peak, str(abs(round(E_peak, 2))))
    if E_peak >= 0:
        plot(t_peak, E_peak, 'r', marker=7)
    if E_peak < 0:
        plot(t_peak, E_peak, 'r', marker=6)
xlim([0, 50])

show()
