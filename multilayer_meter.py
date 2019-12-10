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


t_ref, E_ref = read_1file('./data/Paintmeter/Trazas/1 capa de celofan en metal/ref metal.txt')
t_sam, E_sam = read_1file('./data/Paintmeter/Trazas/1 capa de celofan en metal/sam metal.txt')

delta_t_ref = mean(diff(t_ref))
enlargement = 20 * E_ref.size
E_ref = zero_padding(E_ref, 0, enlargement)
t_ref = concatenate((t_ref, t_ref[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
E_sam = zero_padding(E_sam, 0, enlargement)
t_sam = concatenate((t_sam, t_sam[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
E_sam_max = amax(abs(E_sam))
E_sam_max_idx = where(abs(E_sam) == E_sam_max)[0][0]

nSamp = E_ref.size
nSampPow = nSamp  # nextpow2(nSamp)
f_ref, E_ref_w = fourier_analysis(t_ref, E_ref, nSampPow)
f_sam, E_sam_w = fourier_analysis(t_sam, E_sam, nSampPow)
delta_f_ref = mean(diff(f_ref))


H_w = E_sam_w / E_ref_w
g_filt = gauss_low_filter(f_ref, f_ref[argmax(abs(E_sam_w))], 1.2)
H_w_filt = H_w * g_filt

irf_filt = irfft(H_w_filt, t_ref.size)
irf_filt_max = amax(irf_filt)
irf_filt_max_idx = where(irf_filt == irf_filt_max)[0][0]
t_irf = arange(irf_filt.size) / (irf_filt.size * delta_f_ref)

irf = roll(irf_filt / amax(abs(irf_filt)), E_sam_max_idx - irf_filt_max_idx)
deconv_threshold = 0.25
irf_maxs = where(abs(irf) >= deconv_threshold)[0]
print(irf_maxs)


figure(10)
plot(t_sam, E_sam / amax(E_sam), lw=1)
plot(t_ref, irf, lw=1)
plot(t_ref, 0.25 * t_ref / t_ref, lw=1)
xlim([0, 70])
figure(20)
# plot(f_ref, H_w, 'b-', lw=1)
plot(f_ref, toDb(H_w), 'b-', lw=1)
# plot(f_ref, H_w_filt, 'g-', lw=1)
plot(f_ref, toDb(H_w_filt), 'g-', lw=1)
plot(f_ref, toDb(g_filt), 'r-', lw=1)
plot(f_ref, zeros(f_ref.size), 'y-', lw=0.3)
xlim([0, 2])
ylim([-1, 3])
show()

ref_pulse_idx = centre_loc(E_ref)
# lower bounds for k0, k1, k2, k3 respectively
k_min = [0, t_ref[0], -1, t_ref[ref_pulse_idx]]
# upper bounds for k0, k1, k2, k3 respectively
k_max = [1, t_ref[ref_pulse_idx], 1, t_ref[-1]]

k0 = list()
k1 = list()
k2 = list()
k3 = list()

thick = list()
error_func = list()
repetitions = 10
#
# for i in range(repetitions):
#
#     print('Iteration', i + 1, 'of', repetitions)
#
#     k, fopt = pso(min_function, k_min, k_max,
#                   args=(E_ref, E_sam, t_ref),
#                   swarmsize=1000,
#                   maxiter=2000,
#                   f_ieqcons=constraints,
#                   phig=0.1,
#                   phip=0.1,
#                   minstep=1e-10,
#                   minfunc=1e-10,
#                   debug=False)
#
#
#     # figure(1)
#     # title('Fit')
#     # plot(t_sam, E_sam, lw=1, label='sam')
#     # plot(t_ref, objective_function(k, *(E_ref, E_sam, t_ref)), lw=1, label='fit')
#     # legend()
#
#     delta_t = k[3] - k[1]
#     thickness = c_0 * delta_t * 1e-12 / (2 * 2.6)  # m
#     thickness *= 1e3  # mm
#     thick.append(thickness)
#     error_func.append(fopt)
#     k0.append(k[0])
#     k1.append(k[1])
#     k2.append(k[2])
#     k3.append(k[3])
#     # print('k0 =', k[0])
#     # print('k1 =', k[1])
#     # print('k2 =', k[2])
#     # print('k3 =', k[3])
#     # print('delta t =', delta_t, 'ps')
#     # print('d =', thickness, 'mm')
#
# figure(1)
# plot(arange(repetitions), array(thick))
# figure(2)
# plot(arange(repetitions), array(error_func))
# show()
