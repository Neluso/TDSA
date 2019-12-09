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


t_ref, E_ref = read_1file('./data/Paintmeter/Trazas/6 capas celofan en metal/ref metal celofan.txt')
t_sam, E_sam = read_1file('./data/Paintmeter/Trazas/6 capas celofan en metal/sam metal celofan.txt')

nSamp = E_ref.size
nSampPow = t_ref.size  # nextpow2(nSamp)
f_ref, E_ref_w = fourier_analysis(t_ref, E_ref, nSampPow)
f_sam, E_sam_w = fourier_analysis(t_sam, E_sam, nSampPow)
H_w = E_sam_w / E_ref_w
H_w_filt = lin_low_filter(H_w, 1/25, 1)
irf = ifft(H_w, t_ref.size)
irf_filt = ifft(H_w_filt, t_ref.size)


# figure(10)
# plot(f_ref, E_ref_w)
# plot(f_sam, E_sam_w)
# figure(20)
# plot(t_ref, irf)
# plot(t_ref, irf_filt)
# show()


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
repetitions = 100

for i in range(repetitions):
    
    print('Iteration', i + 1, 'of', repetitions)

    k, fopt = pso(min_function, k_min, k_max,
                  args=(E_ref, E_sam, t_ref),
                  swarmsize=1000,
                  maxiter=2000,
                  f_ieqcons=constraints,
                  phig=0.1,
                  phip=0.1,
                  minstep=1e-10,
                  minfunc=1e-10,
                  debug=False)


    # figure(1)
    # title('Fit')
    # plot(t_sam, E_sam, lw=1, label='sam')
    # plot(t_ref, objective_function(k, *(E_ref, E_sam, t_ref)), lw=1, label='fit')
    # legend()

    delta_t = k[3] - k[1]
    thickness = c_0 * delta_t * 1e-12 / (2 * 2.6)  # m
    thickness *= 1e3  # mm
    thick.append(thickness)
    error_func.append(fopt)
    k0.append(k[0])
    k1.append(k[1])
    k2.append(k[2])
    k3.append(k[3])
    # print('k0 =', k[0])
    # print('k1 =', k[1])
    # print('k2 =', k[2])
    # print('k3 =', k[3])
    # print('delta t =', delta_t, 'ps')
    # print('d =', thickness, 'mm')

figure(1)
plot(arange(repetitions), array(thick))
figure(2)
plot(arange(repetitions), array(error_func))
show()
