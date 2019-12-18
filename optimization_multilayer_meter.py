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
from scipy.optimize import differential_evolution


def objective_function(k, *args):
    k0, k1, k2, k3, k4, k5, k6, k7 = k
    E_val, E_meas, t_val = args
    idx_k1 = centre_loc(E_val) - int(round(k1 / mean(diff(t_val))))
    idx_k3 = centre_loc(E_val) - int(round(k3 / mean(diff(t_val))))
    idx_k5 = centre_loc(E_val) - int(round(k5 / mean(diff(t_val))))
    idx_k7 = centre_loc(E_val) - int(round(k7 / mean(diff(t_val))))
    return k0 * roll(E_val, - idx_k1) + k2 * roll(E_val, - idx_k3) + k4 * roll(E_val, - idx_k5) + k6 * roll(E_val, - idx_k7)


def min_function(k, *args):
    E_val, E_meas, t_val = args
    return sum((objective_function(k, *args) - E_meas)**2)


def constraints(k, *args):
    k0, k1, k2, k3 = k
    return [k0, abs(k0) - abs(k2), - k0 - k2 + 1, k3 - k1, abs(k2)]


t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_w_coat/ref metal wcoat_avg_f.txt')
t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_w_coat/sam metal wcoat1_avg_f.txt')
# t_ref, E_ref = read_1file('./data/old/Paintmeter/Trazas/6 capas celofan en metal/ref metal celofan.txt')
# t_sam, E_sam = read_1file('./data/old/Paintmeter/Trazas/6 capas celofan en metal/sam metal celofan.txt')

delta_t_ref = mean(diff(t_ref))
enlargement = 0 * E_ref.size
E_ref = zero_padding(E_ref, 0, enlargement)
t_ref = concatenate((t_ref, t_ref[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
E_sam = zero_padding(E_sam, 0, enlargement)
t_sam = concatenate((t_sam, t_sam[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
E_sam_max = amax(abs(E_sam))
E_sam_max_idx = where(abs(E_sam) == E_sam_max)[0][0]


f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
delta_f_ref = mean(diff(f_ref))


ref_pulse_idx = centre_loc(E_ref)
# lower bounds for k0, k1, k2, k3 respectively
k_min = [0, t_ref[0], -1, t_ref[ref_pulse_idx]]
# upper bounds for k0, k1, k2, k3 respectively
k_max = [1, t_ref[ref_pulse_idx], 1, t_ref[-1]]

# k_bounds = [(-1, 1), (t_ref[0], t_ref[ref_pulse_idx]), (-1, 1), (t_ref[ref_pulse_idx], t_ref[-1])]
k_bounds = [(-1, 1), (t_ref[0], t_ref[-1]),  # k0, k1
            (-1, 1), (t_ref[0], t_ref[-1]),  # k2, k3
            (-1, 1), (t_ref[0], t_ref[-1]),  # k4, k5
            (-1, 1), (t_ref[0], t_ref[-1])   # k6, k7
            ]


res = differential_evolution(min_function,
                             k_bounds,
                             args=(E_ref, E_sam, t_ref)
                             )

print(res)
k = res.x
pulses = list()
for i in range(int(round(k.size / 2))):
    pulses.append((k[2 * i + 1], k[2 * i]))
pulses = sorted(pulses)
delta_t = abs(pulses[0][0] - pulses[1][0])
thickness = c_0 * delta_t * 1e-12 / (2 * 1.51 * cos(pi / 6))  # m
thickness *= 1e3  # mm
print(thickness, 'mm')

plot(t_sam, E_sam, lw=1)
E_sim = objective_function(k, E_ref, E_sam, t_ref)
plot(t_ref, E_sim, lw=1)
show()
