from TDS_constants import *
from DSP_functions import *
from pyswarm import pso
from read_data import read_1file
from tqdm import trange
from tkinter.ttk import *
from tkinter import *


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
    return [k0, abs(k0) - abs(k2), - k0 - k2 + 1, k3 - k1]



t_ref, E_ref = read_1file('./data/Paintmeter/Trazas/1 capa de celofan en vidrio/ref metal.txt')
t_sam, E_sam = read_1file('./data/Paintmeter/Trazas/1 capa de celofan en vidrio/ref vidrio.txt')


ref_pulse_idx = centre_loc(E_ref)
# lower bounds for k0, k1, k2, k3 respectively
k_min = [0, t_ref[0], -1, t_ref[ref_pulse_idx]]
# upper bounds for k0, k1, k2, k3 respectively
k_max = [1, t_ref[ref_pulse_idx], 0, t_ref[-1]]

k0 = list()
k1 = list()
k2 = list()
k3 = list()

thick = list()
repetitions = 100
popup = Toplevel()
popup.geometry('390x65')
popup.title('Working')
Label(popup).grid(row=0, column=0)
Label(popup).grid(row=1, column=0)
progress_bar = Progressbar(popup, orient="horizontal", length=380, mode="determinate", maximum=repetitions, value=0)
progress_bar.grid(row=1, column=1)
popup.pack_slaves()

for i in range(repetitions):
    popup.update()
    progress_bar['value'] += 1

    k, fopt = pso(min_function, k_min, k_max,
                  args=(E_ref, E_sam, t_ref),
                  swarmsize=5000,
                  maxiter=200,
                  f_ieqcons=constraints,
                  phig=2,
                  phip=2,
                  minstep=1e-15,
                  minfunc=1e-15,
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
plot(arange(repetitions), array(k0))
figure(3)
plot(arange(repetitions), array(k1))
figure(4)
plot(arange(repetitions), array(k2))
figure(5)
plot(arange(repetitions), array(k3))
show()