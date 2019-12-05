from read_data import read_data
from numpy import *
from DSP_functions import *
from TDS_constants import *
from aux_functions import *
from scipy.optimize import minimize
from scipy.interpolate import Akima1DInterpolator
from aux_functions import *
from DSP_functions import *
from pyswarm import pso


def objective_function(k0, k1, k2, k3, E_val_interpolated, t_val):
    idx_k1 = centre_loc(E_val_interpolated(t_val)) - int(round(k1 / mean(diff(t_val))))
    obj_fun = k0 * roll(E_val_interpolated(t_val), - idx_k1)
    idx_k3 = centre_loc(E_val_interpolated(t_val)) - int(round(k3 / mean(diff(t_val))))
    obj_fun += k2 * roll(E_val_interpolated(t_val), - idx_k3)
    return obj_fun


def min_function(params, *args):
    k0, k1, k2, k3 = params
    E_ref, E_meas, t_ref, t_meas = args
    return sum((objective_function(k0, k1, k2, k3, E_ref, t_ref) - E_meas(t_meas))**2)


def cons_0(x, *args):  # k0 > 0
    return x[0]


def cons_1(x, *args):  # k0 > k2 --- k0 - k2
    return abs(x[0]) - abs(x[2])


def cons_2(x, *args):  # k0 + k2 = 1 --- k0 + k2 - 1
    return x[0] + x[2] - 1


def cons_2_bis(x, *args):  # k0 + k2 = 1 --- k0 + k2 - 1 : bis constraint to achieve equality in pyswarm
    return - x[0] - x[2] + 1


def cons_3(x, *args):  # k3 > k1 --- k3 - k1 > 0
    return x[3] - x[1]


def constraints(x, *args):
    return cons_1(x, *args), cons_2_bis(x, *args), cons_3(x, *args), cons_2(x, *args),


fh = open('./data/Paintmeter/Trazas/1 capa de celofan en metal/ref metal.txt')
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
t_ref = array(t_ref)
E_ref = array(E_ref)

fh = open('./data/Paintmeter/Trazas/1 capa de celofan en metal/sam metal.txt')
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
t_sam = array(t_sam)
E_sam = array(E_sam)


nSamp = E_ref.size
nSamp_pow = nextpow2(nSamp)


E_ref_int = Akima1DInterpolator(t_ref, E_ref)
E_sam_int = Akima1DInterpolator(t_sam, E_sam)
enlargement = 1
t_ref = mean(diff(t_ref)) * arange(t_ref.size * enlargement) / enlargement
t_sam = mean(diff(t_sam)) * arange(t_sam.size * enlargement) / enlargement


k_min = [0, t_ref[0], -1, t_ref[0]]  # lower bounds fro k0, k1, k2, k3 respectively
k_max = [1, t_ref[-1], 1, t_ref[-1]]  # upper bounds fro k0, k1, k2, k3 respectively


xopt, fopt = pso(min_function, k_min, k_max,
                 args=(E_ref_int, E_sam_int, t_ref, t_sam),
                 swarmsize=200,
                 maxiter=200,
                 ieqcons=[constraints],
                 minstep=1e-10,
                 minfunc=1e-10,
                 debug=True)
k = xopt

figure(1)
title('Fit')
plot(t_sam, E_sam_int(t_sam), lw=1, label='sam')
plot(t_ref, objective_function(k[0], k[1], k[2], k[3], E_ref_int, t_ref), lw=1, label='fit')
legend()


print('k0 =', k[0])
print('k1 =', k[1])
print('k2 =', k[2])
print('k3 =', k[3])
delta_t = k[3] - k[1]
print('delta t =', delta_t, 'ps')
thickness = c_0 * delta_t * 1e-12 / (2 * 2.6)  # m
thickness *= 1e3  # mm
print('d =', thickness, 'mm')

show()
