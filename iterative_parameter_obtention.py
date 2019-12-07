from numpy import *
from scipy.optimize import minimize
from TDS_constants import *
from aux_functions import *
from DSP_functions import *
from matplotlib.pyplot import *
from tqdm import trange
from pyswarm import pso


# function definitions

def fp_m(n, k, L, f, m):
    n_cplx = n - 1j * k
    rho_term = (n_cplx - n_aire) / (n_cplx + n_aire)
    return (rho_term**2 * exp(-2j * n_cplx * f * 2 * pi * L / c_0))**m


def fp_full(n, k, L, f):
    return 1 / (1 - fp_m(n, k, L, f, 1))


def fp_first_der_m(n, k, L, f, m):
    n_cplx = n - 1j * k
    rho_prima_rho_term = (2 * n_aire) / (n_cplx ** 2 - n_aire ** 2)
    fp_val = 2 * m * (rho_prima_rho_term - 1j * 2 * pi * L / c_0)
    fp_val *= fp_m(n, k, L, f, m)
    # fp_val *= fp_full(n, k, L, f)
    return fp_val


def fp_second_der_m(f, m, n, k, L):
    n_cplx = n - 1j * k
    rho_prima_rho_term = (2 * n_aire) / (n_cplx**2 - n_aire**2)
    fp_val = (2 * m * (rho_prima_rho_term - 1j * 2 * pi * f * L / c_0))**2
    fp_val -= 8 * m * (4 * n_cplx * n_aire / (n_cplx**2 - n_aire**2))
    fp_val *= fp_m(f, m, n, k, L)
    return fp_val


def transfer_function(f, n, k, L, fp_echo):
    n_cplx = n - 1j * k
    n_quo = 4 * n_cplx * n_aire / (n_cplx + n_aire)**2
    exp_term = exp(- 1j * (n_cplx - n_aire) * (2 * pi * L) / c_0)
    full_echos = False
    if full_echos:
        fp_val = fp_full(n, k, L, f)
    else:
        fp_val = 0
        for i in range(fp_echo + 1):
            fp_val += fp_m(f, i, n, k, L)
    return n_quo * exp_term * fp_val  # n_quo * exp_term * fp_full(f, n, k, L)


def delta_min(params, *data):  # f and H_w are single value. Frequency and measured measured H_w at f
    n, k, L = params
    f, H_w = data
    delta_min_val = (log(abs(transfer_function(f, n, k, L, 1))) - log(abs(H_w)))**2
    delta_min_val += (angle(transfer_function(f, n, k, L, 1)) - angle(H_w))**2
    return delta_min_val


# Constraints definition

def n_cons(x):  # to be used as ineq in the minimize routine
    n = x[0]
    return n - 1


def k_cons(x):  # to be used as ineq in the minimize routine
    k = x[1]
    return k


def L_cons(x):  # to be used as ineq in the minimize routine
    L = x[2]
    return L


# main script

fh = open('./data/demo_data/test_ref.txt')
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

fh = open('./data/demo_data/test_sam.txt')
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


f_ref, E_ref_w = fourier_analysis(t_ref, E_ref, nSamp_pow)
f_sam, E_sam_w = fourier_analysis(t_sam, E_sam, nSamp_pow)


H_w = E_sam_w / E_ref_w
f_min_idx, f_max_idx = f_min_max_idx(f_ref)

f_ref = f_ref[f_min_idx:f_max_idx]
H_w = H_w[f_min_idx:f_max_idx]


n_0 = 1.1
k_0 = 0.1
L_0 = 1 * 1e-3


n_opt = zeros(f_ref.size)
k_opt = zeros(f_ref.size)
L_opt = zeros(f_ref.size)

tl = 0.0000001  # tolerance
bnds = ((1, None), (tl, None))  #, (tl, None))
cons = ({'type': 'ineq', 'fun': n_cons}, {'type': 'ineq', 'fun': k_cons})  # , {'type': 'ineq', 'fun': L_cons})

for f_idx in trange(f_ref.size):
    # res = minimize(delta_min, array((n_0, k_0)), args=(L_0, f_ref[f_idx], H_w[f_idx]), tol=tl, bounds=bnds)
    xopt, fopt = pso(delta_min, [1, tl, tl], [100, 100, 1000], args=(f_ref[f_idx], H_w[f_idx]),
                     swarmsize=250,
                     maxiter=500,
                     minstep=1e-12,
                     minfunc=1e-12,
                     debug=False
                     )
    n_opt[f_idx] = xopt[0]  # res.x[0]
    k_opt[f_idx] = xopt[1]  # res.x[1]
    # L_opt[f_idx] = res.x[2]

figure(1)
plot(f_ref, n_opt)
xlabel(r'$f\ (Hz)$')
ylabel(r'$n$')
figure(2)
plot(f_ref, k_opt)
xlabel(r'$f\ (Hz)$')
ylabel(r'$k$')
# figure(3)
# plot(f_ref, L_opt)
# xlabel(r'$f\ (Hz)$')
# ylabel(r'$L$')
print(Lxopt[2])
show()