from numpy import *
from scipy.optimize import minimize
from TDS_constants import *
from aux_functions import *
from DSP_functions import *
from matplotlib.pyplot import *
from matplotlib.axes import *
from mpl_toolkits.mplot3d import Axes3D
from tqdm import trange


def fp_m(f, m, n, k, L):
    n_cplx = n - 1j * k
    rho_term = (n_cplx - n_aire) / (n_cplx + n_aire)
    return (rho_term**2 * exp(-2j * n_cplx * f * 2 * pi * L / c_0))**m


def fp_full(f, n, k, L):
    return 1 / (1 - fp_m(f, 1, n, k, L))


def fp_first_der_m(f, m, n, k, ):
    n_cplx = n - 1j * k
    rho_prima_rho_term = (2 * n_aire) / (n_cplx ** 2 - n_aire ** 2)
    fp_val = 2 * m * (rho_prima_rho_term - 1j * 2 * pi * L / c_0)
    fp_val *= fp_m(f, m, n, k, L)
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
    fp_val = 0
    for i in range(fp_echo + 1):
        fp_val += fp_m(f, i, n, k, L)
    return n_quo * exp_term * fp_val  # n_quo * exp_term * fp_full(f, n, k, L)


def delta_min(n, k, L, f, H_w):  # f and H_w are sigle value. Frequency and measured measured H_w at f
    delta_min_val = ones([n.size, k.size])
    for n_idx in range(n.size):
        for k_idx in range(k.size):
            delta_min_val[n_idx][k_idx] = (log(abs(transfer_function(f, n[n_idx], k[k_idx], L, 1))) - log(abs(H_w)))**2
            delta_min_val[n_idx][k_idx] += (angle(transfer_function(f, n[n_idx], k[k_idx], L, 1)) - angle(H_w))**2
    return delta_min_val


# def grad_delta_min(n, k, L, f, H_w):
#     gradient()
#     return 0
#
#
# def hessian_delta_min(n, k, L, f, H_w):  # f and H_w are sigle value
#     hess = zeros([n.size, k.size])  # row, column
#     for n_idx in range(n.size):
#         for k_idx in range(k.size):
#             hess[n_idx][k_idx] = delta_min(n[n_idx], k[k_idx], L, f, H_w)
#     return 0


delta_x = 0.5
x_max = 10
n_x = delta_x * arange(delta_x, x_max / delta_x)
k_x = delta_x * arange(delta_x, x_max / delta_x)
L_x = delta_x * arange(delta_x, x_max / delta_x)
n_0 = 1
k_0 = delta_x
L = delta_x


# ----------------
# testing purposes
# ----------------

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


f_ref, E_ref_w = fourier_analysis(t_ref, E_ref, nSamp_pow)
f_sam, E_sam_w = fourier_analysis(t_sam, E_sam, nSamp_pow)


H_w = E_sam_w / E_ref_w
f_min_idx, f_max_idx = f_min_max_idx(f_ref)

f_ref = f_ref[f_min_idx:f_max_idx]
H_w = H_w[f_min_idx:f_max_idx]


delta_min_matrix = ones((f_ref.size, n_x.size, k_x.size))

for f_idx in trange(f_ref.size):
    delta_min_matrix[f_idx] = delta_min(n_x, k_x, 1, f_ref[f_idx], H_w[f_idx])


print(delta_min_matrix.shape)
