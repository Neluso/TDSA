from DSP_functions import *
from TDS_constants import *
from numpy import *
from aux_functions import *
from read_data import *
from matplotlib.pyplot import *
from scipy.stats import cauchy
from pyswarm import pso
from tqdm import *
from scipy.optimize import minimize


def thickness_val(n, order_m, f_ref, f_sam):
    return c_0 * (order_m + 0.5) / (n - n_aire) * (1 / f_ref - 1 / f_sam)


def thickness_sigma2(n, sig_n, order_m, f_ref, sig_f_ref, f_sam, sig_f_sam):  # sigs in relative error (0 to 1)
    aux1 = ((1 / f_ref - 1 / f_sam) * c_0 * (order_m + 0.5) / (n - n_aire)**2)**2 * (sig_n * n)**2
    aux2 = (c_0 * (order_m + 0.5) / (n - n_aire))**2 * ((sig_f_ref / f_ref)**2 + (sig_f_sam / f_sam)**2)
    return aux1 + aux2


def objetive_func(x, loc, scale, a, m):
    return cauchy.pdf(x, loc, scale) + a * x + m


def cauchy_fit(params, *args):
    loc, scale, a, m = params
    x, y = args
    return sum((cauchy.pdf(x, loc, scale) + a * x + m - y)**2)


t_ref, E_ref = read_slow_data('./data/notch/ref_notch.txt')
t_sam, E_sam = read_slow_data('./data/notch/sam_notch.txt')

E_ref = - flip(E_ref)
E_sam = - flip(E_sam)


delta_t_ref = mean(diff(t_ref))
enlargement = 0 * E_ref.size
E_ref = zero_padding(E_ref, 0, enlargement)
t_ref = concatenate((t_ref, t_ref[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
E_sam = zero_padding(E_sam, 0, enlargement)
t_sam = concatenate((t_sam, t_sam[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))

f_ref, E_ref_w = fourier_analysis(t_ref, E_ref, E_ref.size)
f_sam, E_sam_w = fourier_analysis(t_sam, E_sam, E_sam.size)


E_max = amax(abs(E_ref_w))
E_ref_w, E_sam_w = prettyfy(E_ref_w, E_max), prettyfy(E_sam_w, E_max)
f_min_idx, f_max_idx = f_min_max_idx(f_ref, 0.225e-12, 0.45e-12)
x_v = f_ref[f_min_idx:f_max_idx]
y_v = - E_ref_w[f_min_idx:f_max_idx]
x_p = f_sam[f_min_idx:f_max_idx]
y_p = - E_sam_w[f_min_idx:f_max_idx]


param_min = (0.225, 0, 0, -10)
param_max = (0.4, 10, 50, 10)
debug_pso = False
f_va_fits = list()
f_vp_fits = list()

# res_v, fopt = pso(cauchy_fit, param_min, param_max,
#                   swarmsize=500,
#                   maxiter=200,
#                   args=(x_v, y_v),
#                   debug=debug_pso)
#
# res_p, fopt = pso(cauchy_fit, param_min, param_max,
#                   swarmsize=500,
#                   maxiter=200,
#                   args=(x_p, y_p),
#                   debug=debug_pso)

res_v = minimize(cauchy_fit, array((0.3, 0.01, 24, -5)), args=(x_v, y_v), method='Nelder-Mead', options={'maxiter':2000})
res_p = minimize(cauchy_fit, array((0.3, 0.03, 24, 0)), args=(x_p, y_p), method='Nelder-Mead')

print(res_v)
print(res_p)


# sigma_va = res_v2.jac
# sigma_vp = res_p2.jac
# sigma_va = dot(sigma_va.T, sigma_va) * sqrt(res_v2.fun / (x_v.size - 4))
# sigma_vp = dot(sigma_vp.T, sigma_vp) * sqrt(res_p2.fun / (x_p.size - 4))

f_va = res_v.x[0]
f_vp = res_p.x[0]
sigma_va = res_v.x[1]
sigma_vp = res_p.x[1]


f_va *= 1e12
sigma_va *= 1e12
f_vp *= 1e12
sigma_vp *= 1e12


n_p = 1.45
d_p = thickness_val(n_p, 1, f_va, f_vp)
sigma_d_p = thickness_sigma2(n_p, 0.1, 1, f_va, sigma_va / f_va, f_vp, sigma_vp / f_vp)
print('va = ', f_va, sigma_va, sigma_va / f_va)
print('vp = ', f_vp, sigma_vp, sigma_vp / f_vp)
print('Paper thickness =', d_p * 1e3, 'mm')
print('Thickness error =', 1e3 * sqrt(sigma_d_p), 'mm')

# figure(1)
# plot(t_ref, E_ref, lw=1)
# plot(t_sam, E_sam, lw=1)
# figure(2)
# plot(f_ref, E_ref_w, lw=1)
# plot(f_sam, E_sam_w, lw=1)
# figure(3)
# plot(x_v, y_v)
# plot(x_p, y_p)
# plot(x_v, cauchy.pdf(x_v, res_v[0], res_v[1]) + res_v[2] * x_v + res_v[3])
# plot(x_p, cauchy.pdf(x_p, res_p[0], res_p[1]) + res_p[2] * x_p + res_p[3])
# show()
