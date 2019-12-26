from scipy.optimize import differential_evolution, minimize
from numpy.random import rand
from TDSA import *


# function definitions

def fp_m(n, k, L, f, m):
    n_cplx = n + 1j * k
    rho_term = (n - n_air) / (n + n_air)
    return ((rho_term**2) * exp(2j * n_cplx * f * 2 * pi * L / c_0))**m


def fp_full(n, k, L, f):
    return 1 / (1 - fp_m(n, k, L, f, 1))


def transfer_function(f, n, k, L, fp_echo):
    n_cplx = n + 1j * k
    n_quo = 4 * n * n_air / (n + n_air) ** 2
    exp_term = exp(1j * (n_cplx - n_air) * (2 * pi * L) / c_0)
    full_echos = False
    if full_echos:
        fp_val = fp_full(n, k, L, f)
    else:
        fp_val = 0
        for i in range(fp_echo + 1):
            fp_val += fp_m(n, k, L, f, i)
    return n_quo * exp_term * fp_val


def delta_min(params, *args):  # f and H_w are single value. Frequency and measured measured H_w at f
    n, k, L = params
    f, H_w = args
    delta_min_val = (log(abs(transfer_function(f, n, k, L, 2))) - log(abs(H_w)))**2
    delta_min_val += (unwrap(angle(transfer_function(f, n, k, L, 1))) - unwrap(angle(H_w)))**2
    return sum(delta_min_val)


# main script
t_ref, E_ref = read_1file('./data/demo_data/test_ref.txt')
t_sam, E_sam = read_1file('./data/demo_data/test_sam.txt')
t_ref *= 1e-12
t_sam *= 1e-12

f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
f_min_idx, f_max_idx = f_min_max_idx(f_ref)

f_ref = f_ref[f_min_idx:f_max_idx]
E_ref_w = E_ref_w[f_min_idx:f_max_idx]
E_sam_w = E_sam_w[f_min_idx:f_max_idx]

H_w = E_sam_w / E_ref_w

tl = 1e-8  # tolerance
L_0 = 1.95*1e-3
n_0 = 1 + c_0 * (t_sam[centre_loc(E_sam)] - t_ref[centre_loc(E_ref)]) / L_0
k_0 = - c_0 * log(abs(amax(E_sam_w)) / abs(amax(E_ref_w))) / (L_0 * 2 * pi * f_ref[argmax(E_ref_w)])
params_0 = array([n_0, k_0, L_0])
bnds = [(0.99*n_0, 1.01*n_0), (0.99*k_0, 1.01*k_0), (0.9*L_0, 1.1*L_0)]
# initial_values = list()
# particle_number = 15*3
# for i in range(particle_number):
#     initial_values.append(array([n_0 + -0.5 * (-1 + 2 * rand()),
#                                  k_0 + 0.01 * (-1 + 2 * rand()),
#                                  L_0 + 0.5*1e-5 * (-1 + 2 * rand())
#                                  ]))
# initial_values = array(initial_values)
n_opt = zeros(f_ref.size)

k_opt = zeros(f_ref.size)
L_opt = zeros(f_ref.size)
for i in trange(f_ref.size):
    res = differential_evolution(delta_min,
                             bnds,
                             args=(f_ref, H_w),
                             strategy='best1bin',
                             popsize=150,
                             # init=initial_values,
                             maxiter=2000,
                             polish=True
                             )
    # res = minimize(delta_min, params_0, args=(f_ref, H_w), method='Nelder-Mead')
    n_opt[i] = res.x[0]
    k_opt[i] = res.x[1]
    L_opt[i] = res.x[2]


# Plots
# # n, alpha, n_avg = jepsen_index(t_ref, E_ref, t_sam, E_sam, 1.95*1e-3)  # 2*res.x[2])
# # n = n[f_min_idx:f_max_idx]
# # alpha = 0.01 * alpha[f_min_idx:f_max_idx]
f_ref *= 1e-12


fig, axs = subplots(2)
fig.suptitle('Abs/Phase')
axs[0].plot(f_ref, abs(H_w), lw=1)
axs[0].plot(f_ref, abs(transfer_function(f_ref, n_opt, k_opt, L_opt, 1)), lw=1)
axs[0].plot(f_ref, abs(transfer_function(f_ref, n_0, k_0, L_0, 1)), lw=1)
axs[0].set_ylabel(r'$\rho$')
axs[0].xaxis.set_visible(False)
axs[0].set_xlim([f_ref[0], f_ref[-1]])
axs[1].plot(f_ref, unwrap(angle(H_w)), lw=1)
axs[1].plot(f_ref, unwrap(angle(transfer_function(f_ref, n_opt, k_opt, L_opt, 1))), lw=1)
axs[1].plot(f_ref, unwrap(angle(transfer_function(f_ref, n_0, k_0, L_0, 1))), lw=1)
axs[1].set_ylabel(r'$\phi$')
axs[1].set_xlim([f_ref[0], f_ref[-1]])
xlabel(r'$f\ (THz)$')


fig, axs = subplots(3)
fig.suptitle('n_opt, k_opt, L_opt')
axs[0].plot(f_ref, n_opt, lw=1)
axs[0].set_ylabel('n')
axs[0].xaxis.set_visible(False)
axs[0].set_xlim([f_ref[0], f_ref[-1]])
axs[1].plot(f_ref, k_opt, lw=1)
axs[1].set_ylabel('k')
axs[1].xaxis.set_visible(False)
axs[1].set_xlim([f_ref[0], f_ref[-1]])
axs[2].plot(f_ref, L_opt, lw=1)
axs[2].set_ylabel('L')
axs[2].set_xlim([f_ref[0], f_ref[-1]])
xlabel(r'$f\ (THz)$')


fig, axs = subplots(2)
fig.suptitle('Re/Im')
axs[0].plot(f_ref, real(H_w), lw=1)
axs[0].plot(f_ref, real(transfer_function(f_ref, res.x[0], res.x[1], res.x[2], 1)), lw=1)
axs[1].plot(f_ref, imag(H_w), lw=1)
axs[1].plot(f_ref, imag(transfer_function(f_ref, res.x[0], res.x[1], res.x[2], 1)), lw=1)
xlabel(r'$f\ (THz)$')

# fig, axs = subplots(2)
# fig.suptitle('Optical parameters')
# axs[0].plot(f_ref, n, lw=1)
# axs[0].set_ylabel(r'$n$')
# axs[1].plot(f_ref, alpha, lw=1)
# axs[1].set_ylabel(r'$\alpha \ (cm^{-1})$')
# xlabel(r'$f\ (THz)$')

show()
