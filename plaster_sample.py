from TDSA import *
from scipy.optimize import curve_fit, differential_evolution
from scipy.signal.windows import tukey
from time import time_ns
from genetic_denoising import genetic_deno


# function definitions
def ct(n_1, n_2):
    return 2 * n_1 / (n_1 + n_2)


def cr(n_1, n_2):
    return (n_2 - n_1) / (n_1 + n_2)


def phase_factor(n, thick, freq):  # theta in radians
    omg = 2 * pi * freq
    phi = omg * thick / c_0
    return np.exp(- 1j * n * phi)


def fabry_perot(n_in, n_out, n_l, d_l, freq):
    return 1 / (1 - cr(n_in, n_l) * cr(n_l, n_out) * phase_factor(n_l, 2*d_l, freq))


def epsilon(e_s, e_inf, tau, freq):  # Debye model
    omg = 2 * pi * freq
    e_w = e_inf + (e_s - e_inf) / (1 + 1j * omg * tau)
    return e_w


def nk_from_eps(e_s, e_inf, tau, freq):
    e_w = epsilon(e_s, e_inf, tau, freq)
    n = sqrt((abs(e_w) + real(e_w)) / 2)
    k = sqrt((abs(e_w) - real(e_w)) / 2)
    return n, k


def H_sim(freq, params):
    d_air, e_s_1, e_inf_1, tau_1, d_1, e_s_2, e_inf_2, tau_2, d_2 = params
    n_1, k_1 = nk_from_eps(e_s_1, e_inf_1, tau_1, freq)
    n_2, k_2 = nk_from_eps(e_s_2, e_inf_2, tau_2, freq)

    H_0 = phase_factor(n_air_cplx, d_air, freq)

    # Layer 1
    H_1 = ct(n_air_cplx, n_1 - 1j * k_1) * ct(n_1 - 1j * k_1, n_2 - 1j * k_2) * phase_factor(n_1 - 1j * k_1, d_1, freq)
    H_1 *= fabry_perot(n_air_cplx, n_2 - 1j * k_2, n_1 - 1j * k_1, d_1, freq)

    #Layer 2
    # H_2 = ct(n_1 - 1j * k_1, n_2 - 1j * k_2) * ct(n_2 - 1j * k_2, n_3 - 1j * k_3) * phase_factor(n_2 - 1j * k_2, d_2, freq)
    H_2 = H_1 * ct(n_2 - 1j * k_2, n_air_cplx) * phase_factor(n_2 - 1j * k_2, d_2, freq)
    H_2 *= fabry_perot(n_1 - 1j * k_1, n_air_cplx, n_2 - 1j * k_2, d_2, freq)
    
    return H_0 * H_1 * H_2


def cost_function(params, *args):
    E_sam, E_ref_w, freqs = args
    H_teo = H_sim(freqs, params)
    E_sam_teo_w = E_ref_w * H_teo
    E_sam_teo = irfft(E_sam_teo_w, n=E_sam.size)
    return sum((E_sam_teo - E_sam) ** 2)


# Plaster
t_ref, E_ref = read_1file('./data/2_layer/ref_plaster.txt')
t_sam, E_sam = read_1file('./data/2_layer/sam_plaster.txt')
t_vid, E_vid = read_1file('./data/2_layer/vid_plaster.txt')

n, alpha_f, n_avg = jepsen_index(t_ref, E_ref, t_vid, E_vid, 2e-3)

t_ref *= 1e-12
t_sam *= 1e-12
t_vid *= 1e-12

f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
f_vid, E_vid_w = fourier_analysis(t_vid, E_vid)

H_w = E_sam_w / E_ref_w


alpha = 1
beta = 10
k_bounds = [  # calibration
    (-5e-3, 5e-3),  # air thickness
    (alpha, beta), (alpha, beta), (0, 1e-13), (1e-6, 3000e-6),
    (alpha, beta), (alpha, beta), (0, 1e-13), (1e-6, 3000e-6)
]
t1 = time_ns()
res = differential_evolution(cost_function,
                             k_bounds,
                             args=(E_sam, E_ref_w, f_ref),
                             popsize=45,
                             maxiter=2000,
                             disp=True,  # step cost_function value
                             polish=True
                             )
t2 = time_ns()
# print()
# d_air = res.x[0]
# print(res)
n1, k1 = nk_from_eps(res.x[1], res.x[2], res.x[3], f_ref)
# print('White coat -', 'n:', round(mean(n1), 2), 'k:', round(mean(k1), 2), 'd:', round(res.x[4] * 1e6, 0), 'um')
# print('\t\t e_s:', round(res.x[2], 2), 'e_inf', round(res.x[1], 2), 'tau', res.x[3])
n2, k2 = nk_from_eps(res.x[5], res.x[6], res.x[7], f_ref)
# print('Green coat -', 'n:', round(mean(n2), 2), 'k:', round(mean(k2), 2), 'd:', round(res.x[8] * 1e6, 0), 'um')
# print('\t\t e_s:', res.x[6], 'e_inf:', res.x[5], 'tau:', res.x[7])
# print('Total:', round((res.x[4] + res.x[8]) * 1e6, 0), 'um')
# print()
# print('Air - d: ', round(d_air * 1e6, 2), 'um')
secs = (t2-t1)*1e-9
if secs < 3600:
    print('Processing time (mm:ss):', strftime('%M:%S', gmtime(secs)))
else:
    print('Processing time (hh:mm:ss):', strftime('%H:%M:%S', gmtime(secs)))
print()


H_fit = H_sim(f_ref, res.x)
E_fit = irfft(H_fit * E_ref_w)

f_min_idx, f_max_idx = f_min_max_idx(f_ref, 0.05, 1.5)


t_sam *= 1e12
t_ref *= 1e12
f_ref *= 1e-12
f_sam *= 1e-12

f_ref = f_ref[f_min_idx:f_max_idx]
H_w = H_w[f_min_idx:f_max_idx]
H_fit = H_fit[f_min_idx:f_max_idx]


figure(1)
plot(t_sam, E_sam, label='sam')
plot(t_sam, E_fit, label='fit')
xlim([t_sam[0], t_sam[-1]])
xlabel('t (ps)')
legend()

n1 = n1[f_min_idx:f_max_idx]
k1 = k1[f_min_idx:f_max_idx]
n2 = n2[f_min_idx:f_max_idx]
k2 = k2[f_min_idx:f_max_idx]
# figure(2)
# plot(f_ref, abs(H_w), label='sam')
# plot(f_ref, abs(H_fit), label='fit')
# xlim([f_ref[0], f_ref[-1]])
# ylabel('Abs')
# xlabel('f (THz)')
# legend()
#
# figure(3)
# plot(f_ref[f_min_idx:f_max_idx], unwrap(angle(H_w))[f_min_idx:f_max_idx], label='sam')
# plot(f_ref[f_min_idx:f_max_idx], unwrap(angle(H_fit))[f_min_idx:f_max_idx], label='fit')
# xlim([f_ref[0], f_ref[-1]])
# ylabel('Arg')
# xlabel('f (THz)')
# legend()
#
figure(4)
plot(f_ref, n1, label='pls')
plot(f_ref, n2, label='vid')
title('n')
xlim([f_ref[0], f_ref[-1]])
ylabel('n')
xlabel('f (THz)')
legend()
#
# figure(5)
# plot(f_ref, k1, label='pls')
# plot(f_ref, k2, label='vid')
# title('k')
# xlim([f_ref[0], f_ref[-1]])
# ylabel('k')
# xlabel('f (THz)')
# legend()

show()
