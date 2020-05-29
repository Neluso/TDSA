from TDSA import *
from numpy.fft import *
from numpy.fft import fft as fft_func
from time import time_ns
import lmfit


# constant definitions
deg_in = 30  # incidence angle in degrees
snell_sin = n_air * sin(deg_in * pi / 180)
layers = 3
# epsilon_model = 'debye'
# epsilon_model = 'cole'
# n_subs = 1.17 - 0.0 * 1j  # substrate refractive index -- cork
# n_subs = 1.17e20 - 0.0 * 1j  # substrate refractive index -- metal
n_subs = 1.25 - 0.000 * 1j  # substrate refractive index -- cork 2.0


def theta(n):
    return arcsin(snell_sin / real(n))


def ct2(n_l, n_l_1):
    n_l *= cos(theta(n_l))
    n_l_1 *= cos(theta(n_l_1))
    return 4 * n_l * n_l_1 / (n_l + n_l_1)**2


def cr_l_1_l(n_l, n_l_1):  # from n_l-1 to n_l
    n_l_1 *= cos(theta(n_l_1))
    n_l *= cos(theta(n_l))
    return - (n_l_1 - n_l) / (n_l_1 + n_l)
    # return (n_l_1 - n_l) / (n_l_1 + n_l)


def phase_factor(n, k, thick, freq):  # theta in radians
    omg = 2 * pi * freq
    thick *= cos(theta(n))
    phi = 2 * omg * thick / c_0
    return exp(- 1j * n * phi) * exp(- k * phi)


def epsilon(e_s, e_inf, tau, freq):  # Debye model
    omg = 2 * pi * freq
    e_w = e_inf + (e_s - e_inf) / (1 + 1j * omg * tau)
    return e_w


def epsilon_cc(e_s, e_inf, tau, alpha, freq):  # Cole-Cole mode
    omg = 2 * pi * freq
    e_w = e_inf + (e_s - e_inf) / (1 + (1j * omg * tau)**(alpha))
    return e_w


def nk_from_eps(e_inf, e_s, tau, freq):
    e_w = epsilon(e_s, e_inf, tau, freq)
    n = sqrt((abs(e_w) + real(e_w)) / 2)
    k = sqrt((abs(e_w) - real(e_w)) / 2)
    return n, k


def H_sim(freq, d_air, n_1, k_1, d_1, n_2, k_2, d_2, n_3, k_3, d_3):

    H_i = cr_l_1_l(n_subs, n_3 - 1j * k_3) * ones(freq.size)
    rlm1l = cr_l_1_l(n_3 - 1j * k_3, n_2 - 1j * k_2)
    tt = ct2(n_3 - 1j * k_3, n_3 - 1j * k_2)
    exp_phi = phase_factor(n_3, k_3, d_3, freq)
    H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)

    rlm1l = cr_l_1_l(n_2 - 1j * k_2, n_1 - 1j * k_1)
    tt = ct2(n_2 - 1j * k_2, n_1 - 1j * k_1)
    exp_phi = phase_factor(n_2, k_2, d_2, freq)
    H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)

    rlm1l = cr_l_1_l(n_1 - 1j * k_1, n_air_cplx)
    tt = ct2(n_1 - 1j * k_1, n_air_cplx)
    exp_phi = phase_factor(n_1, k_1, d_1, freq)
    H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)

    return phase_factor(n_air, 0, d_air, freq) * H_i


def resid(params, *args):
    d_air = params['d_air'].value
    # e_inf_1 = params['e_inf_1'].value
    # e_s_1 = params['e_s_1'].value
    # tau_1 = params['tau_1'].value
    d_1 = params['d_1'].value
    # e_inf_2 = params['e_inf_2'].value
    # e_s_2 = params['e_s_2'].value
    # tau_2 = params['tau_2'].value
    d_2 = params['d_2'].value
    # e_inf_3 = params['e_inf_3'].value
    # e_s_3 = params['e_s_3'].value
    # tau_3 = params['tau_3'].value
    d_3 = params['d_3'].value
    freq, E_r_w, E_meas, n1, k1, n2, k2, n3, k3 = args
    # n_1, k_1 = nk_from_eps(e_inf_1, e_s_1, tau_1, freq)
    # n_2, k_2 = nk_from_eps(e_inf_2, e_s_2, tau_2, freq)
    # n_3, k_3 = nk_from_eps(e_inf_3, e_s_3, tau_3, freq)
    # H_0 = H_sim(freq, d_air, n_1, k_1, d_1, n_2, k_2, d_2, n_3, k_3, d_3)
    H_0 = H_sim(f_ref, d_air, n1, k1, d_1, n2, k2, d_2, n3, k3, d_3)
    E_0 = irfft(H_0 * E_r_w)
    # E_0 = ifft(H_0 * E_r_w)
    return E_0 - E_meas


def fit_status(params, iter, resid, *args, **kws):

    d_air = params['d_air'].value
    # e_inf_1 = params['e_inf_1'].value
    # e_s_1 = params['e_s_1'].value
    # tau_1 = params['tau_1'].value
    d_1 = params['d_1'].value
    # e_inf_2 = params['e_inf_2'].value
    # e_s_2 = params['e_s_2'].value
    # tau_2 = params['tau_2'].value
    d_2 = params['d_2'].value
    # e_inf_3 = params['e_inf_3'].value
    # e_s_3 = params['e_s_3'].value
    # tau_3 = params['tau_3'].value
    d_3 = params['d_3'].value
    freq, E_r_w, E_meas, n1, k1, n2, k2, n3, k3 = args
    # n_1, k_1 = nk_from_eps(e_inf_1, e_s_1, tau_1, freq)
    # n_2, k_2 = nk_from_eps(e_inf_2, e_s_2, tau_2, freq)
    # n_3, k_3 = nk_from_eps(e_inf_3, e_s_3, tau_3, freq)
    H_0 = H_sim(freq, d_air, n1, k1, d_1, n2, k2, d_2, n3, k3, d_3)
    E_0 = irfft(H_0 * E_r_w)
    # E_0 = ifft(H_0 * E_r_w)
    J = sum((E_0 - E_meas)**2)
    # print(e_inf_1, e_s_1, tau_1, d_1, e_inf_2, e_s_2, tau_2, d_2, e_inf_3, e_s_3, tau_3, d_3)
    print('Iteration:', iter, 'J(x, theta) =', J)


# t_ref, E_ref = read_1file('./data/sim_resources/transmision_ref.txt')
t_ref, E_ref = read_1file('./data/sim_resources/refletion_ref.txt')
# t_ref, E_ref = read_1file('./data/sim_resources/mod_transmision_ref.txt')
# t_ref2, E_ref2 = read_1file('./data/sim_resources/mod_refletion_ref.txt')

t_ref *= 1e-12

f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
f_min, f_max = f_min_max_idx(f_ref, 0.2, 1.1)

f_ref *= 1e-12
n1 = 3 * ones(f_ref.size)
k1 = 0.00 * f_ref
n2 = 1.5 * ones(f_ref.size)
k2 = 0.00 * f_ref
n3 = 3 * ones(f_ref.size)
k3 = 0.00 * f_ref
f_ref *= 1e12

dair = 0  # 20e-8
d1 = 50e-6
d2 = 50e-6
d3 = 50e-6

H_0 = H_sim(f_ref, dair, n1, k1, d1, n2, k2, d2, n3, k3, d3)
figure(20)
plot(f_ref, abs(H_0))
xlim([0, 1e12])
figure(21)
plot(f_ref, angle(H_0))
xlim([0, 1e12])
E_sim = irfft(rfft(E_ref) * H_0)
# E_sim = ifft(fft(E_ref) * H_0)
# plot(t_ref, E_sim)
# noise = zeros(E_sim.size)
# for i in range(E_sim.size):
#     aux_snr = abs(E_sim[i]) / amax(abs(E_sim))
#     if aux_snr <= 0.1:
#         noise[i] = 0.01 * E_sim[i] * random.rand() / abs(E_sim[i])
# plot(t_ref, E_sim + noise)
# show()
# quit()
f_sim, E_sim_w = fourier_analysis(t_ref, E_sim)
H_w = E_sim_w / E_ref_w

# E_sim += noise


theta_params = lmfit.Parameters()
theta_params.add('d_air', 0, min=0, max=1e-4)
# theta_params.add('d_1', 50e-6, min=0, max=1e-3)
# theta_params.add('d_2', 12e-6, min=0, max=1e-3)
# theta_params.add('d_3', 5e-6, min=0, max=1e-3)
# for DE
# theta_params.add('d_1', 50e-6, min=1, max=1e-4)
# theta_params.add('d_2', 50e-6, min=1, max=1e-4)
# theta_params.add('d_3', 10e-6, min=1, max=1e-4)

# this s*** ain't fitting, they hatin'
# decentering = 20e-6
# theta_params.add('d_1', d1-decentering, min=0, max=1e-3)
# theta_params.add('d_2', d2-decentering, min=0, max=1e-3)
# theta_params.add('d_3', d3+decentering, min=0, max=1e-3)

theta_params.add('d_1', 50e-6, min=1e-6, max=1e-3)
theta_params.add('d_2', 50e-6, min=1e-6, max=1e-3)
theta_params.add('d_3', 50e-6, min=1e-6, max=1e-3)

max_err = 0  # -0.05  # in % fraction
n1 *= 1 + max_err
n2 *= 1 + max_err
n3 *= 1 + max_err

print('Fitting')


t1 = time_ns()
res = lmfit.minimize(resid, theta_params,
                     args=(f_ref, E_ref_w, E_sim, n1, k1, n2, k2, n3, k3),
                     # method='differential_evolution',
                     method='bfgs',
                     iter_cb=fit_status
                     )
t2 = time_ns()
theta_params_fit = res.params
d_air = theta_params_fit['d_air'].value
d_1 = theta_params_fit['d_1'].value
d_2 = theta_params_fit['d_2'].value
d_3 = theta_params_fit['d_3'].value

print('Results:')
print('White coat --- expected:', round(d1*1e6, 0), ' --- fitted:', round(d_1*1e6, 0),
      ' --- deviation:', abs(round(d1*1e6, 0) - round(d_1*1e6, 0)) / round(d1*1e6, 0))
print('Green coat --- d:', round(d2*1e6, 0), ' --- fitted:', round(d_2*1e6, 0),
      ' --- deviation:', abs(round(d2*1e6, 0) - round(d_2*1e6, 0)) / round(d2*1e6, 0))
print('Primer --- d:', round(d3*1e6, 0), ' --- fitted:', round(d_3*1e6, 0),
      ' --- deviation:', abs(round(d3*1e6, 0) - round(d_3*1e6, 0)) / round(d3*1e6, 0))
print_time_ns(t1, t2)


H_teo = H_sim(f_ref, d_air, n1, k1, d_1, n2, k2, d_2, n3, k3, d_3)
E_fit = irfft(H_teo*E_ref_w)
# E_fit = ifft(H_teo*E_ref_w)
t_ref *= 1e12

figure(1)
# plot(t_ref, E_ref, lw=1, label='ref')
plot(t_ref, E_sim, lw=1, label='sim')
# plot(t_ref, E_fit, lw=1, label='fit')
legend()
xlim([t_ref[0], t_ref[-1]])

# figure(2)
# plot(f_ref[f_min:f_max], abs(H_w)[f_min:f_max])
# plot(f_ref[f_min:f_max], abs(H_teo)[f_min:f_max])
#
# figure(3)
# plot(f_ref[f_min:f_max], unwrap(angle(H_w))[f_min:f_max])
# plot(f_ref[f_min:f_max], unwrap(angle(H_teo))[f_min:f_max])


show()
