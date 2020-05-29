from TDSA import *
from scipy.optimize import differential_evolution, NonlinearConstraint, LinearConstraint, curve_fit
from scipy.signal.windows import tukey
from time import time_ns
import lmfit
import os


# constant definitions
deg_in = 30  # incidence angle in degrees
snell_sin = n_air * sin(deg_in * pi / 180)
layers = 3
# epsilon_model = 'debye'
# epsilon_model = 'cole'
# n_subs = 1.17 - 0.0 * 1j  # substrate refractive index -- cork
n_subs = 1.17e20 - 0.0 * 1j  # substrate refractive index -- metal
# n_subs = 1.25 - 0.0 * 1j  # substrate refractive index -- cork 2.0


reg_rel = 0.1
reg_matrix = identity(1 + layers * 4)
for i in range(reg_matrix.shape[0]):
    for j in range(reg_matrix.shape[1]):
        if i == j:
            reg_matrix[i, j] = 1
        elif i == 0:
            reg_matrix[i, j] = 0
        elif j == 0:
            reg_matrix[i, j] = 0
        else:
            reg_matrix[i, j] = reg_rel
reg_param = 0.1
reg_matrix *= reg_param



# debug variables
debug_value_1 = list()
debug_value_2 = list()
debug_value_3 = list()


# function definitions
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


def nk_from_eps_cc(e_s, e_inf, tau, alpha, freq):
    e_w = epsilon_cc(e_s, e_inf, tau, alpha, freq)
    n = sqrt((abs(e_w) + real(e_w)) / 2)
    k = sqrt((abs(e_w) - real(e_w)) / 2)
    return n, k


def H_sim(freq, d_air, e_inf_1, e_s_1, tau_1, d_1, e_inf_2, e_s_2, tau_2, d_2, e_inf_3, e_s_3, tau_3, d_3):
    n_1, k_1 = nk_from_eps(e_inf_1, e_s_1, tau_1, freq)
    n_2, k_2 = nk_from_eps(e_inf_2, e_s_2, tau_2, freq)
    n_3, k_3 = nk_from_eps(e_inf_3, e_s_3, tau_3, freq)

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
    e_inf_1 = params['e_inf_1'].value
    e_s_1 = params['e_s_1'].value
    tau_1 = params['tau_1'].value
    d_1 = params['d_1'].value
    e_inf_2 = params['e_inf_2'].value
    e_s_2 = params['e_s_2'].value
    tau_2 = params['tau_2'].value
    d_2 = params['d_2'].value
    e_inf_3 = params['e_inf_3'].value
    e_s_3 = params['e_s_3'].value
    tau_3 = params['tau_3'].value
    d_3 = params['d_3'].value
    freq, E_r_w, E_meas, H_w = args
    H_0 = H_sim(freq, d_air, e_inf_1, e_s_1, tau_1, d_1, e_inf_2, e_s_2, tau_2, d_2, e_inf_3, e_s_3, tau_3, d_3)
    E_0 = ifft(H_0 * E_r_w)
    return E_0 - E_meas + abs(H_0) - abs(H_w) + unwrap(angle(H_0)) - unwrap(angle(H_w))


def fit_status(params, iter, resid, *args, **kws):
    d_air = params['d_air'].value
    e_inf_1 = params['e_inf_1'].value
    e_s_1 = params['e_s_1'].value
    tau_1 = params['tau_1'].value
    d_1 = params['d_1'].value
    e_inf_2 = params['e_inf_2'].value
    e_s_2 = params['e_s_2'].value
    tau_2 = params['tau_2'].value
    d_2 = params['d_2'].value
    e_inf_3 = params['e_inf_3'].value
    e_s_3 = params['e_s_3'].value
    tau_3 = params['tau_3'].value
    d_3 = params['d_3'].value
    freq, E_r_w, E_meas, H_w = args
    H_0 = H_sim(freq, d_air, e_inf_1, e_s_1, tau_1, d_1, e_inf_2, e_s_2, tau_2, d_2, e_inf_3, e_s_3, tau_3, d_3)
    E_0 = ifft(H_0 * E_r_w)
    J = sum((E_0 - E_meas)**2)
    # print(e_inf_1, e_s_1, tau_1, d_1, e_inf_2, e_s_2, tau_2, d_2, e_inf_3, e_s_3, tau_3, d_3)
    print('Iteration:', iter, 'J(x, theta) =', J)


# fs = arange(100)/10
# e_inf = 10
# e_s = 20
# tau = 0.01
# n_test, k_test = nk_from_eps(e_inf, e_s, tau, fs)
# plot(fs, n_test)
# plot(fs, sqrt(e_s) * ones(fs.size))
# show()
# quit()


# Main script
t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_w_coat/ref metal wcoat_avg_f.txt')
t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_w_coat/sam metal wcoat3_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_g_coat/ref metal gcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_g_coat/sam metal gcoat1_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/cork_w_coat/ref metal cork wcoat_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/cork_w_coat/ref cork wcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/cork_w_coat/sam cork wcoat1_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_180881_fecha_24_11_2017/metal_w_coat/ref metal wcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_180881_fecha_24_11_2017/metal_w_coat/sam metal wcoat 1_avg_f.txt')

print('DSP')
delta_t_ref = mean(diff(t_ref))
ref_pulse_idx = centre_loc(E_ref)
window = tukey(2 * ref_pulse_idx)
window = zero_padding(window, 0, E_ref.size - window.size)
# E_ref *= window
# E_sam *= window
# enlargement = 0 * E_ref.size
# E_ref = zero_padding(E_ref, 0, enlargement)
# t_ref = concatenate((t_ref, t_ref[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
# E_sam = zero_padding(E_sam, 0, enlargement)
# t_sam = concatenate((t_sam, t_sam[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
# E_sam_max = amax(abs(E_sam))
#
# window = zero_padding(window, 0, enlargement)

t_ref *= 1e-12
t_sam *= 1e-12


f_ref, E_ref_w = fourier_analysis_comp(t_ref, E_ref)
f_sam, E_sam_w = fourier_analysis_comp(t_sam, E_sam)
delta_f_ref = mean(diff(f_ref))
f_min, f_max = f_min_max_idx(f_ref, 0.2, 1.1)
# # f_min, f_max = f_min_max_idx(f_ref, 0, f_ref[-1]*1e-12)
# # f_min, f_max = f_min_max_idx(f_ref, 0, 10)
# f_ref = f_ref[f_min:f_max]
# E_ref_w = E_ref_w[f_min:f_max]
# f_sam = f_sam[f_min:f_max]
# E_sam_w = E_sam_w[f_min:f_max]

H_w = E_sam_w / E_ref_w
filt = wiener_filter(E_ref_w, beta=2e-3)
# plot(f_ref, toDb(wiener_filter(E_ref_w, beta=1e-1)), label='1')
# plot(f_ref, toDb(wiener_filter(E_ref_w, beta=1e-2)), label='2')
# plot(f_ref, toDb(wiener_filter(E_ref_w, beta=1e-3)), label='3')
# plot(f_ref, toDb(wiener_filter(E_ref_w, beta=2e-4)), label='3.5')
# plot(f_ref, toDb(wiener_filter(E_ref_w, beta=1e-4)), label='4')
# legend()
# show()
# quit()

# E_sam_w *= filt
# H_w *= filt
# E_ref_w *= filt


theta_params = lmfit.Parameters()
theta_params.add('d_air', 20e-6, min=0, max=1e-4)

eps_max = 14  # 14
eps_min = 5  # 5
discrp = 3

theta_params.add('discr1', 0, min=-discrp, max=discrp)
theta_params.add('discr2', 0, min=-discrp, max=discrp)
theta_params.add('discr3', 0, min=-discrp, max=discrp)
theta_params.add('discr4', 0, min=-discrp, max=discrp)

theta_params.add('e_inf_1', 9, min=eps_min, max=eps_max)
theta_params.add('e_s_1', 9, min=eps_min, max=eps_max)
theta_params.add('tau_1', 1e-12, min=1e-15, max=1e-9)

theta_params.add('e_inf_2', 9, min=eps_min, max=eps_max, expr='discr1 + e_inf_1')
theta_params.add('e_s_2', 9, min=eps_min, max=eps_max, expr='discr2 + e_s_1')
theta_params.add('tau_2', 1e-12, min=1e-15, max=1e-9)

theta_params.add('e_inf_3', 9, min=eps_min, max=eps_max, expr='discr3 + e_inf_1')
theta_params.add('e_s_3', 9, min=eps_min, max=eps_max, expr='discr4 + e_s_1')
theta_params.add('tau_3', 1e-12, min=1e-15, max=1e-9)

# theta_params.add('d_1', 80e-6, min=80e-6, max=81e-6)
# theta_params.add('d_2', 20e-6, min=19e-6, max=20e-6)
# theta_params.add('d_3', 5e-6, min=4e-6, max=5e-6)


# # pseudo_calibration
# theta_params.add('e_inf_1', 7.24, min=7, max=8)
# theta_params.add('e_s_1', 12.2, min=12, max=13)
# theta_params.add('tau_1', 7.23e-13, min=7e-13, max=8.1e-13)
#
# theta_params.add('e_inf_2', 5, min=5, max=6)
# theta_params.add('e_s_2', 14, min=14, max=15)
# theta_params.add('tau_2', 2.08e-12, min=2e-12, max=3e-12)
#
# theta_params.add('e_inf_3', 11.24, min=11, max=12)
# theta_params.add('e_s_3', 8.2, min=8, max=9)
# theta_params.add('tau_3', 1.08e-12, min=1e-12, max=2e-12)
# #
theta_params.add('d_1', 80e-6, min=0)  # , min=80e-6, max=81e-6)
theta_params.add('d_2', 20e-6, min=0)  # , min=19e-6, max=20e-6)
theta_params.add('d_3', 5e-6, min=0)  # , min=4e-6, max=5e-6)


print('Fitting')


t1 = time_ns()
res = lmfit.minimize(resid, theta_params,
                     args=(f_ref, E_ref_w, E_sam, H_w),
                     # method='differential_evolution',
                     method='bfgs',
                     iter_cb=fit_status
                     )
t2 = time_ns()
theta_params_fit = res.params
d_air = theta_params_fit['d_air'].value
e_inf_1 = theta_params_fit['e_inf_1'].value
e_s_1 = theta_params_fit['e_s_1'].value
tau_1 = theta_params_fit['tau_1'].value
d_1 = theta_params_fit['d_1'].value
e_inf_2 = theta_params_fit['e_inf_2'].value
e_s_2 = theta_params_fit['e_s_2'].value
tau_2 = theta_params_fit['tau_2'].value
d_2 = theta_params_fit['d_2'].value
e_inf_3 = theta_params_fit['e_inf_3'].value
e_s_3 = theta_params_fit['e_s_3'].value
tau_3 = theta_params_fit['tau_3'].value
d_3 = theta_params_fit['d_3'].value

n_1, k_1 = nk_from_eps(e_inf_1, e_s_1, tau_1, f_ref)
n_2, k_2 = nk_from_eps(e_inf_2, e_s_2, tau_2, f_ref)
n_3, k_3 = nk_from_eps(e_inf_3, e_s_3, tau_3, f_ref)


# plot(f_ref, n_1, f_ref, n_2, f_ref, n_3)
# show()
# quit()

print('Results:')
print('White coat --- n:', round(mean(n_1), 2), 'k:', round(mean(k_1), 3), 'd:', round(d_1*1e6, 0))
print('e_inf:', e_inf_1, ', e_s:', e_s_1, ', tau:', tau_1)
print('Green coat --- n:', round(mean(n_2), 2), 'k:', round(mean(k_2), 3), 'd:', round(d_2*1e6, 0))
print('e_inf:', e_inf_2, ', e_s:', e_s_2, ', tau:', tau_2)
print('Primer --- n:', round(mean(n_3), 2), 'k:', round(mean(k_3), 3), 'd:', round(d_3*1e6, 0))
print('e_inf:', e_inf_3, ', e_s:', e_s_3, ', tau:', tau_3)
print_time_ns(t1, t2)

# print(str(e_inf_1) + ';' + str(e_inf_2) + ';' + str(e_inf_3) + ';' + str(e_s_1) + ';' + str(e_s_2) + ';' + str(e_s_3) + ';' + str(tau_1) + ';' + str(tau_2) + ';' + str(tau_3))

H_teo = H_sim(f_ref, d_air, e_inf_1, e_s_1, tau_1, d_1, e_inf_2, e_s_2, tau_2, d_2, e_inf_3, e_s_3, tau_3, d_3)

figure(1)
plot(t_sam, E_sam)
plot(t_sam, irfft(H_teo*E_ref_w))

figure(2)
plot(f_ref[f_min:f_max], abs(H_w)[f_min:f_max])
plot(f_ref[f_min:f_max], abs(H_teo)[f_min:f_max])

figure(3)
plot(f_ref[f_min:f_max], unwrap(angle(H_w))[f_min:f_max])
plot(f_ref[f_min:f_max], unwrap(angle(H_teo))[f_min:f_max])

figure(4)
plot(f_ref[f_min:f_max], n_1[f_min:f_max], f_ref[f_min:f_max], n_2[f_min:f_max], f_ref[f_min:f_max], n_3[f_min:f_max])
figure(5)
plot(f_ref[f_min:f_max], k_1[f_min:f_max], f_ref[f_min:f_max], k_2[f_min:f_max], f_ref[f_min:f_max], k_3[f_min:f_max])


# print(theta_params_fit)


show()
