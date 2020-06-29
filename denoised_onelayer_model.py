from TDSA import *
from scipy.optimize import curve_fit
from scipy.signal.windows import tukey
from time import time_ns
from genetic_denoising import genetic_deno

# constant definitions
deg_in = 30  # incidence angle in degrees
snell_sin = n_air * sin(deg_in * pi / 180)
layers = 3
# epsilon_model = 'debye'
# epsilon_model = 'cole'
n_subs = 1.17 - 0.0 * 1j  # substrate refractive index -- cork
# n_subs = 1.17e20 - 0.0 * 1j  # substrate refractive index -- metal
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
    return 4 * n_l * n_l_1 / (n_l + n_l_1) ** 2


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
    e_w = e_inf + (e_s - e_inf) / (1 + (1j * omg * tau) ** (alpha))
    return e_w


def nk_from_eps(e_s, e_inf, tau, freq):
    e_w = epsilon(e_s, e_inf, tau, freq)
    n = sqrt((abs(e_w) + real(e_w)) / 2)
    k = sqrt((abs(e_w) - real(e_w)) / 2)
    return n, k


def nk_from_eps_cc(e_s, e_inf, tau, alpha, freq):
    e_w = epsilon_cc(e_s, e_inf, tau, alpha, freq)
    n = sqrt((abs(e_w) + real(e_w)) / 2)
    k = sqrt((abs(e_w) - real(e_w)) / 2)
    return n, k


def H_sim(freq, d_air, n_1, k_1, d_1):
    k_1 *= freq * 1e-12

    H_i = cr_l_1_l(n_subs, n_1 - 1j * k_1) * ones(freq.size)
    rlm1l = cr_l_1_l(n_1 - 1j * k_1, n_air_cplx)
    tt = ct2(n_1 - 1j * k_1, n_air_cplx)
    exp_phi = phase_factor(n_1, k_1, d_1, freq)
    H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)

    return phase_factor(n_air, 0, d_air, freq) * H_i


# Main script
# Main script
# Boleto 176054
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_w_coat/ref metal wcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_w_coat/sam metal wcoat1_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_g_coat/ref metal gcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_g_coat/sam metal gcoat1_avg_f.txt')
t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/cork_w_coat/ref metal cork wcoat_avg_f.txt')
t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/cork_w_coat/sam cork wcoat1_avg_f.txt')

# Boleto 180881
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_180881_fecha_24_11_2017/metal_w_coat/ref metal wcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_180881_fecha_24_11_2017/metal_w_coat/sam metal wcoat 3_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_180881_fecha_24_11_2017/cork_w_coat/ref metal cork_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_180881_fecha_24_11_2017/cork_w_coat/sam cork wcoat 1_avg_f.txt')

# Boleto 177910
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_177910_fecha_04_12_2017/metal_w_coat/ref metal wcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_177910_fecha_04_12_2017/metal_w_coat/sam metal wcoat 1_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_177910_fecha_04_12_2017/cork_w_coat/ref metal cork_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_177910_fecha_04_12_2017/cork_w_coat/sam cork wcoat 1_avg_f.txt')

t_ref *= 1e-12
t_sam *= 1e-12

plot(t_sam, E_sam, label='org')
print('Denoising')
# E_sam = genetic_deno(t_ref, E_ref, t_sam, E_sam)

print('DSP')
delta_t_ref = mean(diff(t_ref))
ref_pulse_idx = centre_loc(E_ref)
window = tukey(2 * ref_pulse_idx)
window = zero_padding(window, 0, E_ref.size - window.size)
E_ref *= window
E_sam *= window
enlargement = 0 * E_ref.size
E_ref = zero_padding(E_ref, 0, enlargement)
t_ref = concatenate((t_ref, t_ref[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
E_sam = zero_padding(E_sam, 0, enlargement)
t_sam = concatenate((t_sam, t_sam[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
E_sam_max = amax(abs(E_sam))

window = zero_padding(window, 0, enlargement)


f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
delta_f_ref = mean(diff(f_ref))
# f_min, f_max = f_min_max_idx(f_ref, 0.2, 0.6)
# # f_min, f_max = f_min_max_idx(f_ref, 0, f_ref[-1]*1e-12)
# # f_min, f_max = f_min_max_idx(f_ref, 0, 10)
# f_ref = f_ref[f_min:f_max]
# E_ref_w = E_ref_w[f_min:f_max]
# f_sam = f_sam[f_min:f_max]
# E_sam_w = E_sam_w[f_min:f_max]

H_w = E_sam_w / E_ref_w
filt = wiener_filter(E_ref_w, beta=1e-2)
# E_sam_w *= filt
# H_w *= filt
# E_ref_w *= filt


# p0 = array((20e-6, 2.0, 0.1, 60e-6, 2.0, 0.1, 30e-6, 2.0, 0.1, 10e-6))
# p0 = array((20e-6, 3.0, 0.1, 60e-6, 3.0, 0.1, 20e-6, 3.0, 0.1, 5e-6))
# lwr_bnds = array((0, 1, 0, 1e-6, 1, 0, 1e-6, 1, 0, 1e-6))
# hgr_bnds = array((100e-6, 5, 10, 100e-6, 5, 10, 100e-6, 5, 10, 100e-6))

p0 = array((20e-6, 3, 0.1, 80e-6))
lwr_bnds = array((0, 2, 0, 80e-6))
hgr_bnds = array((100e-6, 5, 1, 81e-6))

# print(H_sim(f_ref, 20e-6, 2.0, 0.1, 50e-6, 2.0, 0.1, 50e-6, 2.0, 0.1, 50e-6))
# quit()

print('Fitting')


def E_sim(times, d_air, n_1, k_1, d_1):
    H_0 = H_sim(f_ref, d_air, n_1, k_1, d_1)
    E_0 = irfft(H_0 * E_ref_w)
    print('cost = ', sum(abs(E_0 - E_sam)))
    return E_0


t1 = time_ns()
popt, pcov = curve_fit(E_sim, t_sam, E_sam, p0,
                       bounds=(lwr_bnds, hgr_bnds),
                       method='dogbox'
                       )
t2 = time_ns()
print('popt =', popt)
print('pcov =', pcov)

print('Results:')
print('All-in-1 coat --- n:', round(popt[1], 2), 'k:', round(popt[2], 2), 'd:', round(popt[3] * 1e6, 0))

figure(1)
plot(t_sam, E_sam, label='dns')
plot(t_sam, E_sim(t_sam, popt[0], popt[1], popt[2], popt[3]), label='fit')
legend()

figure(2)
plot(f_ref, abs(H_w))
plot(f_ref, abs(H_sim(f_ref, popt[0], popt[1], popt[2], popt[3])))

figure(3)
plot(f_ref, unwrap(angle(H_w)))
plot(f_ref, unwrap(angle(H_sim(f_ref, popt[0], popt[1], popt[2], popt[3]))))
show()
