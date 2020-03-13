from TDSA import *
from scipy.optimize import differential_evolution, NonlinearConstraint, LinearConstraint
from scipy.signal.windows import tukey
from time import time_ns, strftime, gmtime
import math

# constant definitions
deg_in = 30  # incidence angle in degrees
snell_sin = n_air * sin(deg_in * pi / 180)
layers = 3
# epsilon_model = 'debye'
# epsilon_model = 'cole'
# n_subs = 1.17 - 0.0 * 1j  # substrate refractive index -- cork
n_subs = 1.17e20 - 0.0 * 1j  # substrate refractive index -- metal
# n_subs = 1.25 - 0.0 * 1j  # substrate refractive index -- cork 2.0


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


# def phase_factor(n, k, thick_n, thick_k, freq):  # theta in radians
#     omg = 2 * pi * freq
#     thick_n *= cos(theta(n))
#     phi_n = 2 * omg * thick_n / c_0
#     thick_k *= cos(theta(n))
#     phi_k = 2 * omg * thick_k / c_0
#     return exp(- 1j * phi_n) * exp(- k * phi_k)


def phase_factor(n, thick_n, thick_k, freq):  # theta in radians
    omg = 2 * pi * freq
    thick_n *= cos(theta(n))
    phi_n = 2 * omg * thick_n / c_0
    thick_k *= cos(theta(n))
    phi_k = 2 * omg * thick_k / c_0
    return exp(- 1j * phi_n) * exp(- phi_k)


def epsilon(e_s, e_inf, tau, freq):  # Debye model
    omg = 2 * pi * freq
    e_w = e_inf + (e_s - e_inf) / (1 + 1j * omg * tau)
    return e_w


def epsilon_cc(e_s, e_inf, tau, alpha, freq):  # Cole-Cole mode
    omg = 2 * pi * freq
    e_w = e_inf + (e_s - e_inf) / (1 + (1j * omg * tau)**(alpha))
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


def H_sim(ns, ks, thicks_n, thicks_k, d_air, freq):
    # most inner layer, in contact with substrate
    H_i = cr_l_1_l(n_subs, ns[-1] - 1j * ks[-1]) * ones(freq.size)
    rlm1l = cr_l_1_l(ns[-1] - 1j * ks[-1], ns[-2] - 1j * ks[-2])
    tt = ct2(ns[-1] - 1j * ks[-1], ns[-2] - 1j * ks[-2])
    # exp_phi = phase_factor(ns[-1], ks[-1], thicks_n[-1], thicks_k[-1], freq)
    exp_phi = phase_factor(ns[-1], thicks_n[-1], thicks_k[-1], freq)
    H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)
    
    # inner layers, no contact with air nor substrate, indexes from 'layer_qty - 2' to '1'
    for i in range(1, layers - 1):
        i = layers - i - 1
        rlm1l = cr_l_1_l(ns[i] - 1j * ks[i], ns[i - 1] - 1j * ks[i - 1])
        tt = ct2(ns[i] - 1j * ks[i], ns[i - 1] - 1j * ks[i - 1])
        # exp_phi = phase_factor(ns[i], ks[i], thicks_n[i], thicks_k[i], freq)
        exp_phi = phase_factor(ns[i], thicks_n[i], thicks_k[i], freq)
        H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)
    
    # most outer layer, in contact with air
    rlm1l = cr_l_1_l(ns[0] - 1j * ks[0], n_air_cplx)
    tt = ct2(ns[0] - 1j * ks[0], n_air_cplx)
    # exp_phi = phase_factor(ns[0], ks[0], thicks_n[0], thicks_k[0], freq)
    exp_phi = phase_factor(ns[0], thicks_n[0], thicks_k[0], freq)
    H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)
    return exp(- 1j * 2 * 2 * pi * freq * d_air / c_0) * H_i


def cost_function(params, *args):
    d_air = params[0]
    params = params[1:]
    ns = list()
    ks = list()
    thicks_n = list()
    thicks_k = list()
    E_sam, E_ref_w, freqs = args
    for i in range(layers):
        e_inf = params[5 * i]
        e_s = params[5 * i + 1]
        tau = params[5 * i + 2]
        thicks_n.append(params[5 * i + 3])
        thicks_k.append(params[5 * i + 4])
        n, k = nk_from_eps(e_inf, e_s, tau, freqs)  # debye model
        ns.append(n)
        ks.append(k)
    ns = array(ns)
    ks = array(ks)
    thicks_n = array(thicks_n)
    thicks_k = array(thicks_k)
    H_teo = H_sim(ns, ks, thicks_n, thicks_k, d_air, freqs)
    E_sam_teo_w = E_ref_w * H_teo
    E_sam_teo = irfft(E_sam_teo_w, n=E_sam.size)
    return sum((E_sam_teo - E_sam)**2)  # + sum(params**2)


# Main script
# Boleto 176054
t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_w_coat/ref metal wcoat_avg_f.txt')
t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_w_coat/sam metal wcoat1_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_g_coat/ref metal gcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_g_coat/sam metal gcoat1_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/cork_w_coat/ref metal cork wcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/cork_w_coat/sam cork wcoat2_avg_f.txt')

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

t_ref *= 1e-12
t_sam *= 1e-12


f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
# f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)


alpha = 1
beta = 100
k_bounds = [  # calibration
    (0, 100e-6),  # air thickness
    (alpha, beta), (alpha, beta), (1e-14, 1e-12), (0, 100e-6), (0, 100e-6),
    (alpha, beta), (alpha, beta), (1e-14, 1e-12), (0, 100e-6), (0, 100e-6),
    (alpha, beta), (alpha, beta), (1e-14, 1e-12), (0, 100e-6), (0, 100e-6)
]

print('Fitting')
t1 = time_ns()
res = differential_evolution(cost_function,
                             k_bounds,
                             args=(E_sam, E_ref_w, f_ref),
                             # popsize=150,
                             # maxiter=2000,
                             disp=True,  # step cost_function value
                             polish=True
                             )
t2 = time_ns()
print()
d_air = res.x[0]
print(res)
res.x = res.x[1:]
# res.x[2] /= res.x[0]
# res.x[3] /= res.x[1]
# res.x[6] /= res.x[4]
# res.x[7] /= res.x[5]
print('Results:')
print('White coat -', 'n:', round(res.x[0], 2), 'k:', round(res.x[1], 2), 'd_n:', round(res.x[2] * 1e6, 2), 'um', 'd_k:', round(res.x[3] * 1e6, 2), 'um')
print('Green coat -', 'n:', round(res.x[4], 2), 'k:', round(res.x[5], 2), 'd_n:', round(res.x[6] * 1e6, 2), 'um', 'd_k:', round(res.x[7] * 1e6, 2), 'um')
if layers == 3:
    # res.x[]
    print('Primer coat -', 'n:', round(res.x[8], 2), 'k:', round(res.x[9], 2), 'd:', round(res.x[10] * 1e6, 2), 'um', 'd_k:', round(res.x[11] * 1e6, 2), 'um')
    print('Total:', round((res.x[2] + res.x[6] + res.x[10]) * 1e6, 2), 'um')
    # ns = array([res.x[0], res.x[3], res.x[6]])
    # ks = array([res.x[1], res.x[4], res.x[7]])
    # thikcs = array([res.x[2], res.x[5], res.x[8]])
else:
    print('Total:', round((res.x[2] + res.x[6]) * 1e6, 2), 'um')
    # ns = array([res.x[0], res.x[3]])
    # ks = array([res.x[1], res.x[4]])
    # thikcs = array([res.x[2], res.x[5]])
print()
secs = (t2-t1)*1e-9
if secs < 3600:
    print('Processing time (mm:ss):', strftime('%M:%S', gmtime(secs)))
else:
    print('Processing time (mm:ss):', strftime('%H:%M:%S', gmtime(secs)))
print()