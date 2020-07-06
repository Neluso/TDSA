from TDSA import *
from scipy.optimize import differential_evolution
from scipy.signal.windows import tukey
from time import time_ns, strftime, gmtime

# constant definitions
deg_in = 30  # incidence angle in degrees
snell_sin = n_air * sin(deg_in * pi / 180)
layers = 3
# epsilon_model = 'debye'
# epsilon_model = 'cole'
n_subs = 1.17 - 0.0 * 1j  # substrate refractive index -- cork
# n_subs = 1.17e20 - 0.0 * 1j  # substrate refractive index -- metal
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


def phase_factor(n, k, thick, freq):  # theta in radians
    omg = 2 * pi * freq
    thick *= cos(theta(n))
    phi_n = 2 * omg * thick / c_0
    phi_k = 2 * omg * thick / c_0
    return exp(- 1j * phi_n) * exp(- k * phi_k)


def epsilon(e_s, e_inf, tau, freq):  # Debye model
    omg = 2 * pi * freq
    e_w = e_inf + (e_s - e_inf) / (1 + 1j * omg * tau)
    return e_w


def nk_from_eps(e_s, e_inf, tau, freq):
    e_w = epsilon(e_s, e_inf, tau, freq)
    n = sqrt((abs(e_w) + real(e_w)) / 2)
    k = sqrt((abs(e_w) - real(e_w)) / 2)
    return n, k


def H_sim(ns, ks, thicks, d_air, freq):
    # most inner layer, in contact with substrate
    H_i = cr_l_1_l(n_subs, ns[-1] - 1j * ks[-1]) * ones(freq.size)
    rlm1l = cr_l_1_l(ns[-1] - 1j * ks[-1], ns[-2] - 1j * ks[-2])
    tt = ct2(ns[-1] - 1j * ks[-1], ns[-2] - 1j * ks[-2])
    exp_phi = phase_factor(ns[-1], ks[-1], thicks[-1], freq)
    H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)
    
    # inner layers, no contact with air nor substrate, indexes from 'layer_qty - 2' to '1'
    for i in range(1, layers - 1):
        i = layers - i - 1
        rlm1l = cr_l_1_l(ns[i] - 1j * ks[i], ns[i - 1] - 1j * ks[i - 1])
        tt = ct2(ns[i] - 1j * ks[i], ns[i - 1] - 1j * ks[i - 1])
        exp_phi = phase_factor(ns[i], ks[i], thicks[i], freq)
        H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)
    
    # most outer layer, in contact with air
    rlm1l = cr_l_1_l(ns[0] - 1j * ks[0], n_air_cplx)
    tt = ct2(ns[0] - 1j * ks[0], n_air_cplx)
    exp_phi = phase_factor(ns[0], ks[0], thicks[0], freq)
    H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)
    return exp(- 1j * 2 * 2 * pi * freq * d_air / c_0) * H_i


def cost_function(params, *args):
    d_air = params[0]
    params = params[1:]
    ns = list()
    ks = list()
    thicks = list()
    E_sam, E_ref_w, freqs = args
    for i in range(layers):
        e_inf = params[4 * i]
        e_s = params[4 * i + 1]
        tau = params[4 * i + 2]
        thicks.append(params[4 * i + 3])
        n, k = nk_from_eps(e_inf, e_s, tau, freqs)  # debye model
        ns.append(n)
        ks.append(k)
    ns = array(ns)
    ks = array(ks)
    thicks = array(thicks)
    H_teo = H_sim(ns, ks, thicks, d_air, freqs)
    E_sam_teo_w = E_ref_w * H_teo
    E_sam_teo = irfft(E_sam_teo_w, n=E_sam.size)
    aus, E_sam_w = fourier_analysis(freqs, E_sam)
    f_min_idx, f_max_idx = f_min_max_idx(freqs, 0.1, 1.1)
    H_w = E_sam_w / E_ref_w
    H_w = H_w[f_min_idx:f_max_idx]
    H_teo = H_teo[f_min_idx:f_max_idx]
    return sum((E_sam_teo - E_sam)**2) + sum((abs(H_teo) - abs(H_w))**2) + sum((unwrap(angle(H_teo)) - unwrap(angle(H_w)))**2)


# Main script
# Boleto 176054
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_w_coat/ref metal wcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_w_coat/sam metal wcoat3_avg_f.txt')
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_g_coat/ref metal gcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/metal_g_coat/sam metal gcoat1_avg_f.txt')
t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/cork_w_coat/ref metal cork wcoat_avg_f.txt')
t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_176054_fecha_15_06_2018/cork_w_coat/sam cork wcoat1_avg_f.txt')

# Boleto 180881
# t_ref, E_ref = read_1file('./data/muestras_airbus_boleto_180881_fecha_24_11_2017/metal_w_coat/ref metal wcoat_avg_f.txt')
# t_sam, E_sam = read_1file('./data/muestras_airbus_boleto_180881_fecha_24_11_2017/metal_w_coat/sam metal wcoat 1_avg_f.txt')
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
# E_ref *= window
# E_sam *= window
enlargement = 0 * E_ref.size
E_ref = zero_padding(E_ref, 0, enlargement)
t_ref = concatenate((t_ref, t_ref[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
E_sam = zero_padding(E_sam, 0, enlargement)
t_sam = concatenate((t_sam, t_sam[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
E_sam_max = amax(abs(E_sam))

t_ref *= 1e-12
t_sam *= 1e-12


f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)

H_w = E_sam_w / E_ref_w


# alpha = 1
# beta = 100
# k_bounds = [  # calibration
#     (-100e-6, 100e-6),  # air thickness
#     (alpha, beta), (alpha, beta), (1e-14, 1e-12), (10e-6, 100e-6),
#     (alpha, beta), (alpha, beta), (1e-14, 1e-12), (10e-6, 100e-6),
#     (alpha, beta), (alpha, beta), (1e-14, 1e-12), (1e-6, 20e-6)
# ]

alpha = 10
beta = 15
# k_bounds = [  # calibration
#     (-100e-6, 100e-6),  # air thickness
#     (alpha, beta), (alpha, beta), (1e-15, 2e-12), (1e-6, 100e-6),
#     (alpha, beta), (alpha, beta), (1e-15, 2e-12), (1e-6, 100e-6),
#     (alpha, beta), (alpha, beta), (1e-15, 2e-12), (1e-6, 100e-6)
# ]

k_bounds = [  # calibration
    (-100e-6, 100e-6),  # air thickness
    (alpha, beta), (alpha, beta), (0, 1e-13), (1e-6, 100e-6),
    (alpha, beta), (alpha, beta), (0, 1e-13), (1e-6, 100e-6),
    (alpha, beta), (alpha, beta), (0, 1e-13), (1e-6, 100e-6)
]

print('Fitting')
t1 = time_ns()
res = differential_evolution(cost_function,
                             k_bounds,
                             args=(E_sam, E_ref_w, f_ref),
                             # popsize=45,
                             # maxiter=2000,
                             disp=True,  # step cost_function value
                             polish=True
                             )
t2 = time_ns()
print()
d_air = res.x[0]
print(res)
n1, k1 = nk_from_eps(res.x[1], res.x[2], res.x[3], f_ref)
print('White coat -', 'n:', round(mean(n1), 2), 'k:', round(mean(k1), 2), 'd:', round(res.x[4] * 1e6, 0), 'um')
print('\t\t e_s:', round(res.x[2], 2), 'e_inf', round(res.x[1], 2), 'tau', res.x[3])
n2, k2 = nk_from_eps(res.x[5], res.x[6], res.x[7], f_ref)
print('Green coat -', 'n:', round(mean(n2), 2), 'k:', round(mean(k2), 2), 'd:', round(res.x[8] * 1e6, 0), 'um')
print('\t\t e_s:', res.x[6], 'e_inf:', res.x[5], 'tau:', res.x[7])
n3, k3 = nk_from_eps(res.x[9], res.x[10], res.x[11], f_ref)
print('Primer coat -', 'n:', round(mean(n3), 2), 'k:', round(mean(k3), 2), 'd:', round(res.x[12] * 1e6, 0), 'um')
print('\t\t e_s:', res.x[10], 'e_inf:', res.x[9], 'tau:', res.x[11])
print('Total:', round((res.x[4] + res.x[8] + res.x[12]) * 1e6, 0), 'um')
print()
print('Air - d: ', round(d_air * 1e6, 2), 'um')
secs = (t2-t1)*1e-9
if secs < 3600:
    print('Processing time (mm:ss):', strftime('%M:%S', gmtime(secs)))
else:
    print('Processing time (hh:mm:ss):', strftime('%H:%M:%S', gmtime(secs)))
print()


ns = list()
ks = list()
thicks = list()
params = res.x[1:]
for i in range(layers):
    e_inf = params[4 * i]
    e_s = params[4 * i + 1]
    tau = params[4 * i + 2]
    thicks.append(params[4 * i + 3])
    n, k = nk_from_eps(e_inf, e_s, tau, f_ref)  # debye model
    ns.append(n)
    ks.append(k)
ns = array(ns)
ks = array(ks)
thicks = array(thicks)

H_fit = H_sim(ns, ks, thicks, d_air, f_ref)
E_fit = irfft(H_fit * E_ref_w)

f_min_idx, f_max_idx = f_min_max_idx(f_ref, 0.05, 1.5)


t_sam *= 1e12
t_ref *= 1e12
f_ref *= 1e-12
f_sam *= 1e-12

f_ref = f_ref[f_min_idx:f_max_idx]
H_w = H_w[f_min_idx:f_max_idx]
H_fit = H_fit[f_min_idx:f_max_idx]
n1 = n1[f_min_idx:f_max_idx]
k1 = k1[f_min_idx:f_max_idx]
n2 = n2[f_min_idx:f_max_idx]
k2 = k2[f_min_idx:f_max_idx]
n3 = n3[f_min_idx:f_max_idx]
k3 = k3[f_min_idx:f_max_idx]


figure(1)
plot(t_sam, E_sam, label='sam')
plot(t_sam, E_fit, label='fit')
xlim([t_sam[0], t_sam[-1]])
xlabel('t (ps)')
legend()

figure(2)
plot(f_ref, abs(H_w), label='sam')
plot(f_ref, abs(H_fit), label='fit')
xlim([f_ref[0], f_ref[-1]])
ylabel('Abs')
xlabel('f (THz)')
legend()

figure(3)
plot(f_ref[f_min_idx:f_max_idx], unwrap(angle(H_w))[f_min_idx:f_max_idx], label='sam')
plot(f_ref[f_min_idx:f_max_idx], unwrap(angle(H_fit))[f_min_idx:f_max_idx], label='fit')
xlim([f_ref[0], f_ref[-1]])
ylabel('Arg')
xlabel('f (THz)')
legend()

figure(4)
plot(f_ref, n1, label='w')
plot(f_ref, n2, label='g')
if layers == 3:
    plot(f_ref, n3, label='p')
title('n')
xlim([f_ref[0], f_ref[-1]])
ylabel('n')
xlabel('f (THz)')
legend()

figure(5)
plot(f_ref, k1, label='w')
plot(f_ref, k2, label='g')
if layers == 3:
    plot(f_ref, k3, label='p')
title('k')
xlim([f_ref[0], f_ref[-1]])
ylabel('k')
xlabel('f (THz)')
legend()

show()
