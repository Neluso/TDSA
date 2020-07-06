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
# n_subs = 1.17 - 0.0 * 1j  # substrate refractive index -- cork
n_subs = 1.17e20 - 0.0 * 1j  # substrate refractive index -- metal
# n_subs = 1.25 - 0.0 * 1j  # substrate refractive index -- cork 2.0


# function definitions
def theta(n):
    return arcsin(snell_sin / real(n))


def ct(ni, nt):
    ni *= cos(theta(ni))
    nt *= cos(theta(nt))
    return 2 * ni / (ni + nt)


def cr(ni, nt):  # from n_l-1 to n_l
    ni *= cos(theta(ni))
    nt *= cos(theta(nt))
    return (ni - nt) / (ni + nt)


def phase_factor(n, k, thick, freq):  # theta in radians
    omg = 2 * pi * freq
    thick *= cos(theta(n))
    phi_n = omg * thick / c_0
    phi_k = omg * thick / c_0
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


# def H_sim(ns, ks, thicks, d_air, freq):
#     # most inner layer, in contact with substrate
#     H_i = cr_l_1_l(n_subs, ns[-1] - 1j * ks[-1]) * ones(freq.size)
#     rlm1l = cr_l_1_l(ns[-1] - 1j * ks[-1], ns[-2] - 1j * ks[-2])
#     tt = ct2(ns[-1] - 1j * ks[-1], ns[-2] - 1j * ks[-2])
#     exp_phi = phase_factor(ns[-1], ks[-1], thicks[-1], freq)
#     H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)
#
#     # inner layers, no contact with air nor substrate, indexes from 'layer_qty - 2' to '1'
#     for i in range(1, layers - 1):
#         i = layers - i - 1
#         rlm1l = cr_l_1_l(ns[i] - 1j * ks[i], ns[i - 1] - 1j * ks[i - 1])
#         tt = ct2(ns[i] - 1j * ks[i], ns[i - 1] - 1j * ks[i - 1])
#         exp_phi = phase_factor(ns[i], ks[i], thicks[i], freq)
#         H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)
#
#     # most outer layer, in contact with air
#     rlm1l = cr_l_1_l(ns[0] - 1j * ks[0], n_air_cplx)
#     tt = ct2(ns[0] - 1j * ks[0], n_air_cplx)
#     exp_phi = phase_factor(ns[0], ks[0], thicks[0], freq)
#     H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)
#     return exp(- 1j * 2 * 2 * pi * freq * d_air / c_0) * H_i


def H_sim(ns, ks, thicks, d_air, freq):
    
    # First layer through
    H_t = list()
    ta0 = ct(n_air_cplx, ns[0] - 1j * ks[0])
    t01 = ct(ns[0] - 1j * ks[0], ns[1] - 1j * ks[1])
    r01 = cr(ns[0] - 1j * ks[0], ns[1] - 1j * ks[1])
    r0a = cr(ns[0] - 1j * ks[0], n_air_cplx)
    exp_phi = phase_factor(ns[0], ks[0], thicks[0], freq)
    H_t.append(phase_factor(n_air, 0, d_air, freq) * ta0 * exp_phi * t01 / (1 - r01 * r0a * exp_phi))
    
    # Mid layers through
    for i in range(1, layers - 1):
        tlf = ct(ns[i] - 1j * ks[i], ns[i+1] - 1j * ks[i+1])  # l = layer, f = following
        rlf = cr(ns[i] - 1j * ks[i], ns[i+1] - 1j * ks[i+1])
        rlp = cr(ns[i] - 1j * ks[i], ns[i-1] - 1j * ks[i-1])  # p = previous
        exp_phi = phase_factor(ns[0], ks[0], thicks[0], freq)
        H_t.append(H_t[-1] * exp_phi * tlf / (1 - rlf * rlp * exp_phi))

    # Inner layer only reflects
    rpl = cr(ns[-2] - 1j * ns[-2], ns[-1] - 1j * ns[-1])
    tlp = ct(ns[-1] - 1j * ns[-1], ns[-2] - 1j * ns[-2])
    rls = cr(ns[-1] - 1j * ns[-1], n_subs)  # s = substrate
    exp_phi = phase_factor(ns[-1], ks[-1], 2 * thicks[-1], freq)
    H_r = rpl + H_t[-1] * tlp * rls * exp_phi / (1 + rpl * rls * exp_phi)
    
    # Mid layers reflection
    for i in range(1, layers - 1):
        j = layers - i - 1
        rpl = cr(ns[j-1] - 1j * ns[j-1], ns[j] - 1j * ns[j])
        tlp = ct(ns[j] - 1j * ns[j], ns[j-1] - 1j * ns[j-1])
        rlf = cr(ns[j] - 1j * ns[j], ns[j+1] - 1j * ns[j+1])
        exp_phi = phase_factor(ns[j], ks[j], 2 * thicks[j], freq)
        H_r = rpl + H_t[j-1] * tlp * H_r * exp_phi / (1 + rpl * H_r * exp_phi)
        
    # First layer reflection
    ta0 = ct(n_air_cplx, ns[0] - 1j * ks[0])
    t0a = ct(ns[0] - 1j * ks[0], n_air_cplx)
    r0a = cr(ns[0] - 1j * ks[0], n_air_cplx)
    exp_phi = phase_factor(ns[0], ks[0], 2 * thicks[0], freq)
    H_r = r0a + ta0 * t0a * H_r * exp_phi / (1 + r0a * H_r * exp_phi)
    
    return H_r * phase_factor(n_air, 0, d_air, freq)


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
    return sum((E_sam_teo - E_sam)**2)  # + sum(params**2)


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
#     (-200e-6, 200e-6),  # air thickness
#     (alpha, beta), (alpha, beta), (0, 1), (1e-6, 1000e-6),
#     (alpha, beta), (alpha, beta), (0, 1), (1e-6, 1000e-6),
#     (alpha, beta), (alpha, beta), (0, 1), (1e-6, 1000e-6)
# ]

alpha = 10
beta = 12
k_bounds = [  # calibration
    (-100e-6, 100e-6),  # air thickness
    (alpha, beta), (alpha, beta), (0, 1), (1e-6, 100e-6),
    (alpha, beta), (alpha, beta), (0, 1), (1e-6, 100e-6),
    (alpha, beta), (alpha, beta), (0, 1), (1e-6, 100e-6)
]

print('Fitting')
t1 = time_ns()
res = differential_evolution(cost_function,
                             k_bounds,
                             args=(E_sam, E_ref_w, f_ref),
                             popsize=60,
                             maxiter=2000,
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
