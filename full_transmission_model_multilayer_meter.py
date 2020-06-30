from TDSA import *
from scipy.optimize import differential_evolution
from scipy.signal.windows import tukey
from time import time_ns, strftime, gmtime
from genetic_denoising import genetic_deno

# constant definitions
deg_in = 30  # incidence angle in degrees
snell_sin = n_air * sin(deg_in * pi / 180)
layers = 3
# epsilon_model = 'debye'
# epsilon_model = 'cole'
# n_subs = 1.17 - 0.0 * 1j  # substrate refractive index -- cork
# n_subs = 1.17e20 - 0.0 * 1j  # substrate refractive index -- metal
# n_subs = 1.25 - 0.0 * 1j  # substrate refractive index -- cork 2.0


# function definitions
def theta(n):
    return arcsin(snell_sin / real(n))


def ct(n_1, n_2):
    return 2 * n_1 / (n_1 + n_2)


def cr(n_1, n_2):
    return (n_2 - n_1) / (n_1 + n_2)


def ct2(n_l, n_l_1):
    n_l *= cos(theta(n_l))
    n_l_1 *= cos(theta(n_l_1))
    return 4 * n_l * n_l_1 / (n_l + n_l_1)**2


def cr_l_1_l(n_l, n_l_1):  # from n_l-1 to n_l
    n_l_1 *= cos(theta(n_l_1))
    n_l *= cos(theta(n_l))
    return - (n_l_1 - n_l) / (n_l_1 + n_l)


def phase_factor(n, thick, freq):  # theta in radians
    omg = 2 * pi * freq
    thick *= cos(theta(n))
    phi = omg * thick / c_0
    return exp(- 1j * phi)


def epsilon(e_s, e_inf, tau, freq):  # Debye model
    omg = 2 * pi * freq
    e_w = e_inf + (e_s - e_inf) / (1 + 1j * omg * tau)
    return e_w


def nk_from_eps(e_s, e_inf, tau, freq):
    e_w = epsilon(e_s, e_inf, tau, freq)
    n = sqrt((abs(e_w) + real(e_w)) / 2)
    k = sqrt((abs(e_w) - real(e_w)) / 2)
    return n, k


def fabry_perot(n_in, n_out, n_l, d_l, freq):
    return 1 / (1 - cr(n_in, n_l) * cr(n_l, n_out) * phase_factor(n_l, 2*d_l, freq))


def H_sim(freq, params):  # d_air, d_subs, n_1, k_1, d_1, n_2, k_2, d_2, n_3, k_3, d_3):
    # d_air, d_subs, e_inf_1, e_s_1, tau_1, d_1, e_inf_2, e_s_2, tau_2, d_2, e_inf_3, e_s_3, tau_3, d_3 = params
    d_air, d_subs, e_inf_1, e_s_1, tau_1, d_1, e_inf_2, e_s_2, tau_2, d_2, e_inf_3, e_s_3, tau_3, d_3, n_subs = params

    n_1, k_1 = nk_from_eps(e_inf_1, e_s_1, tau_1, freq)  # debye model
    n_2, k_2 = nk_from_eps(e_inf_2, e_s_2, tau_2, freq)  # debye model
    n_3, k_3 = nk_from_eps(e_inf_3, e_s_3, tau_3, freq)  # debye model
    
    # k_1 *= freq * 1e-12
    # k_2 *= freq * 1e-12
    # k_3 *= freq * 1e-12
    
    # d_air1 = - d_subs - d_1 - d_2 - d_3
    
    H_0 = phase_factor(n_air_cplx, d_air, freq)
    
    # Layer 1
    H_1 = ct(n_air_cplx, n_1 - 1j * k_1) * ct(n_1 - 1j * k_1, n_2 - 1j * k_2) * phase_factor(n_1 - 1j * k_1, d_1, freq)
    H_1 *= fabry_perot(n_air_cplx, n_2 - 1j * k_2, n_1 - 1j * k_1, d_1, freq)

    # Layer 2
    # H_2 = ct(n_1 - 1j * k_1, n_2 - 1j * k_2) * ct(n_2 - 1j * k_2, n_3 - 1j * k_3) * phase_factor(n_2 - 1j * k_2, d_2, freq)
    H_2 = H_1 * ct(n_2 - 1j * k_2, n_3 - 1j * k_3) * phase_factor(n_2 - 1j * k_2, d_2, freq)
    H_2 *= fabry_perot(n_1 - 1j * k_1, n_3 - 1j * k_3, n_2 - 1j * k_2, d_2, freq)
    
    # Layer 3
    # H_3 = ct(n_2 - 1j * k_2, n_3 - 1j * k_3) * ct(n_3 - 1j * k_3, n_subs) * phase_factor(n_3 - 1j * k_3, d_3, freq)
    H_3 = H_2 * ct(n_3 - 1j * k_3, n_subs) * phase_factor(n_3 - 1j * k_3, d_3, freq)
    H_3 *= fabry_perot(n_2 - 1j * k_2, n_subs, n_3 - 1j * k_3, d_3, freq)
    
    # Substrate
    # H_subs = ct(n_3 - 1j * k_3, n_subs) * ct(n_subs, n_air_cplx) * phase_factor(n_subs, d_subs, freq)
    H_subs = H_3 * ct(n_subs, n_air_cplx) * phase_factor(n_subs, d_subs, freq)
    # H_subs *= fabry_perot(n_3 - 1j * k_3, n_air_cplx, n_subs, d_subs, freq)
    
    # return H_1 * H_2 * H_3 * H_subs
    return H_0 * H_1 * H_2 * H_3 * H_subs


def cost_function(params, *args):
    E_sam, E_ref_w, freqs = args
    H_teo = H_sim(freqs, params)  # d_air, d_subs, n_1, k_1, d_1, n_2, k_2, d_2, n_3, k_3, d_3)
    E_sam_teo_w = E_ref_w * H_teo
    E_sam_teo = irfft(E_sam_teo_w, n=E_sam.size)
    return sum((E_sam_teo - E_sam) ** 2)  # + sum(params**2)


# Main script
# Boleto 176054
t_ref, E_ref = read_1file('./data/airbus_transmision/ref_176054.txt')
t_sam, E_sam = read_1file('./data/airbus_transmision/sam1_176054.txt')

# Boleto 180881
# t_ref, E_ref = read_1file('./data/airbus_transmision/ref_180881.txt')
# t_sam, E_sam = read_1file('./data/airbus_transmision/sam3_180881.txt')

# Boleto 177910
# t_ref, E_ref = read_1file('./data/airbus_transmision/ref_177910.txt')
# t_sam, E_sam = read_1file('./data/airbus_transmision/sam1_177910.txt')


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
E_sam_raw = E_sam

t_ref *= 1e-12
t_sam *= 1e-12

# E_sam = genetic_deno(t_ref, E_ref, t_sam, E_sam_raw)

f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)

H_w = E_sam_w / E_ref_w
alpha = 3
beta = 15
k_bounds = [  # calibration
    (-3e-3, 0),  # 3e-3),  # air thickness
    (2.5e-3, 3.5e-3),  # substrate thickness
    (alpha, beta), (alpha, beta), (1e-14, 1e-11), (20e-6, 100e-6),
    (alpha, beta), (alpha, beta), (1e-14, 1e-11), (20e-6, 100e-6),
    (alpha, beta), (alpha, beta), (1e-14, 1e-11), (1e-6, 20e-6)
    , (1, 1.5)  # substrate index
]

print('Fitting')
t1 = time_ns()
res = differential_evolution(cost_function,
                             k_bounds,
                             args=(E_sam, E_ref_w, f_ref),
                             # popsize=150,
                             maxiter=2000,
                             disp=True,  # step cost_function value
                             polish=True
                             )
t2 = time_ns()
print()
d_air = res.x[0]
d_subs = res.x[1]
print(res)
n, k = nk_from_eps(res.x[2], res.x[3], res.x[4], f_ref)
figure(4)
plot(f_ref, n, label='w')
figure(5)
plot(f_ref, k, label='w')
print('White coat -', 'n:', round(mean(n), 2), 'k:', round(mean(k), 2), 'd:', round(res.x[5] * 1e6, 2), 'um')
print('\t\t e_s:', res.x[3], 'e_inf', res.x[2], 'tau', res.x[4])
n, k = nk_from_eps(res.x[6], res.x[7], res.x[8], f_ref)
figure(4)
plot(f_ref, n, label='g')
figure(5)
plot(f_ref, k, label='g')
print('Green coat -', 'n:', round(mean(n), 2), 'k:', round(mean(k), 2), 'd:', round(res.x[9] * 1e6, 2), 'um')
print('\t\t e_s:', res.x[7], 'e_inf', res.x[6], 'tau', res.x[8])
n, k = nk_from_eps(res.x[10], res.x[11], res.x[12], f_ref)
figure(4)
plot(f_ref, n, label='p')
figure(5)
plot(f_ref, k, label='p')
print('Primer coat -', 'n:', round(mean(n), 2), 'k:', round(mean(k), 2), 'd:', round(res.x[13] * 1e6, 2), 'um')
print('\t\t e_s:', res.x[11], 'e_inf', res.x[10], 'tau', res.x[12])
print('Total:', round((res.x[5] + res.x[9] + res.x[13]) * 1e6, 2), 'um')
print()
print('Air - d: ', round(d_air * 1e3, 2), 'mm')
print('Subs - n:', res.x[-1], 'd: ', round(d_subs * 1e3, 2), 'mm')
secs = (t2-t1)*1e-9
if secs < 3600:
    print('Processing time (mm:ss):', strftime('%M:%S', gmtime(secs)))
else:
    print('Processing time (mm:ss):', strftime('%H:%M:%S', gmtime(secs)))
print()


H_fit = H_sim(f_ref, res.x)
E_fit = irfft(H_fit * E_ref_w)

figure(1)
# plot(t_ref, E_ref, label='ref')
plot(t_sam, E_sam_raw, label='raw')
plot(t_sam, E_sam, label='sam')
plot(t_sam, E_fit, label='fit')
legend()

figure(2)
plot(f_ref, abs(H_w))
plot(f_ref, abs(H_fit))

figure(3)
plot(f_ref, unwrap(angle(H_w)))
plot(f_ref, unwrap(angle(H_fit)))

figure(4)
title('n')
xlim([0, 2e12])
legend()
figure(5)
title('k')
xlim([0, 2e12])
legend()

show()
