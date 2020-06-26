from TDSA import *
from scipy.optimize import curve_fit
from scipy.signal.windows import tukey
from time import time_ns
from genetic_denoising import genetic_deno

# constant definitions
n_subs = 1.17 - 0.0 * 1j  # substrate refractive index -- cork
# n_subs = 1.17e20 - 0.0 * 1j  # substrate refractive index -- metal
# n_subs = 1.25 - 0.0 * 1j  # substrate refractive index -- cork 2.0


# debug variables
debug_value_1 = list()
debug_value_2 = list()
debug_value_3 = list()


# function definitions
def ct(n_1, n_2):
    return 2 * n_1 / (n_1 + n_2)


def cr(n_1, n_2):
    return (n_2 - n_1) / (n_1 + n_2)


def phase_factor(n, thick, freq):  # theta in radians
    omg = 2 * pi * freq
    phi = 2 * omg * thick / c_0
    return np.exp(- 1j * n * phi)


def fabry_perot(n_in, n_out, n_l, d_l, freq):
    return 1 / (1 - cr(n_in, n_l) * cr(n_l, n_out) * phase_factor(n_l, d_l, freq))


def H_sim(freq, d_air, d_subs, n_1, k_1, d_1, n_2, k_2, d_2, n_3, k_3, d_3):
    k_1 *= freq * 1e-12
    k_2 *= freq * 1e-12
    k_3 *= freq * 1e-12

    H_0 = phase_factor(n_air_cplx, d_air, freq)

    # Layer 1
    H_1 = ct(n_air_cplx, n_1 - 1j * k_1) * ct(n_1 - 1j * k_1, n_2 - 1j * k_2) * phase_factor(n_1 - 1j * k_1, d_1, freq)
    H_1 *= fabry_perot(n_air_cplx, n_2 - 1j * k_2, n_1 - 1j * k_1, d_1, freq)

    #Layer 2
    H_2 = ct(n_1 - 1j * k_1, n_2 - 1j * k_2) * ct(n_2 - 1j * k_2, n_3 - 1j * k_3) * phase_factor(n_2 - 1j * k_2, d_2, freq)
    H_2 *= fabry_perot(n_1 - 1j * k_1, n_3 - 1j * k_3, n_2 - 1j * k_2, d_2, freq)

    # Layer 3
    H_3 = ct(n_2 - 1j * k_2, n_3 - 1j * k_3) * ct(n_3 - 1j * k_3, n_subs) * phase_factor(n_3 - 1j * k_3, d_3, freq)
    H_3 *= fabry_perot(n_2 - 1j * k_2, n_subs, n_3 - 1j * k_3, d_3, freq)

    # Substrate
    H_subs = ct(n_3 - 1j * k_3, n_subs) * ct(n_subs, n_air_cplx) * phase_factor(n_subs, d_subs, freq)
    # H_subs *= fabry_perot(n_3 - 1j * k_3, n_air_cplx, n_subs, d_subs, freq)

    return H_0 * H_1 * H_2 * H_3 * H_subs


# Main script
# Main script
# Boleto 176054
t_ref, E_ref = read_1file('./data/airbus_transmision/ref_176054.txt')
t_sam, E_sam = read_1file('./data/airbus_transmision/sam1_176054.txt')

# Boleto 180881
# t_ref, E_ref = read_1file('./data/airbus_transmision/ref_180881.txt')
# t_sam, E_sam = read_1file('./data/airbus_transmision/sam3_180881.txt')

# Boleto 177910
# t_ref, E_ref = read_1file('./data/airbus_transmision/ref_177910.txt')
# t_sam, E_sam = read_1file('./data/airbus_transmision/sam3_177910.txt')


t_ref *= 1e-12
t_sam *= 1e-12

# plot(t_sam, E_sam, label='org')
# print('Denoising')
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


p0 = array((20e-6, 3e-3, 2.0, 0.1, 60e-6, 2.0, 0.1, 30e-6, 2.0, 0.1, 10e-6))
lwr_bnds = array((0, 0, 1e-3, 0, 1e-6, 1, 0, 1e-6, 1, 0, 1e-6))
hgr_bnds = array((100e-6, 5e-3, 5, 10, 100e-6, 5, 10, 100e-6, 5, 10, 100e-6))

# p0 = array((20e-6, 3e-3, 3, 0.1, 80e-6, 3, 0.1, 19e-6, 3, 0.1, 4e-6))
# lwr_bnds = array((0, 1e-3, 2, 0, 80e-6, 2, 0, 19e-6, 2, 0, 4e-6))
# hgr_bnds = array((100e-6, 5e-3, 5, 1, 81e-6, 5, 1, 20e-6, 5, 1, 5e-6))

# print(H_sim(f_ref, 20e-6, 2.0, 0.1, 50e-6, 2.0, 0.1, 50e-6, 2.0, 0.1, 50e-6))
# quit()

print('Fitting')


def E_sim(times, d_air, d_subs, n_1, k_1, d_1, n_2, k_2, d_2, n_3, k_3, d_3):
    H_0 = H_sim(f_ref, d_air, d_subs, n_1, k_1, d_1, n_2, k_2, d_2, n_3, k_3, d_3)
    E_0 = irfft(H_0 * E_ref_w)
    return E_0


t1 = time_ns()
popt, pcov = curve_fit(E_sim, t_sam, E_sam, p0,
                       bounds=(lwr_bnds, hgr_bnds),
                       method='dogbox'
                       )
t2 = time_ns()
# print('popt =', popt)
# print('pcov =', pcov)

print('Results:')
print('White coat --- n:', round(popt[2], 2), 'k:', round(popt[3], 2), 'd:', round(popt[4] * 1e6, 0))
print('Green coat --- n:', round(popt[5], 2), 'k:', round(popt[6], 2), 'd:', round(popt[7] * 1e6, 0))
print('Primer --- n:', round(popt[8], 2), 'k:', round(popt[9], 2), 'd:', round(popt[10] * 1e6, 0))
print('Cork --- d:', round(popt[1], 2) * 1e3)

figure(1)
plot(t_ref, E_ref, label='ref')
plot(t_sam, E_sam, label='sam') 
plot(t_sam, E_sim(t_sam, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10]), label='fit')
legend()

figure(2)
plot(f_ref, abs(H_w))
plot(f_ref, abs(H_sim(f_ref, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10])))

figure(3)
plot(f_ref, unwrap(angle(H_w)))
plot(f_ref, unwrap(angle(H_sim(f_ref, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10]))))
show()
