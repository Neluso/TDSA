import numpy as np

import genetic_denoising
from TDSA import *


deg_in = 45  # incidence angle in degrees
snell_sin = n_air * sin(deg_in * pi / 180)
# n_subs = 1.17 - 0.0 * 1j  # substrate refractive index -- cork
n_subs = 1e20 - 0.0 * 1j  # substrate refractive index -- metal
# n_subs = n_air_cplx


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
    return (n_l_1 - n_l) / (n_l_1 + n_l)


def phase_factor(n, thick, freq):  # theta in radians
    omg = 2 * pi * freq
    thick *= cos(theta(n))
    phi = 2 * omg * thick / c_0
    return exp(- 1j * n * phi)


def H_sim(freq, n_l, n_l_1, thick, H_subs):
    H_i = H_subs  # cr_l_1_l(n_subs, n - 1j * k)
    rlm1l = cr_l_1_l(n_l, n_l_1)
    tt = ct2(n_l, n_l_1)
    exp_phi = phase_factor(n_l, thick, freq)
    H_i = rlm1l + (tt * H_i * exp_phi) / (1 + rlm1l * H_i * exp_phi)
    return H_i


t_ref, E_ref = read_1file('ref1.txt')
t_ref_lock, E_ref_lock = read_slow_data('ref1_lockin.txt')
t_sam, E_sam = read_1file('sam1.txt')
t_sam_lock, E_sam_lock = read_slow_data('sam1_lockin.txt')


if t_ref.size < t_sam.size:
    E_sam = interp(t_ref, t_sam, E_sam)
elif t_sam.size < t_ref.size:
    E_ref = interp(t_sam, t_ref, E_ref)

t_ref *= 1e-12
t_ref_lock *= 1e-12
t_sam *= 1e-12
t_sam_lock *= 1e-12

# E_ref_lock = - smooth(E_ref_lock, span=1)


# plot(E_ref_lock)
# plot(E_sam_lock)

# E_ref_lock = genetic_denoising.genetic_deno(t_ref_lock, E_ref_lock, t_ref_lock, E_ref_lock)
# E_sam_lock = genetic_denoising.genetic_deno(t_ref_lock, E_ref_lock, t_sam_lock, E_sam_lock)

# plot(E_ref_lock)
# plot(E_sam_lock)
# show()
# quit()

f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
f_ref_lock, E_ref_lock_w = fourier_analysis(t_ref_lock, E_ref_lock)
f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
f_sam_lock, E_sam_lock_w = fourier_analysis(t_sam_lock, E_sam_lock)


block_size = 200  # points
block_freq = 1.2  # THz
filt = ones(E_ref_w.shape)
block_band = zeros(block_size)
block_freq *= 1e12
block_freq_idx = where(f_ref >= block_freq)[0][0]
filt[block_freq_idx:block_freq_idx + block_size] = block_band


# figure(1)
# plot(irfft(E_ref_w))
# plot(irfft(E_ref_w * filt))
# # figure(2)
# # plot(unwrap(angle(E_ref_w)))
# # plot(unwrap(angle(E_sam_w)))
# show()
# quit()


# if f_ref.size < f_sam.size:
#     E_sam_w = interp(f_ref, f_sam, E_sam_w)
#     f_good = f_ref
# elif t_sam.size < t_ref.size:
#     E_ref_w = interp(f_sam, f_ref, E_ref_w)
#     f_good = f_sam

# samp_int = float(mean(diff(f_good)))
# t_good = rfftfreq(f_good.size, d=samp_int)


H_w = E_sam_w / E_ref_w
H_w_lock = E_sam_lock_w / E_ref_lock_w

# plot(E_ref_lock)
# show()
# quit()
H_w_lock_corr = correlate(E_sam_lock, E_ref_lock[1000:1200], 'same')
plot(H_w_lock_corr)
show()
quit()
# plot(abs(H_w))
# show()
# quit()

# H_w_lock *= wiener_filter(E_ref_lock_w, beta=1e-9)
SNR = abs(E_ref_lock_w) / 0.02
# H_w_lock *= wiener_filter(E_ref_lock_w, beta=1/SNR)


n_1 = 1.5
thick_1 = 30e-6
n_2 = 2.6
thick_2 = 180e-6


H_teo_no_adiab = cr_l_1_l(n_subs, n_2)
H_teo_no_adiab = H_sim(f_ref_lock, n_2, n_1, thick_2, H_teo_no_adiab)
H_teo_no_adiab = H_sim(f_ref_lock, n_1, n_air, thick_1, H_teo_no_adiab)
# H_teo_no_adiab *= phi_air
E_sim_nad = irfft(H_teo_no_adiab * E_ref_lock_w)
H_w_lock_corr_nad = correlate(E_sim_nad, E_ref_lock[1000:1200], 'same')
# t_ref, E_ref = read_1file('./celo_vid/ref1.txt')
# t_sam, E_sam = read_1file('./celo_vid/sam1.txt')
SNR = abs(E_ref_lock_w) / 0.02
G = conjugate(H_w_lock) / (abs(H_w_lock)**2 + 1 / SNR)
E_sam_aux = E_sam
E_sim_aux = smooth(E_sam_aux, span=1000)
plot(E_sam)
plot(E_sim_aux, lw=0.8)
show()
h_t = irfft(H_w)
for i in range(10):
    print(np.dot(h_t, E_sam) / np.dot(h_t, E_sim_aux))
    E_sim_aux = E_sim_aux * (np.dot(h_t, E_sam) / np.dot(h_t, E_sim_aux))
    plot(E_sim_aux, lw=0.8)
# plot(toDb(G * E_sam_lock_w))
# plot(toDb(E_sam_lock_w))
# plot(toDb(G))
# plot(SNR)
# plot(irfft(H_w_lock), lw=0.8)
# plot(irfft(H_w_lock * wiener_filter(E_ref_lock_w)), lw=0.8)
# plot(irfft(H_w_lock * wiener_filter(E_ref_lock_w, beta=1/SNR)), lw=0.8)
# figure(1)
# plot(toDb(E_ref_lock_w))
# plot(toDb(E_sam_lock_w))
# plot(E_sam_lock / amax(abs(E_sam_lock)))
# plot(E_sim_nad/ amax(abs(E_sim_nad)))
# figure(2)
# plot(irfft(H_teo_no_adiab))
# plot(irfft(H_w_lock))  # / amax(abs(H_w_lock))))
# plot(H_w_lock_corr)
# plot(correlate(E_sim_nad, E_ref_lock, 'full'))
show()
