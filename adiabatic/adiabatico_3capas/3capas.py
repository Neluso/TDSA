import numpy as np
import genetic_denoising
from TDSA import *


deg_in = 0  # incidence angle in degrees
snell_sin = n_air * sin(deg_in * pi / 180)
# n_subs = 1.17 - 0.0 * 1j  # substrate refractive index -- cork
n_subs = 1e20 - 0.0 * 1j  # substrate refractive index -- metal
# n_subs = n_air_cplx


# function definitions
def theta(n):
    return arcsin(snell_sin / real(n))


def t_i_j(n_i, n_j):
    return (2 * n_i) / (n_i + n_j)


def r_i_j(n_i, n_j):
    return (n_i - n_j) / (n_i + n_j)


def phase_factor(n_i, thick_i, freq):  # theta in radians
    omg = 2 * pi * freq
    phi = omg * thick_i / c_0
    return exp(1j * n_i * phi)


def fabry_perot(n_i, thick_i, n_im1, n_ip1, freq):
    fp = 1 - r_i_j(n_im1, n_i) * r_i_j(n_i, n_ip1) * phase_factor(n_i, 2 * thick_i, freq)
    return 1 / fp


def P_i(n_i, thick_i, freq):
    phi_fac_1 = phase_factor(n_i, thick_i, freq)
    phi_fac_2 = phase_factor(n_i, - thick_i, freq)
    return array(((phi_fac_1, 0), (0, phi_fac_2)))


def D_i_j(n_i, n_j):
    t_12 = t_i_j(n_i, n_j)
    inv_t_12 = 1 / t_12
    return inv_t_12 * array(((1, r_i_j(n_i, n_j)), (r_i_j(n_i, n_j), 1)))


def H_sim(freq, n_i, n_im1, n_ip1, thick_i, H_prev):  # in first input: t_i_j(n_im1, n_i)
    H_i = H_prev * t_i_j(n_i, n_ip1) * phase_factor(n_i, thick_i, freq) * fabry_perot(n_i, thick_i, n_im1, n_ip1, freq)
    return H_i


def full_H(freq, n_s, thick_s):
    n_1, n_2, n_3, n_4, n_5 = n_s
    thick_1, thick_2, thick_3, thick_4, thick_5 = thick_s
    H_i = H_sim(freq, n_1, n_air, n_2, thick_1, t_i_j(n_air, n_1))
    H_i = H_sim(freq, n_2, n_1, n_3, thick_2, H_i)
    H_i = H_sim(freq, n_3, n_2, n_4, thick_3, H_i)
    H_i = H_sim(freq, n_4, n_3, n_5, thick_4, H_i)
    H_i = H_sim(freq, n_5, n_4, n_air, thick_5, H_i)
    return H_i


t_ref, E_ref = read_1file('./20210628_adiabatic/ref2.txt')
t_sam, E_sam = read_1file('./20210628_adiabatic/sam2.txt')
delta_t_ref = mean(diff(t_ref))
enlargement = 3 * E_ref.size
# ref_pulse_idx = centre_loc(E_ref)
E_ref = zero_padding(E_ref, 0, enlargement)
t_ref = concatenate((t_ref, t_ref[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
E_sam = zero_padding(E_sam, 0, enlargement)
t_sam = concatenate((t_sam, t_sam[-1] * ones(enlargement) + delta_t_ref * arange(1, enlargement + 1)))
t_ref *= 1e-12
t_sam *= 1e-12
# E_ref /= amax(abs(E_ref))
# E_sam /= amax(abs(E_sam))
# plot(E_ref)
# plot(E_sam)
E_sam = smooth(E_sam, span=12)
# plot(E_sam)
f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)
n_1 = 1.55 - 1j * 0.01 * f_ref * 1e-12
n_2 = 1.45 - 1j * 0.01 * f_ref * 1e-12
n_s = [1.55, 1.525, 1.5, 1.475, 1.45]
thick_s = [1000e-6, 330e-6, 330e-6, 330e-6, 1000e-6]


plot(abs(irfft(E_sam_w / E_ref_w * wiener_filter(E_ref_w, beta=1e-4))))
# plot(abs(E_sam_w / E_ref_w))
H_teo = full_H(f_ref, n_s, thick_s)
# plot(abs(H_teo)) # * wiener_filter(E_ref_w, beta=1e-4)))
# plot(E_ref / amax(abs(E_ref)))
# plot(irfft(E_ref_w * full_H(f_ref, n_s, [100e-6, 33e-6, 33e-6, 33e-6, 100e-6])))
# plot(irfft(E_ref_w * full_H(f_ref, n_s, [50e-6, 17e-6, 17e-6, 17e-6, 50e-6])))
# plot(irfft(E_ref_w * full_H(f_ref, n_s, [30e-6, 10e-6, 10e-6, 10e-6, 30e-6])))
show()
