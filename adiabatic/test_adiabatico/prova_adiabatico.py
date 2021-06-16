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


t_ref, E_ref = read_1file('ref3.txt')
t_sam, E_sam = read_1file('sam3.txt')

E_ref = smooth(E_ref, span=1)
E_sam = smooth(E_sam, span=1)

if t_ref.size < t_sam.size:
    E_sam = interp(t_ref, t_sam, E_sam)
elif t_sam.size < t_ref.size:
    E_ref = interp(t_sam, t_ref, E_ref)

t_ref *= 1e-12
t_sam *= 1e-12


f_ref, E_ref_w = fourier_analysis(t_ref, E_ref)
f_sam, E_sam_w = fourier_analysis(t_sam, E_sam)


# if f_ref.size < f_sam.size:
#     E_sam_w = interp(f_ref, f_sam, E_sam_w)
#     f_good = f_ref
# elif t_sam.size < t_ref.size:
#     E_ref_w = interp(f_sam, f_ref, E_ref_w)
#     f_good = f_sam

# samp_int = float(mean(diff(f_good)))
# t_good = rfftfreq(f_good.size, d=samp_int)


H_w = E_sam_w / E_ref_w


# n_1 = 1.55
# thick_1 = 750e-6
# n_2 = 1.65
# thick_2 = 750e-6


# H_teo_no_adiab = cr_l_1_l(n_subs, n_2)
# H_teo_no_adiab = H_sim(f_ref, n_2, n_1, thick_2, H_teo_no_adiab)
# H_teo_no_adiab = H_sim(f_ref, n_1, n_air, thick_1, H_teo_no_adiab)
# # H_teo_no_adiab *= phi_air
# E_sim_nad = irfft(H_teo_no_adiab * E_ref_w)
# figure(1)
plot(E_sam)
# figure(2)
plot(irfft(H_w))
# plot(irfft(H_teo_no_adiab))
show()
